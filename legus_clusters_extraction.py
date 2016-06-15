#!/usr/bin/env python
"""
Code used to extract a candidate list of clusters in LEGUS images
based on a magnitude cut and a concentration index cut.
Read the tutorial [legus_clusters_extraction.pdf]
to understand the installation,
dependencies, and use of this code.
@author: Leonardo UBEDA
@organization: Space Telescope Science Institute
@team: Legacy ExtraGalactic UV Survey


@change: 18 Jan 2016 Jenna Ryon
            In Step 3: Photometry of sources located outside of FOV are given 66.666 (instead of 99.999)
            In Step 5: Remove duplicates from catalog
            In Steps 5 & 6: Change photometry of sources with undefined CI-based aperture corrections to 44.444


@major updates: 23June2015 Jenna Ryon, Dave Cook, Angela Adamo
            Significant overhaul of the structure and output


@change: 17Dec2014 DaveC
		   added minor axis ticks
                   changed /usr/bin/sed to sed
                   added a 'userannulus' to be +3 pixels from useraperture
                   	dwarfs require a larger aperture
                   find correct aperture correction for Linux users
                   	glob gave an incorrect order for filters
                   Add centroid function to step 6. The manual process
                   	does not neccessarily center on the object.

@change: 26NOV2014 Leonardo version 3.2
                   added step 6 to perform photometry on additional clusters.


@change: 01NOV2014 David version 3.1
                   find the correct pix.mag files.
                   glob gave either 4px or 20px as the first element depending on the filter


@change: 15OCT2014 Leonardo, Angela, Jenna Linda version 3.0
                            included filter F606W, and made multiple changes
                            suggested by Angela, Jenna, and Linda
                            email october 14
@change: 14OCT2014 Leonardo, Jenna version 2.2
                            add change suggested by Jenna
                            email october 13
@change: 01OCT2014 Leonardo version 2.1
                            use Angela's coding for photometry
@change: 10SEP2014 Leonardo version 2.0
                            translate Angela's code into Python
@change: 06AUG2014 Angela version 1.0

"""


# Import packages

import numpy as np
import time, shutil
import sys
import string
import matplotlib.pyplot as plt
import os, glob, pyfits, pdb, math
from pyraf import iraf
import pywcs



# Location of target directory (should be current directory)
target_dir = os.getcwd()

# Location of /data directory and this code
pydir = '/'.join(target_dir.split('/')[:-1])

# Open necessary packages in iraf
iraf.noao(_doprint=0)
iraf.digiphot(_doprint=0)
iraf.obsutil(_doprint=0)
iraf.daophot(_doprint=0)

# Define input file
infile = 'legus_clusters_extraction.input'

# Verify that input file exists and load it in
if os.path.exists(target_dir + '/' + infile) == False:
    print ''
    print 'File ', infile, ' could not be found in ', target_dir
    sys.exit('quitting now...')

input = np.loadtxt(target_dir + '/' + infile, unpack=True, skiprows=0, delimiter='#', dtype='str')

# Read in and print contents of input file
print ''
print '________________________________________________________________________________________'
print ''

print 'target name : ', input[0]
target = input[0].strip()

print 'target distance : ', input[1]
distance = input[1].strip()

print 'step 1 : ', input[2].strip()
flagstep1 = input[2].strip()

if flagstep1 not in ('yes', 'no', 'YES', 'NO'):
    sys.exit('Wrong input in line 3')

print 'sextractor image : ', input[3].strip()
seximage = input[3].strip()

print 'step 2 : ', input[4].strip()
flagstep2 = input[4].strip()

if flagstep2 not in ('yes', 'no', 'YES', 'NO'):
    sys.exit('Wrong input in line 4')

print 'aperture radius : ', input[5].strip()
useraperture = input[5].strip()


if useraperture not in ['4.0','5.0','6.0']:
    print 'Chosen aperture radius is not supported.'
    print 'Please choose 4.0, 5.0, or 6.0 pixels.'
    sys.exit('Quitting...')

print 'step 3 : ', input[6].strip()
flagstep3 = input[6].strip()

if flagstep3 not in ('yes', 'no', 'YES', 'NO'):
    sys.exit('Wrong input in line 6')

print 'CI : ', input[7].strip()
ci = input[7].strip()

print 'step 4 : ', input[8].strip()
flagstep4 = input[8].strip()

if flagstep4 not in ('yes', 'no', 'YES', 'NO'):
    sys.exit('Wrong input in line 7')

print 'upper limit for avg apcorr : ', input[9].strip()
uplim_apcor = float(input[9].strip())

print 'lower limit for avg apcorr : ', input[10].strip()
lolim_apcor = float(input[10].strip())

print 'step 5 : ',  input[11].strip()
flagstep5 = input[11].strip()

if flagstep5 not in ('yes', 'no', 'YES', 'NO'):
    sys.exit('Wrong input in line 12')

print 'step 6 : ', input[12].strip()
flagstep6 = input[12].strip()

if flagstep6 not in ('yes', 'no', 'YES', 'NO'):
    sys.exit('Wrong input in line 11')

print 'additional sources coordinates file: ', input[13].strip()
inputlist = input[13].strip()

print 'smoothed image for cluster centering: ', input[14].strip()
smooth_image = input[14].strip()

print ''
print '________________________________________________________________________________________'


#####  Setup a few things  #####


# Verify that the target has an entry in the galactic extinction table
extinct = pydir + '/data/legus_galactic_extinction.tab'
target_ext = np.loadtxt(extinct, unpack=True, skiprows=2, dtype='str', usecols=(0,))

id = (target_ext == target)

if ~np.any(id):
    print 'Please update table ' + extinct + ' to include LEGUS target ' + target
    print ''
    sys.exit('Quitting now...')


# Verify that we are executing only one step
steplist = [flagstep1, flagstep2, flagstep3, flagstep4, flagstep5, flagstep6]

if (steplist.count('no') + steplist.count('NO')) == 6:
    print 'Nothing to do.'
    sys.exit('Quitting now...')

if (steplist.count('yes') + steplist.count('YES')) >= 2:
    print 'Please select only one step to execute. Check your input file.'
    sys.exit('Quitting now...')

# Set science sky annulus and width here so that they can be changed
# without hunting for them
annulus_sci = 7.0
dannulus_sci = 1.0


# Set some IRAF phot parameters that will be common to all photometry runs
iraf.unlearn('datapars')
iraf.datapars.scale = 1.0
iraf.datapars.fwhmpsf = 2.0
iraf.datapars.sigma = 0.01
iraf.datapars.readnoise = 5.0
iraf.datapars.itime = 1.

iraf.unlearn('centerpars')
iraf.centerpars.cbox = 1
iraf.centerpars.cmaxiter = 3
iraf.centerpars.maxshift = 1

iraf.unlearn('fitskypars')
iraf.fitskypars.salgori = 'mode'

iraf.unlearn('photpars')

iraf.unlearn('phot')
iraf.phot.interactive = 'no'
iraf.phot.verbose = 'yes'
iraf.phot.verify = 'no'




##############   STEP 1   ##############

#####  Create directories and move files  #####

if  flagstep1 == 'yes':

    print 'Creating directories ...'

    if os.path.exists(target_dir + '/s_extraction') == False:
        os.makedirs(target_dir + '/s_extraction')

    if os.path.exists(target_dir + '/photometry') == False:
        os.makedirs(target_dir + '/photometry')

    if os.path.exists(target_dir + '/plots') == False:
        os.makedirs(target_dir + '/plots')

    if os.path.exists(target_dir + '/img') == False:
        os.makedirs(target_dir + '/img')

    if os.path.exists(target_dir + '/init') == False:
        print 'No init/ directory. Please create.'
        sys.exit()

    # Only move images to /img directory if they aren't there already
    im_check = glob.glob(target_dir + '/img/*.fits')

    if len(im_check) == 0:
        print 'Moving fits files ...'

        imlist = glob.glob(target_dir + '/*.fits')

        for impath in imlist:
            source = impath
            image = impath.split('/')[-1]
            destination = target_dir + '/img/' + image
            os.rename(source, destination)

            # Update header of each image slightly
            pf = pyfits.open(destination, mode='update')
            pf[0].header['BUNIT']= 'ELECTRONS/S'
            pf.close()


    source =  target_dir + '/init/output.param'
    destination = target_dir + '/s_extraction/output.param'
    shutil.copyfile(source, destination)

    source =  target_dir + '/init/R2_wl_aa.config'
    destination = target_dir + '/s_extraction/R2_wl_aa.config'
    shutil.copyfile(source, destination)

    source =  target_dir + '/init/default.nnw'
    destination = target_dir + '/s_extraction/default.nnw'
    shutil.copyfile(source, destination)


#####  Source extraction  #####

    print 'Executing SExtractor on user selected image : ', seximage

    # Verify that file exists
    if os.path.exists(target_dir + '/img/' + seximage) == False:
        print 'File ' + seximage + ' could not be found in ' + target_dir + '/img/'
        sys.exit('Quitting now...')

    # Run sextractor
    os.chdir(target_dir + '/s_extraction')
    command = 'sex ' + target_dir + '/img/' + seximage + '  -c R2_wl_aa.config'
    os.system(command)
    os.chdir(target_dir)

    # Read in results and make regions file of objects sextracted
    xx, yy = np.loadtxt(target_dir + '/s_extraction/R2_wl_dpop_detarea.cat', unpack=True,
                        skiprows=5, usecols=(0,1))

    outputfile = target_dir + '/s_extraction/catalog_ds9_sextractor.reg'

    file = open(outputfile, 'w')
    file.write('global color=blue width=5 font="helvetica 15 normal roman" highlite=1 \n')
    file.write('image\n')

    for i in range(len(xx)):
        newline = 'circle(' + str(xx[i]) + ',' + str(yy[i]) +  ',7) \n'
        file.write(newline)
    file.close()


#####  Pick reference filter  #####

    # Find reference image, f555w or f606w
    image1 = glob.glob(target_dir + '/img/*f555w_sci.fits')
    image2 = glob.glob(target_dir + '/img/*f606w_sci.fits')

    if len(image1) + len(image2) == 0:
        sys.exit('reference image is missing')

    if len(image1) == 1:
        ref_filter = 'f555w'
        image = image1

    if len(image2) == 1:
        ref_filter = 'f606w'
        image = image2

    # Save filter of reference image in a text file for future use
    f = open(target_dir + '/s_extraction/ref_filter.txt', 'w')
    f.write(ref_filter)
    f.close()

    print 'reference filter is : ' + ref_filter

    # Create empty files for user to populate with isolated stars and clusters
    f = open(target_dir + '/photometry/isolated_stars.coo','a')
    f.close()

    f = open(target_dir + '/photometry/isolated_clusters.coo','a')
    f.close()

    print ''
    print 'Check catalog_ds9_sextractor.reg in the /s_extraction directory for'
    print 'the quality of source extraction.'
    print ''
    print 'Populate the following files with isolated stars and clusters: '
    print 'photometry/isolated_stars.coo  --->  xy coordinates of isolated stars'
    print 'photometry/isolated_clusters.coo  --->  xy coordinates of isolated clusters'
    print ''

    sys.exit('STEP 1 is finished')



################   STEP 2   ################


#####  CI and growth curve photometry of isolated sources  #####

if  flagstep2 == 'yes':

    # Preparation for running growth curve photometry on reference filter
    ref_filter = np.loadtxt(target_dir + '/s_extraction/ref_filter.txt', unpack=True,  dtype='str')
    # Remove string from array
    ref_filter = ref_filter[()]

    inst_zp, filter_zp, zp_zp = np.loadtxt(pydir + '/data/legus_zeropoints.tab', unpack=True, dtype='str')

    imlist = glob.glob(target_dir + '/img/*_sci.fits')

    # Finds the full path of reference image from the filter
    ref_image = [image for image in imlist if ref_filter in image][0]

    # Set the necessary variables for photometry on the reference image
    ref_exptime = pyfits.getheader(ref_image)['EXPTIME']
    ref_inst = pyfits.getheader(ref_image)['INSTRUME']
    ref_inst = ref_inst.lower()

    match = (inst_zp == ref_inst) & (filter_zp == ref_filter)
    ref_zp = zp_zp[match]

    # ref_zp is a string within an array, so need to turn into a float
    ref_zp = float(ref_zp[0])

    # Define coordinate and output photometry file paths
    stars_coo = target_dir + '/photometry/isolated_stars.coo'
    stars_mag = target_dir + '/photometry/isolated_stars.mag'

    clusters_coo = target_dir + '/photometry/isolated_clusters.coo'
    clusters_mag = target_dir + '/photometry/isolated_clusters.mag'

    # Remove phot files if they already exist
    if os.path.exists(stars_mag) == True:
        os.remove(stars_mag)

    if os.path.exists(clusters_mag) == True:
        os.remove(clusters_mag)


    # Set iraf parameters for reference filter
    iraf.datapars.epadu = ref_exptime

    iraf.centerpars.calgorithm = 'centroid'

    iraf.fitskypars.annulus = 21.
    iraf.fitskypars.dannulu = 1.

    iraf.photpars.apertures = '1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0'
    iraf.photpars.zmag = ref_zp

    # For isolated stars
    iraf.phot(ref_image, stars_coo, stars_mag)

    # For isolated clusters
    iraf.phot(ref_image, clusters_coo, clusters_mag)

    # Move phot results to new files, remove INDEFs
    stars_growth = target_dir + '/photometry/isolated_stars_growth.mag'
    clusters_growth = target_dir + '/photometry/isolated_clusters_growth.mag'

    cmd = 'grep "*" ' + stars_mag + ' > ' + stars_growth
    os.system(cmd)

    cmd = 'sed -i.bak "s/INDEF/99.999/g" ' + stars_growth
    os.system(cmd)

    cmd = 'grep "*" ' + clusters_mag + ' > ' + clusters_growth
    os.system(cmd)

    cmd = 'sed -i.bak "s/INDEF/99.999/g" ' + clusters_growth
    os.system(cmd)

    # Remove .bak files to prevent confusion
    bak_stars = stars_growth + '.bak'
    bak_clusters = clusters_growth + '.bak'
    os.remove(bak_stars)
    os.remove(bak_clusters)



#####  Histogram of CI values for isolated stars and clusters  #####


    aper_st, mag_st = np.loadtxt(stars_growth, unpack=True, usecols=(0,4))

    # Pick out photometry from apertures of 1 and 3 pixels
    # Stars
    mag_st_1 = mag_st[aper_st == 1.0]
    mag_st_3 = mag_st[aper_st == 3.0]
    ci_st = mag_st_1 - mag_st_3

    # Clusters
    aper_cl, mag_cl = np.loadtxt(clusters_growth, unpack=True, usecols=(0,4))

    mag_cl_1 = mag_cl[aper_cl == 1.0]
    mag_cl_3 = mag_cl[aper_cl == 3.0]
    ci_cl = mag_cl_1 - mag_cl_3

    # Create histogram, figure size in inches
    fig = plt.figure(figsize = (7,7))

    # Draw histogram of CI values
    plt.hist(ci_st, bins = 24, color='red', lw=2., histtype='step', label='stars')
    plt.hist(ci_cl, bins = 24, color='blue', lw=2., histtype='step', label='clusters')

    # Define axes
    plt.xlabel('CI = mag(1px) - mag(3px)')
    plt.ylabel('N')
    plt.minorticks_on()
    plt.legend()

    # Save plot as a PDF file
    fig.savefig(target_dir + '/plots/plot_histogram_ci.pdf')



#####  Growth curves of isolated stars and clusters  #####

    # Stars
    aper_st, flux_st = np.loadtxt(stars_growth, unpack=True, usecols=(0,3))

    ratio_st = np.empty(len(aper_st))

    naper = 20
    nstar = len(aper_st)/naper
    aper_ind = naper - 1

    for k in range(nstar):

        for i in range(naper):

            ratio_st[i + k*naper] = flux_st[i + k*naper]/flux_st[aper_ind + k*naper]

    # Clusters
    aper_cl, flux_cl = np.loadtxt(clusters_growth, unpack=True, usecols=(0,3))

    ratio_cl = np.empty(len(aper_cl))

    naper = 20
    nclus = len(aper_cl)/naper
    aper_ind = naper - 1

    for k in range(nclus):

        for i in range(naper):

            ratio_cl[i + k*naper] = flux_cl[i + k*naper]/flux_cl[aper_ind + k*naper]


    # Find median ratio at each aperture between all the stars and all the clusters
    med_st = np.empty(naper)

    for i in range(naper):

        med_st[i] = np.median(ratio_st[i::naper])


    med_cl = np.empty(naper)

    for i in range(naper):

        med_cl[i] = np.median(ratio_cl[i::naper])


    # Plot growth curves
    fig = plt.figure(figsize = (7,7))

    aper_x = np.arange(naper) + 1

    for i in range(nstar):

        ratio_y = ratio_st[i*naper:(i + 1)*naper]
        plt.plot(aper_x, ratio_y, 'y-')
        plt.annotate(str(i + 1), xy=(8.0, ratio_y[7]),
            horizontalalignment='left', verticalalignment='top', fontsize=6)

    for i in range(nclus):

        ratio_y = ratio_cl[i*naper:(i + 1)*naper]
        plt.plot(aper_x, ratio_y, 'c-')
        plt.annotate(str(i + 1), xy=(12.0, ratio_y[11]),
            horizontalalignment='left', verticalalignment='top', fontsize=6)

    plt.plot(aper_x, med_st, 'r-' , linewidth=4.0)
    plt.plot(aper_x, med_cl, 'b-' , linewidth=4.0)
    plt.hlines(0.5, 0, 20, color='black', linewidth=2, zorder=10)
    plt.vlines(4, 0, 1.1, color='black', linewidth=2, linestyle='dashed', zorder=10)
    plt.vlines(5, 0, 1.1, color='black', linewidth=2, linestyle='dashed', zorder=10)
    plt.vlines(6, 0, 1.1, color='black', linewidth=2, linestyle='dashed', zorder=10)

    plt.ylabel('Normalized Flux ' + ref_filter.upper())
    plt.xlabel('Radius (pix)')
    plt.xlim(1,20)
    plt.minorticks_on()

    fig.savefig(target_dir + '/plots/plot_growth_curve.pdf')



#####  CI photometry for full Sextractor catalog in reference frame  #####

    # Define coordinate and output photometry file paths
    fullcat_coo = target_dir + '/s_extraction/R2_wl_dpop_detarea.cat'
    fullcat_mag = target_dir + '/photometry/ci_fullcat_ref.mag'

    # Remove output file if it already exists
    if os.path.exists(fullcat_mag) == True:
        os.remove(fullcat_mag)

    # Run photometry
    iraf.datapars.epadu = ref_exptime

    iraf.centerpars.calgorithm = 'centroid'

    iraf.fitskypars.annulus = annulus_sci
    iraf.fitskypars.dannulu = dannulus_sci

    iraf.photpars.apertures = '1.0,3.0'
    iraf.photpars.zmag = ref_zp

    iraf.phot(ref_image, fullcat_coo, fullcat_mag)

    # Move phot results to new files, remove INDEFs
    fullcat_mag_short = target_dir + '/photometry/ci_fullcat_ref_short.mag'

    cmd = 'grep "*" ' + fullcat_mag + ' > ' + fullcat_mag_short
    os.system(cmd)

    cmd = 'sed -i.bak "s/INDEF/99.999/g" ' + fullcat_mag_short
    os.system(cmd)

    # Remove .bak files to prevent confusion
    bak_fullcat = fullcat_mag_short + '.bak'
    os.remove(bak_fullcat)

    print ''
    print 'Use the plot "plot_growth_curve.pdf" to determine which stars/clusters'
    print 'should be ignored in the isolated stars/clusters catalogs.'
    print 'Comment them out with "#" and run Step 2 again.'
    print ''
    print 'Use the plot "plot_histogram_ci.pdf" to determine an appropriate CI limit.'
    print 'Write this value in line 8 of the input file.'
    print ''

    sys.exit('STEP 2 is finished')



###################    STEP 3    ###################


if  flagstep3 == 'yes':


#####  CI histogram plot of sextractor catalog  #####

    # Calculate CI for full sextractor catalog in reference filter
    fullcat_mag_short = target_dir + '/photometry/ci_fullcat_ref_short.mag'

    aper, mag = np.loadtxt(fullcat_mag_short, unpack=True, usecols=(0,4))

    # Calculate CI for all objects in sextractor catalog
    mag_1 = mag[aper == 1.0]
    mag_3 = mag[aper == 3.0]
    ci_fullcat = mag_1 - mag_3

    # ci is CI limit determined by user, turn from string to float
    ci = float(ci)

    # Create histogram of all CI values from sextractor catalog
    fig = plt.figure(figsize = (7,7))

    # Draw histogram
    num, bins, patches = plt.hist(ci_fullcat, bins = 50, range=(0.1,4.0), color='blue')
    plt.vlines(ci, 0, 10000, color='red', linewidth=2)

    plt.xlabel('CI(1-3 px) (mag)')
    plt.ylabel('N')
    plt.axis([0.1,4,0, max(num) + 0.1*max(num)])
    plt.minorticks_on()

    plt.savefig(target_dir + '/plots/plot_histogram_ci_fullcat.pdf')


    # Save the centered xy coordinates from CI photometry of sextractor catalog
    fullcat_mag = target_dir + '/photometry/ci_fullcat_ref.mag'
    fullcat_center_coo = target_dir + '/photometry/fullcat_ref_center.coo'

    iraf.txdump(fullcat_mag,'XCENTER,YCENTER','yes',Stdout=fullcat_center_coo)

    xc, yc = np.loadtxt(fullcat_center_coo, unpack=True)

    # ci is CI limit determined by user
    xc_cicut = xc[ci_fullcat >= ci]
    yc_cicut = yc[ci_fullcat >= ci]
    ci_cicut = ci_fullcat[ci_fullcat >= ci]

    # Create coordinate file of centered xy coordinates of sources satisfying CI cut
    cicut_center_coo = target_dir + '/photometry/fullcat_ref_center_cicut.coo'
    f = open(cicut_center_coo, 'w')

    for k in range(len(xc_cicut)):
        f.write(str(xc_cicut[k]) + '  ' + str(yc_cicut[k]) + '\n')
    f.close()

    # Create file of CI values of sources satisfying CI cut
    cicut_center_cimag = target_dir + '/photometry/civalues_fullcat_ref_center_cicut.mag'
    f = open(cicut_center_cimag, 'w')

    for k in range(len(xc_cicut)):
        f.write(str(ci_cicut[k]) + '\n')
    f.close()


#####   Photometry on CI-limited catalog in every filter - Science aperture   #####


    # A few setup things
    imlist = glob.glob(target_dir + '/img/*_sci.fits')

    cicut_center_coo = target_dir + '/photometry/fullcat_ref_center_cicut.coo'

    inst_zp, filter_zp, zp_zp = np.loadtxt(pydir + '/data/legus_zeropoints.tab', unpack=True, dtype='str')

    # Make a file containing the filters, instruments, zeropoints, and exptimes for
    # user reference
    header_info = target_dir + '/photometry/header_info.txt'
    file = open(header_info, 'w')

    # Remove output files in case they exist
    for image in imlist:

        filter = image.split('/')[-1].split('_')[2]

        if os.path.exists(target_dir + '/photometry/' + filter + '_4px_mode.mag') == True:
            os.remove(target_dir + '/photometry/' + filter + '_4px_mode.mag')

        if os.path.exists(target_dir + '/photometry/' + filter + '_4px_mode_short.mag') == True:
            os.remove(target_dir + '/photometry/' + filter + '_4px_mode_short.mag')

        if os.path.exists(target_dir + '/photometry/' + filter + '_4px_mode_final.mag') == True:
            os.remove(target_dir + '/photometry/' + filter + '_4px_mode_final.mag')


    # Set some parameters common to all filters
    iraf.centerpars.calgorithm = 'none'  ### no centering

    iraf.fitskypars.annulus = annulus_sci
    iraf.fitskypars.dannulu = dannulus_sci

    # Aperture photometry with an aperture of 4 px

    for image in imlist:

        # Get filter from filename
        filter = image.split('/')[-1].split('_')[2]

        # Set the necessary variables for photometry on the reference image
        exptime = pyfits.getheader(image)['EXPTIME']
        inst = pyfits.getheader(image)['INSTRUME']
        inst = inst.lower()

        match = (inst_zp == inst) & (filter_zp == filter)
        zp = zp_zp[match]

        # zp is a string within an array, so need to turn into a float
        zp = float(zp[0])

        # Save information in header info text file
        file.write(filter + '  ' + inst + '  ' + str(zp) + '  ' + str(exptime) + '\n')

        # Specify names of output files
        cicut_4px = target_dir + '/photometry/' + filter + '_4px_mode.mag'
        cicut_4px_short = target_dir + '/photometry/' + filter + '_4px_mode_short.mag'

        # Do photometry
        iraf.datapars.epadu = exptime

        iraf.photpars.apertures = useraperture
        iraf.photpars.zmag = zp

        iraf.phot(image, cicut_center_coo, cicut_4px)

        iraf.txdump(cicut_4px, 'XCENTER,YCENTER,FLUX,MAG,MERR,MSKY,ID', 'yes', Stdout = cicut_4px_short)

        # Set placeholders for sources outside of FOV and undetected sources
        # For outside of FOV, use 66.666 instead of INDEF
        # For undetected sources, use 99.999 instead of INDEF

        # Sources outside of FOV have exactly zero flux
        x, y, flux, mag, merr, msky, id = np.loadtxt(cicut_4px_short, unpack = True,
                                                        dtype = str)

        flux = flux.astype(float)

        out_fov = (flux == 0.)

        mag[out_fov] = 66.666
        merr[out_fov] = 66.666

        # Undetected sources, those with negative flux or fluxes so small that mag err
        # is INDEF
        neg_flux = (flux < 0.)
        tiny_flux = (flux > 0.) & (merr == 'INDEF')

        mag[neg_flux] = 99.999
        merr[neg_flux] = 99.999

        merr[tiny_flux] = 99.999

        # Save results to new file
        x = x.astype(float)
        y = y.astype(float)
        mag = mag.astype(float)
        merr = merr.astype(float)
        msky = msky.astype(float)
        id = id.astype(int)

        cicut_4px_final = target_dir + '/photometry/' + filter + '_4px_mode_final.mag'

        zip_phot = zip(x, y, mag, merr, msky, id)

        np.savetxt(cicut_4px_final, zip_phot,
                    fmt = '%.3f  %.3f  %.3f  %.3f  %.9f  %i')

        # Get rid of INDEFs
#         cmd = 'sed -i.bak "s/INDEF/99.999/g" ' + cicut_4px_short
#         os.system(cmd)
#
#         bak_cicut_4px_short = cicut_4px_short + '.bak'
#         os.remove(bak_cicut_4px_short)

    file.close()


#####  Photometry on isolated clusters for aperture corrections  #####


    # Remove output files in case they exist
    for image in imlist:

        filter = image.split('/')[-1].split('_')[2]

        if os.path.exists(target_dir + '/photometry/' + filter + '_4px_apcor.mag') == True:
            os.remove(target_dir + '/photometry/' + filter + '_4px_apcor.mag')

        if os.path.exists(target_dir + '/photometry/' + filter + '_4px_apcor_short.mag') == True:
            os.remove(target_dir + '/photometry/' + filter + '_4px_apcor_short.mag')

        if os.path.exists(target_dir + '/photometry/' + filter + '_20px_apcor.mag') == True:
            os.remove(target_dir + '/photometry/' + filter + '_20px_apcor.mag')

        if os.path.exists(target_dir + '/photometry/' + filter + '_20px_apcor_short.mag') == True:
            os.remove(target_dir + '/photometry/' + filter + '_20px_apcor_short.mag')


    # Set some parameters common to all filters
    iraf.centerpars.calgorithm = 'centroid' ###center for aperture corrections


    # Aperture photometry with an aperture of 4 px on isolated clusters chosen by
    # user

    for image in imlist:

        # Get filter from filename
        filter = image.split('/')[-1].split('_')[2]

        # Set the necessary variables for photometry on the reference image
        exptime = pyfits.getheader(image)['EXPTIME']
        inst = pyfits.getheader(image)['INSTRUME']
        inst = inst.lower()

        match = (inst_zp == inst) & (filter_zp == filter)
        zp = zp_zp[match]

        # zp is a string within an array, so need to turn into a float
        zp = float(zp[0])

        # Specify names of output files
        iso_4px = target_dir + '/photometry/' + filter + '_4px_apcor.mag'
        iso_4px_short = target_dir + '/photometry/' + filter + '_4px_apcor_short.mag'

        # Do photometry
        iraf.datapars.epadu = exptime

        iraf.fitskypars.annulus = annulus_sci
        iraf.fitskypars.dannulu = dannulus_sci

        iraf.photpars.apertures = useraperture
        iraf.photpars.zmag = zp

        iraf.phot(image, target_dir + '/photometry/isolated_clusters.coo', iso_4px)

        iraf.txdump(iso_4px, 'XCENTER,YCENTER,MAG,MERR,MSKY,ID', 'yes', Stdout = iso_4px_short)

        # Get rid of INDEFs
        cmd = 'sed -i.bak "s/INDEF/99.999/g" ' + iso_4px_short
        os.system(cmd)

        bak_iso_4px_short = iso_4px_short + '.bak'
        os.remove(bak_iso_4px_short)


        # Aperture photometry with an aperture of 20 px

    for image in imlist:

        # Get filter from filename
        filter = image.split('/')[-1].split('_')[2]

        # Set the necessary variables for photometry on the reference image
        exptime = pyfits.getheader(image)['EXPTIME']
        inst = pyfits.getheader(image)['INSTRUME']
        inst = inst.lower()

        match = (inst_zp == inst) & (filter_zp == filter)
        zp = zp_zp[match]

        # zp is a string within an array, so need to turn into a float
        zp = float(zp[0])

        # Specify names of output files
        iso_20px = target_dir + '/photometry/' + filter + '_20px_apcor.mag'
        iso_20px_short = target_dir + '/photometry/' + filter + '_20px_apcor_short.mag'

        # Do photometry
        iraf.datapars.epadu = exptime

        iraf.fitskypars.annulus = 21.
        iraf.fitskypars.dannulu = 1.

        iraf.photpars.apertures = 20.0
        iraf.photpars.zmag = zp

        iraf.phot(image, target_dir + '/photometry/isolated_clusters.coo', iso_20px)

        iraf.txdump(iso_20px, 'XCENTER,YCENTER,MAG,MERR,MSKY,ID', 'yes', Stdout=iso_20px_short)

        # Get rid of INDEFs
        cmd = 'sed -i.bak "s/INDEF/99.999/g" ' + iso_20px_short
        os.system(cmd)

        bak_iso_20px_short = iso_20px_short + '.bak'
        os.remove(bak_iso_20px_short)


    print ''
    print 'Check the plot "plot_histogram_ci_fullcat.pdf" for appropriate CI limit.'
    print ''

    sys.exit('STEP 3 is finished')



###################   STEP 4   ######################


#####  Calculate average aperture corrections  #####

if  flagstep4 == 'yes':

    imlist = glob.glob(target_dir + '/img/*_sci.fits')


    apcor_text = target_dir + '/photometry/avg_aperture_correction.txt'
    file = open(apcor_text, 'w')

    for image in imlist:

        # Get filter from filename
        filter = image.split('/')[-1].split('_')[2]

        iso_4px_short = target_dir + '/photometry/' + filter + '_4px_apcor_short.mag'
        iso_20px_short = target_dir + '/photometry/' + filter + '_20px_apcor_short.mag'

        x, y, mag_4 = np.loadtxt(iso_4px_short, unpack=True, usecols=(0,1,2))
        mag_20 = np.loadtxt(iso_20px_short, unpack=True, usecols=(2,))

        # Calculate aperture corrections
        apcor = mag_20 - mag_4

        # Limit range of aperture corrections allowed to go into average
        lim = (apcor < uplim_apcor) & (apcor > lolim_apcor)
        apcor_lim = apcor[lim]
        apcor_avg = np.mean(apcor[lim])
        apcor_err = np.std(apcor_lim)/np.sqrt(len(apcor_lim))

        # write out results to file
        file.write(filter + '  ' + str(apcor_avg) + '  ' + str(apcor_err) + '  ' + str(len(apcor_lim)) + '\n')


        print ''
        print 'filter: ' + filter
        print 'apcor: %.3f' % apcor_avg
        print 'apcor error: %.3f' % apcor_err


        # Make a plot for each filter
        fig = plt.figure(figsize = (7,7))

        # Draw histogram
        num, bins, patches = plt.hist(apcor, bins=50, range=(-4.0,1.0), color='red', lw=2, histtype='step')

        plt.vlines(apcor_avg, 0, 1, color='blue', linewidth=5)
        plt.vlines(uplim_apcor, 0, 2, color='black', linestyle='--', linewidth=2)
        plt.vlines(lolim_apcor, 0, 2, color='black', linestyle='--', linewidth=2)

        plt.xlabel('Aperture Correction')
        plt.ylabel('N')
        plt.minorticks_on()
        plt.axis([-4.0,1.0, 0, max(num) + 0.1*max(num)])

        fig.savefig(target_dir + '/plots/plot_apcor_' + filter + '.pdf')

    file.close()


#####  CI-based aperture corrections calculation  #####

    # Do CI photometry in all 5 bands
    imlist = glob.glob(target_dir + '/img/*_sci.fits')

    fullcat_center_cicut = target_dir + '/photometry/fullcat_ref_center_cicut.coo'

    inst_zp, filter_zp, zp_zp = np.loadtxt(pydir + '/data/legus_zeropoints.tab', unpack=True, dtype='str')

    # Load in appropriate CI-apcor relation from Dave
    # Find integer version of aperture radius
    aper_int = '%1i' % float(useraperture)

    instr_rel, c0, c1, c2, c3, range_min, range_max = np.loadtxt(pydir + '/data/ci_apcor_rel_' + aper_int + 'px.dat', unpack=True, dtype='str')

    # Erase CI files if they already exist
    for image in imlist:

        filter = image.split('/')[-1].split('_')[2]

        if os.path.exists(target_dir + '/photometry/' + filter + '_ci.mag') == True:
                os.remove(target_dir + '/photometry/' + filter + '_ci.mag')

        if os.path.exists(target_dir + '/photometry/' + filter + '_ci_short.mag') == True:
                os.remove(target_dir + '/photometry/' + filter + '_ci_short.mag')


    # Set some parameters common to all filters
    iraf.centerpars.calgorithm = 'none'

    for image in imlist:

        # Get filter from filename
        filter = image.split('/')[-1].split('_')[2]

        # Set the necessary variables for photometry on the reference image
        exptime = pyfits.getheader(image)['EXPTIME']
        inst = pyfits.getheader(image)['INSTRUME']
        inst = inst.lower()

        match = (inst_zp == inst) & (filter_zp == filter)
        zp = zp_zp[match]

        # zp is a string within an array, so need to turn into a float
        zp = float(zp[0])

        # Specify names of output files
        cicut_mag = target_dir + '/photometry/' + filter + '_ci.mag'
        cicut_mag_short = target_dir + '/photometry/' + filter + '_ci_short.mag'

        # Do photometry
        iraf.datapars.epadu = exptime

        iraf.fitskypars.annulus = annulus_sci
        iraf.fitskypars.dannulu = dannulus_sci

        iraf.photpars.apertures = '1.0,3.0'
        iraf.photpars.zmag = zp

        iraf.phot(image, fullcat_center_cicut, cicut_mag)

        iraf.txdump(cicut_mag, 'MAG[1],MERR[1],MAG[2],MERR[2]', 'yes', Stdout = cicut_mag_short)

        # Get rid of INDEFs
        cmd = 'sed -i.bak "s/INDEF/99.999/g" ' + cicut_mag_short
        os.system(cmd)

        bak_cicut_mag_short = cicut_mag_short + '.bak'
        os.remove(bak_cicut_mag_short)


        # Calculate apcor in each filter based on CI-apcor relations from Dave

        # Get filter from filename
        filter = image.split('/')[-1].split('_')[2]

        # Set the necessary variables for photometry on the reference image
        inst = pyfits.getheader(image)['INSTRUME']
        inst = inst.lower()

        # Pick right relation for instrument
        ind = (instr_rel == inst)
        c0_inst = float(c0[ind][0])
        c1_inst = float(c1[ind][0])
        c2_inst = float(c2[ind][0])
        c3_inst = float(c3[ind][0])
        range_min_inst = float(range_min[ind][0])
        range_max_inst = float(range_max[ind][0])

        # Read in CI photometry
        cicut_mag_short = target_dir + '/photometry/' + filter + '_ci_short.mag'
        mag_1, merr_1, mag_3, merr_3 = np.loadtxt(cicut_mag_short, unpack=True)

        ci_cicut = mag_1 - mag_3
        cierr_cicut = np.sqrt(merr_1**2 + merr_3**2)

        # Set CIs outside the range min and max from Dave equal to the range min and max
        large = (ci_cicut > range_max_inst)
        small = (ci_cicut < range_min_inst)

        ci_cicut[large] = range_max_inst
        ci_cicut[small] = range_min_inst

        # Calculate apcor based on CIs
        # Use relation: apcor = c0 + c1*CI + c2*CI**2 + c3*CI**3
        apcor_fromci = c0_inst + c1_inst*ci_cicut + c2_inst*ci_cicut**2 + c3_inst*ci_cicut**3

        # Find apcor errors.
        # Read in appropriate artificial clusters apcor standard deviation file
        ci_artif, stdev_apcor = np.loadtxt(pydir + '/data/artificial_stdev_' + inst + '.dat', unpack=True)

        # Pick apcor standard deviation closest to CI of object
        best_stdev_apcor = np.empty(len(ci_cicut))

        for i in np.arange(len(ci_cicut)):

            min_diff_ind = np.abs(ci_artif - ci_cicut[i]).argmin()
            best_stdev_apcor[i] = stdev_apcor[min_diff_ind]


        # If either mag_1 or mag_3 were 99.999, set apcor back to 99.999
        no99 = (mag_1 == 99.999) | (mag_3 == 99.999)

        apcor_fromci[no99] = 99.999
        best_stdev_apcor[no99] = 99.999

        # Save apcors to file
        zip_apcor_fromci = zip(apcor_fromci, best_stdev_apcor)

        file_apcor_fromci = target_dir + '/photometry/apcor_fromci_' + filter + '.mag'
        np.savetxt(file_apcor_fromci, zip_apcor_fromci, fmt='%.3f  %.3f')



    print ''
    print 'Check the aperture correction plots for each filter in the plots/ directory.'
    print 'Update lines 10 and 11 in the input file and rerun Step 4 as necessary.'
    print ''

    sys.exit('STEP 4 is finished')



####################   STEP 5   #####################


if  flagstep5 == 'yes':


#####  Final photometric catalogues - SETUP  #####

    # Read in reference filter
    ref_filter = np.loadtxt(target_dir + '/s_extraction/ref_filter.txt', dtype='str')
    ref_filter = ref_filter[()]

    # Read in CI values in reference filter for CI-limited catalog
    cicut_civalues = target_dir + '/photometry/civalues_fullcat_ref_center_cicut.mag'
    ci_ref = np.loadtxt(cicut_civalues)

    # Read in science photometry for CI-limited catalog (user defined aperture)

    # Get list of files with '4px_mode_final.mag' suffixes
    filelist = glob.glob(target_dir + '/photometry/*_4px_mode_final.mag')

    # Find files corresponding to filter names
    file_2 = [file for file in filelist if 'f275w' in file][0]
    file_3 = [file for file in filelist if 'f336w' in file][0]
    file_8 = [file for file in filelist if 'f814w' in file][0]

    # Special cases
    file_4 = [file for file in filelist if 'f435w' in file or 'f438w' in file][0]
    file_5 = [file for file in filelist if 'f555w' in file or 'f606w' in file][0]

    # Read in the data
    xc5, yc5, mag5, err5, id5 = np.loadtxt(file_5, unpack=True, usecols=(0,1,2,3,5))
    mag2, err2 = np.loadtxt(file_2, unpack=True, usecols=(2,3))
    mag3, err3 = np.loadtxt(file_3, unpack=True, usecols=(2,3))
    mag4, err4 = np.loadtxt(file_4, unpack=True, usecols=(2,3))
    mag8, err8 = np.loadtxt(file_8, unpack=True, usecols=(2,3))


    # Select only sources that have low photometric errors in two neighboring filters
    sel_2n = ((err2 <= 0.3) & (err3 <= 0.3)) | ((err3 <= 0.3) & (err4 <= 0.3)) | ((err4 <= 0.3) & (err5 <= 0.3)) | ((err5 <= 0.3) & (err8 <= 0.3))

    xc5_2n = xc5[sel_2n]
    yc5_2n = yc5[sel_2n]
    mag2_2n = mag2[sel_2n]
    err2_2n = err2[sel_2n]
    mag3_2n = mag3[sel_2n]
    err3_2n = err3[sel_2n]
    mag4_2n = mag4[sel_2n]
    err4_2n = err4[sel_2n]
    mag5_2n = mag5[sel_2n]
    err5_2n = err5[sel_2n]
    mag8_2n = mag8[sel_2n]
    err8_2n = err8[sel_2n]
    ci_ref_2n = ci_ref[sel_2n]

    print 'Length of catalog with duplicates: ', len(xc5_2n)

    # Assign number of filters to each cluster
    # Must have low photometric errors in each filter considered
    nfilt_avg = np.zeros(len(xc5_2n))

    nfilt_4 = ((err2_2n <= 0.3) & (err4_2n <= 0.3) & (err5_2n <= 0.3) &
                (err8_2n <= 0.3)) | ((err3_2n <= 0.3) & (err4_2n <= 0.3) &
                (err5_2n <= 0.3) & (err8_2n <= 0.3))

    nfilt_5 = ((err2_2n <= 0.3) & (err3_2n <= 0.3) & (err4_2n <= 0.3) & (err5_2n <= 0.3) &
                (err8_2n <= 0.3))

    nfilt_avg[nfilt_4] = 4
    nfilt_avg[nfilt_5] = 5


    # Remove duplicate objects (distance between object centers less than 0.5 pix)
    tolerance = 0.5
    for_removal = []

    for k in np.arange(len(xc5_2n) - 1):

        # Calculate distance between kth object and subsequent objects
        d = np.sqrt((xc5_2n[k] - xc5_2n[k + 1:])**2 + (yc5_2n[k] - yc5_2n[k + 1:])**2)

        # If any of the distances calculate is less than the tolerance, save k
        if (d < tolerance).any():

            for_removal.append(k)

    print 'Indices of duplicate objects to be removed: ', for_removal

    # Generate an array of True booleans of same length as data
    bool = np.ones(len(xc5_2n), dtype=bool)

    # Change elements in bool that correspond to indices in for_removal to False
    bool[for_removal] = False

    # Apply bool to data to remove duplicates
    xc5_2n = xc5_2n[bool]
    yc5_2n = yc5_2n[bool]
    mag2_2n = mag2_2n[bool]
    err2_2n = err2_2n[bool]
    mag3_2n = mag3_2n[bool]
    err3_2n = err3_2n[bool]
    mag4_2n = mag4_2n[bool]
    err4_2n = err4_2n[bool]
    mag5_2n = mag5_2n[bool]
    err5_2n = err5_2n[bool]
    mag8_2n = mag8_2n[bool]
    err8_2n = err8_2n[bool]
    ci_ref_2n = ci_ref_2n[bool]
    nfilt_avg = nfilt_avg[bool]

    print 'Length of catalog without duplicates: ', len(xc5_2n)

    # Assign the galactic extinction for the five instrument/filter combinations
    # Read in extinction file as array (a vs u means acs/uvis)
    target_ext, uf2, uf3, af4, uf4, af5, uf5, af6, af8, uf8 = np.loadtxt(pydir + '/data/legus_galactic_extinction.tab', unpack=True, skiprows=2, dtype='str')

    # Pick extinctions for our target only
    id = (target_ext == target)

    imlist = glob.glob(target_dir + '/img/*_sci.fits')

    # Save extinction ordered by filter wavelength
    gal_ext = np.zeros(5)

    for image in imlist:

        # Get filter from filename
        filter = image.split('/')[-1].split('_')[2]

        # Get instrument from header
        instr = pyfits.getheader(image)['INSTRUME']
        instr = instr.lower()

        if (instr == 'wfc3') & (filter == 'f275w'):
            gal_ext[0] = float(uf2[id][0])

        if (instr == 'wfc3') & (filter == 'f336w'):
            gal_ext[1] = float(uf3[id][0])

        if (instr == 'acs') & (filter == 'f435w'):
            gal_ext[2] = float(af4[id][0])

        if (instr == 'wfc3') & (filter == 'f438w'):
            gal_ext[2] = float(uf4[id][0])

        if (instr == 'acs') & (filter == 'f555w'):
            gal_ext[3] = float(af5[id][0])

        if (instr == 'wfc3') & (filter == 'f555w'):
            gal_ext[3] = float(uf5[id][0])

        if (instr == 'acs') & (filter == 'f606w'):
            gal_ext[3] = float(af6[id][0])

        if (instr == 'acs' ) & (filter == 'f814w'):
            gal_ext[4] = float(af8[id][0])

        if (instr == 'wfc3') & (filter == 'f814w'):
            gal_ext[4] = float(uf8[id][0])


    # Make an array of cluster ids
    ids = np.arange(len(xc5_2n)) + 1

    # Convert xy coordinates into RA Dec of reference filter
    ref_image = [image for image in imlist if ref_filter in image][0]

    # Get header from reference image using pyfits
    header_ref = pyfits.getheader(ref_image)

    # Get wcs solution from reference image header
    wcs_ref = pywcs.WCS(header_ref)

    # Calculate RA and Dec for xy coordinates. 1 refers to origin of image in ds9.
    ra_2n, dec_2n = wcs_ref.wcs_pix2sky(xc5_2n, yc5_2n, 1)

    # Calculate apparent magnitude of absolute mag limit (-6)
    print ''
    print 'distance (Mpc) to target = ', distance

    distance = float(distance)

    dist_mod = 5 * math.log(distance,10) + 25.
    m_apparent = dist_mod - 6.



#####  Make Average Aperture Correction Final Catalog  #####

    # Sort aperture correction file by filter wavelength
    filter, avg_apcor, err_avg_apcor, nclusters = np.loadtxt(target_dir + '/photometry/avg_aperture_correction.txt', unpack=True, dtype='str')

    ind = filter.argsort()
    avg_apcor_sort = avg_apcor[ind].astype(float)
    err_avg_apcor_sort = err_avg_apcor[ind].astype(float)

    # Make catalog of all sources that meet CI limit (including magnitude,
    # CI, and nfilt) and a catalog of xy coordinates of these sources

    # Apply aperture corrections and galactic extinction to photometry

    mag2_avg = mag2_2n + avg_apcor_sort[0] - gal_ext[0]
    mag3_avg = mag3_2n + avg_apcor_sort[1] - gal_ext[1]
    mag4_avg = mag4_2n + avg_apcor_sort[2] - gal_ext[2]
    mag5_avg = mag5_2n + avg_apcor_sort[3] - gal_ext[3]
    mag8_avg = mag8_2n + avg_apcor_sort[4] - gal_ext[4]

    err2_avg = np.sqrt(err2_2n**2 + err_avg_apcor_sort[0]**2)
    err3_avg = np.sqrt(err3_2n**2 + err_avg_apcor_sort[1]**2)
    err4_avg = np.sqrt(err4_2n**2 + err_avg_apcor_sort[2]**2)
    err5_avg = np.sqrt(err5_2n**2 + err_avg_apcor_sort[3]**2)
    err8_avg = np.sqrt(err8_2n**2 + err_avg_apcor_sort[4]**2)

    # If anything was set to 66.666 or 99.999 before, set them back to those values

    # These are sources that fall outside the FOV, both mag and merr are set to 66.666
    fov_m2 = (mag2_2n == 66.666)
    fov_m3 = (mag3_2n == 66.666)
    fov_m4 = (mag4_2n == 66.666)
    fov_m5 = (mag5_2n == 66.666)
    fov_m8 = (mag8_2n == 66.666)

    mag2_avg[fov_m2] = 66.666
    mag3_avg[fov_m3] = 66.666
    mag4_avg[fov_m4] = 66.666
    mag5_avg[fov_m5] = 66.666
    mag8_avg[fov_m8] = 66.666

    err2_avg[fov_m2] = 66.666
    err3_avg[fov_m3] = 66.666
    err4_avg[fov_m4] = 66.666
    err5_avg[fov_m5] = 66.666
    err8_avg[fov_m8] = 66.666

    # Sources with negative fluxes have both mag and merr set to 99.999
    neg_m2 = (mag2_2n == 99.999)
    neg_m3 = (mag3_2n == 99.999)
    neg_m4 = (mag4_2n == 99.999)
    neg_m5 = (mag5_2n == 99.999)
    neg_m8 = (mag8_2n == 99.999)

    mag2_avg[neg_m2] = 99.999
    mag3_avg[neg_m3] = 99.999
    mag4_avg[neg_m4] = 99.999
    mag5_avg[neg_m5] = 99.999
    mag8_avg[neg_m8] = 99.999

    # Sources with very small positive fluxes may have only merr set to 99.999
    tiny_m2 = (err2_2n == 99.999)
    tiny_m3 = (err3_2n == 99.999)
    tiny_m4 = (err4_2n == 99.999)
    tiny_m5 = (err5_2n == 99.999)
    tiny_m8 = (err8_2n == 99.999)

    err2_avg[tiny_m2] = 99.999
    err3_avg[tiny_m3] = 99.999
    err4_avg[tiny_m4] = 99.999
    err5_avg[tiny_m5] = 99.999
    err8_avg[tiny_m8] = 99.999


    # Save to compare to automatic_catalog_target.mag
    zip_table2 = zip(ids, xc5_2n, yc5_2n, ra_2n, dec_2n, mag2_avg, err2_avg, mag3_avg, err3_avg, mag4_avg, err4_avg, mag5_avg, err5_avg, mag8_avg, err8_avg, ci_ref_2n, nfilt_avg)

    np.savetxt(target_dir + '/photometry/automatic_catalog_' + target + '_avgapcor.tab', zip_table2, fmt='%i  %.3f  %.3f  %.8f  %.8f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %i')



#####  Make CI-based Aperture Correction Final Catalog  #####

    apcorlist = glob.glob(target_dir + '/photometry/apcor_fromci_*.mag')

    # Find apcor files corresponding to filter names
    apcorfile_2 = [file for file in apcorlist if 'f275w' in file][0]
    apcorfile_3 = [file for file in apcorlist if 'f336w' in file][0]
    apcorfile_8 = [file for file in apcorlist if 'f814w' in file][0]

    # Special cases
    apcorfile_4 = [file for file in apcorlist if 'f435w' in file or 'f438w' in file][0]
    apcorfile_5 = [file for file in apcorlist if 'f555w' in file or 'f606w' in file][0]

    # Read in CI-based aperture corrections
    apcor2, err_apcor2 = np.loadtxt(apcorfile_2, unpack=True)
    apcor3, err_apcor3 = np.loadtxt(apcorfile_3, unpack=True)
    apcor4, err_apcor4 = np.loadtxt(apcorfile_4, unpack=True)
    apcor5, err_apcor5 = np.loadtxt(apcorfile_5, unpack=True)
    apcor8, err_apcor8 = np.loadtxt(apcorfile_8, unpack=True)

    # Must apply low photometric error selection to CI-based apcors too, otherwise array
    # lengths don't match!
    sel_2n = ((err2 <= 0.3) & (err3 <= 0.3)) | ((err3 <= 0.3) & (err4 <= 0.3)) | ((err4 <= 0.3) & (err5 <= 0.3)) | ((err5 <= 0.3) & (err8 <= 0.3))

    # Must also remove duplicates, so apply bool array to aperture corrections as well
    # (see lines 1151-1187)

    apcor2_2n = apcor2[sel_2n]
    apcor3_2n = apcor3[sel_2n]
    apcor4_2n = apcor4[sel_2n]
    apcor5_2n = apcor5[sel_2n]
    apcor8_2n = apcor8[sel_2n]

    apcor2_2n = apcor2_2n[bool]
    apcor3_2n = apcor3_2n[bool]
    apcor4_2n = apcor4_2n[bool]
    apcor5_2n = apcor5_2n[bool]
    apcor8_2n = apcor8_2n[bool]

    err_apcor2_2n = err_apcor2[sel_2n]
    err_apcor3_2n = err_apcor3[sel_2n]
    err_apcor4_2n = err_apcor4[sel_2n]
    err_apcor5_2n = err_apcor5[sel_2n]
    err_apcor8_2n = err_apcor8[sel_2n]

    err_apcor2_2n = err_apcor2_2n[bool]
    err_apcor3_2n = err_apcor3_2n[bool]
    err_apcor4_2n = err_apcor4_2n[bool]
    err_apcor5_2n = err_apcor5_2n[bool]
    err_apcor8_2n = err_apcor8_2n[bool]

    # Apply CI-based apcors to photometry

    mag2_cib = mag2_2n + apcor2_2n - gal_ext[0]
    mag3_cib = mag3_2n + apcor3_2n - gal_ext[1]
    mag4_cib = mag4_2n + apcor4_2n - gal_ext[2]
    mag5_cib = mag5_2n + apcor5_2n - gal_ext[3]
    mag8_cib = mag8_2n + apcor8_2n - gal_ext[4]

    err2_cib = np.sqrt(err2_2n**2 + err_apcor2_2n**2)
    err3_cib = np.sqrt(err3_2n**2 + err_apcor3_2n**2)
    err4_cib = np.sqrt(err4_2n**2 + err_apcor4_2n**2)
    err5_cib = np.sqrt(err5_2n**2 + err_apcor5_2n**2)
    err8_cib = np.sqrt(err8_2n**2 + err_apcor8_2n**2)

    # If magnitudes were 66.666 or 99.999 before, set them back to 66.666 or 99.999

    # These are sources that fall outside the FOV, both mag and merr are set to 66.666
    fov_m2 = (mag2_2n == 66.666)
    fov_m3 = (mag3_2n == 66.666)
    fov_m4 = (mag4_2n == 66.666)
    fov_m5 = (mag5_2n == 66.666)
    fov_m8 = (mag8_2n == 66.666)

    mag2_cib[fov_m2] = 66.666
    mag3_cib[fov_m3] = 66.666
    mag4_cib[fov_m4] = 66.666
    mag5_cib[fov_m5] = 66.666
    mag8_cib[fov_m8] = 66.666

    err2_cib[fov_m2] = 66.666
    err3_cib[fov_m3] = 66.666
    err4_cib[fov_m4] = 66.666
    err5_cib[fov_m5] = 66.666
    err8_cib[fov_m8] = 66.666

    # Do same for sources set to 99.999, but mag and err separately
    # This is because sometimes merr is 99.999, but mag is not
    neg_m2 = (mag2_2n == 99.999)
    neg_m3 = (mag3_2n == 99.999)
    neg_m4 = (mag4_2n == 99.999)
    neg_m5 = (mag5_2n == 99.999)
    neg_m8 = (mag8_2n == 99.999)

    mag2_cib[neg_m2] = 99.999
    mag3_cib[neg_m3] = 99.999
    mag4_cib[neg_m4] = 99.999
    mag5_cib[neg_m5] = 99.999
    mag8_cib[neg_m8] = 99.999

    tiny_m2 = (err2_2n == 99.999)
    tiny_m3 = (err3_2n == 99.999)
    tiny_m4 = (err4_2n == 99.999)
    tiny_m5 = (err5_2n == 99.999)
    tiny_m8 = (err8_2n == 99.999)

    err2_cib[tiny_m2] = 99.999
    err3_cib[tiny_m3] = 99.999
    err4_cib[tiny_m4] = 99.999
    err5_cib[tiny_m5] = 99.999
    err8_cib[tiny_m8] = 99.999

    # If apcorr is 99.999 but photometry is detected (less than 66.666), set to 44.444
    ap_m2 = (apcor2_2n == 99.999) & (mag2_2n < 66.666)
    ap_m3 = (apcor3_2n == 99.999) & (mag3_2n < 66.666)
    ap_m4 = (apcor4_2n == 99.999) & (mag4_2n < 66.666)
    ap_m5 = (apcor5_2n == 99.999) & (mag5_2n < 66.666)
    ap_m8 = (apcor8_2n == 99.999) & (mag8_2n < 66.666)

    mag2_cib[ap_m2] = 44.444
    mag3_cib[ap_m3] = 44.444
    mag4_cib[ap_m4] = 44.444
    mag5_cib[ap_m5] = 44.444
    mag8_cib[ap_m8] = 44.444

    err2_cib[ap_m2] = 44.444
    err3_cib[ap_m3] = 44.444
    err4_cib[ap_m4] = 44.444
    err5_cib[ap_m5] = 44.444
    err8_cib[ap_m8] = 44.444


    # Update number of filters for each cluster
    # Must be detected in each filter considered
    nfilt_cib = np.copy(nfilt_avg)

    nfilt_4 = (((mag2_cib < 44.444) & (mag4_cib < 44.444) & (mag5_cib < 44.444) & (mag8_cib < 44.444)) | ((mag3_cib < 44.444) & (mag4_cib < 44.444) & (mag5_cib < 44.444) & (mag8_cib < 44.444))) & ((nfilt_cib == 4) | (nfilt_cib == 5))

    nfilt_5 = ((mag2_cib < 44.444) & (mag3_cib < 44.444) & (mag4_cib < 44.444) & (mag5_cib < 44.444) & (mag8_cib < 44.444)) & (nfilt_cib == 5)

    nfilt_cib[nfilt_4] = 4
    nfilt_cib[nfilt_5] = 5
    nfilt_cib[~nfilt_5 & ~nfilt_4] = 0

    # Save to compare to automatic_catalog_target.mag
    zip_table3 = zip(ids, xc5_2n, yc5_2n, ra_2n, dec_2n, mag2_cib, err2_cib, mag3_cib, err3_cib, mag4_cib, err4_cib, mag5_cib, err5_cib, mag8_cib, err8_cib, ci_ref_2n, nfilt_cib)

    np.savetxt(target_dir + '/photometry/automatic_catalog_' + target + '_cibased.tab', zip_table3, fmt='%i  %.3f  %.3f  %.8f  %.8f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %i')


#####  Make visual inspection catalog  #####

    # Make coordinate file of clusters that have >= 4 filters, are brighter than
    # magnitude limit (when apcor-CI correction included), and CI bigger than limit

    # Calculate reference filter magnitude with aperture correction APPLIED
    # Aperture correction in this case is already based on the CI, so we use mag5_cib
    ci = float(ci)

    vis_ind2 = (nfilt_avg >= 4) & (mag5_cib <= m_apparent) & (ci_ref_2n >= ci)

    xc5_fin2 = xc5_2n[vis_ind2]
    yc5_fin2 = yc5_2n[vis_ind2]

    zip_fin2 = zip(xc5_fin2, yc5_fin2)

    np.savetxt(target_dir + '/photometry/cat_UBVI_avgapcor_magcut_ci' + str(ci) + '_' + target + '.coo', zip_fin2, fmt='%.3f  %.3f')

    # Make a regions file of final clusters
    outputfile = target_dir + '/photometry/catalog_ds9_clusters_avgapcor_final.reg'

    file = open(outputfile, 'w')
    file.write('global color=blue width=3 font="helvetica 15 normal roman" highlite=1 \n')
    file.write('image\n')

    for i in range(len(xc5_fin2)):

        x = str(xc5_fin2[i])
        y = str(yc5_fin2[i])
        newline = 'circle('+ x + ',' + y +  ',7)\n'
        file.write(newline)

    file.close()



#####  Make readme file  #####

    # Get necessary variables
    ref_filter_up = ref_filter.upper()

    imlist = glob.glob(target_dir + '/img/*_sci.fits')

    filterlist = np.chararray(len(imlist), itemsize=5)
    instlist = np.chararray(len(imlist), itemsize=5)
    i = 0

    for image in imlist:

        # Get filter from filename
        filter = image.split('/')[-1].split('_')[2]

        # Get instrument from header
        instr = pyfits.getheader(image)['INSTRUME']

        filterlist[i] = filter.upper()
        instlist[i] = instr + '/'

        i += 1

    sort_ind = filterlist.argsort()
    filt_inst_sort = instlist[sort_ind] + filterlist[sort_ind]

    print ''
    print filt_inst_sort



    # Make readme string
    readme = """# Photometry is in Vega system, using ZP from the STScI webpage.
# The aperture radius used for photometry is {aperture}.
# Galactic ext from NED and aperture correction is included.
# Distance modulus used {dist_mod:.2f} mag ({distance} Mpc) from NED.
# When 66.666 appears in a photometry column, it means the source was not located
# within the field of view of that frame.
# When 99.999 appears in a photometry column, it means the source was not detected,
# i.e., the flux was negative or extremely small.
# When 44.444 appears in a photometry column (only present in the 'cibased' catalog),
# it means the source was detected, but the CI-based aperture correction was undefined.
# The photometry in the 'avgapcor' catalog uses average aperture corrections
# calculated from isolated clusters in the field. The photometry in the 'cibased'
# catalog uses aperture corrections calculated from a CI-aperture correction
# relation derived from artificial clusters. The number of filters in which each
# source is detected may be different between the two catalogs because CI-based
# aperture corrections are occasionally unreliable for certain filters.

Order of the columns
1. source id
2. x coordinates in the {ref_filter} frame (image aligned and registered, these coordinates are the same in each filter)
3. y coordinates in the {ref_filter} frame (image aligned and registered, these coordinates are the same in each filter)
4. RA decimal coordinates in the {ref_filter} frame (image aligned and registered, these coordinates are the same in each filter)
5. Dec decimal coordinates in the {ref_filter} frame (image aligned and registered, these coordinates are the same in each filter)
6. final total mag in {filt1}
7. final photometric error in {filt1}
8. final total mag in {filt2}
9. final photometric error in {filt2}
10. final total mag in {filt3}
11. final photometric error in {filt3}
12. final total mag in {filt4}
13. final photometric error in {filt4}
14. final total mag in {filt5}
15. final photometric error in {filt5}
16. CI=mag(1px)-mag(3px) measured in the {ref_filter} frame. This catalogue contains only sources with CI >= {ci}.
17. Number of filters. Sources with UBVI or UV-BVI detection have Nflt=4; sources with detection in UV-UBVI have Nflt=5. The remaining sources have Nflt=0 and should be considered with some regard."""

    # Dictionary to format the readme string
    context = {
        "aperture": useraperture,
        "dist_mod": dist_mod,
        "distance": distance,
        "ref_filter": ref_filter_up,
        "filt1": filt_inst_sort[0],
        "filt2": filt_inst_sort[1],
        "filt3": filt_inst_sort[2],
        "filt4": filt_inst_sort[3],
        "filt5": filt_inst_sort[4],
        "ci": ci,
        }

    filename = target_dir + '/photometry/automatic_catalogue_' + target + '.readme'

    with open(filename, 'w') as f:
        f.write(readme.format(**context))
    f.close()


    # Make new final directory
    if os.path.exists(target_dir + '/final') == False:
        os.makedirs(target_dir + '/final')


    # move final files to final directory
    source = target_dir + '/photometry/automatic_catalog_' + target + '_avgapcor.tab'
    destination = target_dir + '/final/automatic_catalog_' + target + '_avgapcor.tab'
    os.rename(source, destination)

    source = target_dir + '/photometry/automatic_catalog_' + target + '_cibased.tab'
    destination = target_dir + '/final/automatic_catalog_' + target + '_cibased.tab'
    os.rename(source, destination)

    source = target_dir + '/photometry/cat_UBVI_avgapcor_magcut_ci' + str(ci) + '_' + target + '.coo'
    destination = target_dir + '/final/cat_UBVI_avgapcor_magcut_ci' + str(ci) + '_' + target + '.coo'
    os.rename(source, destination)

    source = target_dir + '/photometry/catalog_ds9_clusters_avgapcor_final.reg'
    destination = target_dir + '/final/catalog_ds9_clusters_avgapcor_final.reg'
    os.rename(source, destination)

    source = target_dir + '/photometry/automatic_catalogue_' + target + '.readme'
    destination = target_dir + '/final/automatic_catalogue_' + target + '.readme'
    os.rename(source, destination)


    print ''
    print 'Check final catalogs in /final directory.'
    print ''

    sys.exit('STEP 5 is finished')



##########################   STEP 6   ##########################

# we perform the same corrections and procedure
# using the input file listed in line 12 of the
# input file.
# The output files have th suffix _add


if  flagstep6 == 'yes':

#####  photometry of manually added clusters  #####

    # If no inputlist, no need to run Step 6
    if inputlist == 'none':

        print 'No manually added clusters. No need to run Step 6.'
        sys.exit('Quitting...')

    elif inputlist != 'none':

        # Define full path of input list
        inputlist_path = target_dir + '/' + inputlist

    # Set image for centering sources to reference frame if no smoothed image
    if smooth_image == 'none':

        imlist = glob.glob(target_dir + '/img/*_sci.fits')
        ref_filter = np.loadtxt(target_dir + '/s_extraction/ref_filter.txt', dtype='str')
        ref_filter = ref_filter[()]
        ref_image = [image for image in imlist if ref_filter in image][0]

        center_image = ref_image

    elif smooth_image != 'none':

        # Define full path of smooth_image
        center_image = target_dir + '/img/' + smooth_image
        ref_filter = np.loadtxt(target_dir + '/s_extraction/ref_filter.txt', dtype='str')
        ref_filter = ref_filter[()]


    # Make sure inputlist exists in target directory
    if os.path.exists(inputlist_path) == False:
        print 'Image ', inputlist_path ,' could not be found'

    # Make sure center_image exists in /img directory
    if os.path.exists(center_image) == False:
        print 'Image ', center_image ,' could not be found'


    # Remove output files in case they exist

    imlist = glob.glob(target_dir + '/img/*_sci.fits')

    for image in imlist:

        filter = image.split('/')[-1].split('_')[2]


        if os.path.exists(target_dir + '/photometry/' + filter + '_add.mag') == True:
            os.remove(target_dir + '/photometry/' + filter + '_add.mag')

        if os.path.exists(target_dir + '/photometry/' + filter + '_add_short.mag') == True:
            os.remove(target_dir + '/photometry/' + filter + '_add_short.mag')

        if os.path.exists(target_dir + '/photometry/' + filter + '_add_sci.mag') == True:
            os.remove(target_dir + '/photometry/' + filter + '_add_sci.mag')


    # Get necessary variables for reference/smoothed image (center_image)
    cen_exptime = pyfits.getheader(center_image)['EXPTIME']
    cen_inst = pyfits.getheader(center_image)['INSTRUME']
    cen_inst = cen_inst.lower()

    cen_filter = center_image.split('/')[-1].split('_')[2]

    inst_zp, filter_zp, zp_zp = np.loadtxt(pydir + '/data/legus_zeropoints.tab', unpack=True, dtype='str')

    match = (inst_zp == cen_inst) & (filter_zp == cen_filter)
    cen_zp = zp_zp[match]
    cen_zp = float(cen_zp[0])

    # Define output coordinate file
    center_mag = target_dir + '/photometry/center_add.mag'
    center_coo = target_dir + '/photometry/center_add.coo'

    # Run photometry to get centered coordinates in smoothed/reference image

    #get centroid from smoothed images on reference image. no centering in other filters
    iraf.centerpars.calgorithm = 'centroid'

    iraf.fitskypars.annulus = annulus_sci
    iraf.fitskypars.dannulu = dannulus_sci

    iraf.datapars.epadu = cen_exptime

    iraf.photpars.apertures = useraperture
    iraf.photpars.zmag = cen_zp

    iraf.phot(center_image, inputlist_path, center_mag)

    # Just output the coordinates
    iraf.txdump(center_mag, 'XCENTER,YCENTER', 'yes', Stdout=center_coo)


    # Perform photometry on manually added clusters using the centered coordinates
    # in all 5 frames
    iraf.centerpars.calgorithm = 'none'

    iraf.fitskypars.annulus = annulus_sci
    iraf.fitskypars.dannulu = dannulus_sci

    for image in imlist:

        exptime = pyfits.getheader(image)['EXPTIME']
        inst = pyfits.getheader(image)['INSTRUME']
        inst = inst.lower()

        filter = image.split('/')[-1].split('_')[2]

        match = (inst_zp == inst) & (filter_zp == filter)
        zp = zp_zp[match]
        zp = float(zp[0])

        # Define filenames
        add_mag = target_dir + '/photometry/' + filter + '_add.mag'
        add_mag_short = target_dir + '/photometry/' + filter + '_add_short.mag'
        add_mag_sci = target_dir + '/photometry/' + filter + '_add_sci.mag'

        # Do photometry
        iraf.datapars.epadu = exptime

        iraf.photpars.apertures = '1.0,3.0,' + useraperture
        iraf.photpars.zmag = zp

        iraf.phot(image, center_coo, add_mag)

        iraf.txdump(add_mag, 'MAG[1],MERR[1],MAG[2],MERR[2],FLUX[3],MAG[3],MERR[3]', 'yes', Stdout=add_mag_short)


        # Set placeholders for sources outside of FOV and undetected sources for science
        # photometry only
        # For outside of FOV, use 66.666 instead of INDEF
        # For undetected sources, use 99.999 instead of INDEF

        # Sources outside of FOV have exactly zero flux
        flux, mag, merr = np.loadtxt(add_mag_short, unpack = True, usecols=(4,5,6),
                                        dtype = str)

        flux = flux.astype(float)

        out_fov = (flux == 0.)

        mag[out_fov] = 66.666
        merr[out_fov] = 66.666

        # Undetected sources, those with negative flux or fluxes so small that mag err
        # is INDEF
        neg_flux = (flux < 0.)
        tiny_flux = (flux > 0.) & (merr == 'INDEF')

        mag[neg_flux] = 99.999
        merr[neg_flux] = 99.999

        merr[tiny_flux] = 99.999

        # Save results to new file
        mag = mag.astype(float)
        merr = merr.astype(float)

        add_mag_sci = target_dir + '/photometry/' + filter + '_add_sci.mag'

        zip_phot = zip(mag, merr)

        np.savetxt(add_mag_sci, zip_phot, fmt = '%.3f  %.3f')


        # Remove INDEFs from add_mag_short anyway

        cmd = 'sed -i.bak "s/INDEF/99.999/g" ' + add_mag_short
        os.system(cmd)

        bak_add_mag_short = add_mag_short + '.bak'
        os.remove(bak_add_mag_short)


##### Calculate apcor in each filter based on CI-apcor relations from Dave  #####

        # Load in appropriate CI-apcor relation from Dave
        # Find integer version of aperture radius
        aper_int = '%1i' % float(useraperture)
        instr_rel, c0, c1, c2, c3, range_min, range_max = np.loadtxt(pydir + '/data/ci_apcor_rel_' + aper_int + 'px.dat', unpack=True, dtype='str')

        # Pick right relation for instrument
        ind = (instr_rel == inst)
        c0_inst = float(c0[ind][0])
        c1_inst = float(c1[ind][0])
        c2_inst = float(c2[ind][0])
        c3_inst = float(c3[ind][0])
        range_min_inst = float(range_min[ind][0])
        range_max_inst = float(range_max[ind][0])

        # Read in CI photometry
        mag_1, merr_1, mag_3, merr_3 = np.loadtxt(add_mag_short, unpack=True, usecols=(0,1,2,3))

        ci_add = mag_1 - mag_3
        cierr_add = np.sqrt(merr_1**2 + merr_3**2)

        # Set CIs outside the range min and max from Dave equal to the range min and max
        large = (ci_add > range_max_inst)
        small = (ci_add < range_min_inst)

        ci_add[large] = range_max_inst
        ci_add[small] = range_min_inst

        # Calculate apcor based on CIs
        # Use relation: apcor = c0 + c1*CI + c2*CI**2 + c3*CI**3
        apcor_fromci = c0_inst + c1_inst*ci_add + c2_inst*ci_add**2 + c3_inst*ci_add**3

        # Find apcor errors.
        # Read in appropriate artificial clusters apcor standard deviation file
        ci_artif, stdev_apcor = np.loadtxt(pydir + '/data/artificial_stdev_' + inst + '.dat', unpack=True)

        # Pick apcor standard deviation closest to CI of object
        best_stdev_apcor = np.empty(len(ci_add))

        for i in np.arange(len(ci_add)):

            min_diff_ind = np.abs(ci_artif - ci_add[i]).argmin()
            best_stdev_apcor[i] = stdev_apcor[min_diff_ind]

        # If either mag_1 or mag_3 were 99.999, set apcor back to 99.999
        no99 = (mag_1 == 99.999) | (mag_3 == 99.999)

        apcor_fromci[no99] = 99.999
        best_stdev_apcor[no99] = 99.999

        # Save apcors to file
        zip_apcor_fromci = zip(apcor_fromci, best_stdev_apcor)

        file_apcor_fromci = target_dir + '/photometry/apcor_fromci_' + filter + '_add.mag'
        np.savetxt(file_apcor_fromci, zip_apcor_fromci, fmt='%.3f  %.3f')




####  Make final photometric catalogs for added clusters (separate files)  #####


#####  Final photometric catalogues - SETUP  #####

    # Read in photometry for additional catalog, CI and science aperture

    # Get list of files with 'add_sci.mag' suffixes
    filelist = glob.glob(target_dir + '/photometry/*_add_sci.mag')

    # Find files corresponding to filter names
    file_2 = [file for file in filelist if 'f275w' in file][0]
    file_3 = [file for file in filelist if 'f336w' in file][0]
    file_8 = [file for file in filelist if 'f814w' in file][0]

    # Special cases
    file_4 = [file for file in filelist if 'f435w' in file or 'f438w' in file][0]
    file_5 = [file for file in filelist if 'f555w' in file or 'f606w' in file][0]

    # Read in the data
    xc5, yc5 = np.loadtxt(target_dir + '/photometry/center_add.coo', unpack=True)
    mag5_sci, err5_sci = np.loadtxt(file_5, unpack=True)
    mag2_sci, err2_sci = np.loadtxt(file_2, unpack=True)
    mag3_sci, err3_sci = np.loadtxt(file_3, unpack=True)
    mag4_sci, err4_sci = np.loadtxt(file_4, unpack=True)
    mag8_sci, err8_sci = np.loadtxt(file_8, unpack=True)

    # Calculate ci values in reference filter for additional clusters
    cifilelist = glob.glob(target_dir + '/photometry/*_add_short.mag')

    cifile_5 = [file for file in cifilelist if 'f555w' in file or 'f606w' in file][0]

    mag5_a1, mag5_a3 = np.loadtxt(cifile_5, unpack=True, usecols=(0,2))

    ci5 = mag5_a1 - mag5_a3

    # Select only sources that have low photometric errors in two neighboring filters
    sel_2n = ((err2_sci <= 0.3) & (err3_sci <= 0.3)) | ((err3_sci <= 0.3) & (err4_sci <= 0.3)) | ((err4_sci <= 0.3) & (err5_sci <= 0.3)) | ((err5_sci <= 0.3) & (err8_sci <= 0.3))

    xc5_2n = xc5[sel_2n]
    yc5_2n = yc5[sel_2n]
    mag2_2n = mag2_sci[sel_2n]
    err2_2n = err2_sci[sel_2n]
    mag3_2n = mag3_sci[sel_2n]
    err3_2n = err3_sci[sel_2n]
    mag4_2n = mag4_sci[sel_2n]
    err4_2n = err4_sci[sel_2n]
    mag5_2n = mag5_sci[sel_2n]
    err5_2n = err5_sci[sel_2n]
    mag8_2n = mag8_sci[sel_2n]
    err8_2n = err8_sci[sel_2n]
    ci5_2n = ci5[sel_2n]

    # Should be no duplicates in a manually-added catalog, so do not perform same
    # procedure as lines 1151-1187

    # Assign number of filters to each cluster
    # Must have low photometric errors in each filter considered
    nfilt = np.zeros(len(xc5_2n))

    nfilt_4 = ((err2_2n <= 0.3) & (err4_2n <= 0.3) & (err5_2n <= 0.3) &
                (err8_2n <= 0.3)) | ((err3_2n <= 0.3) & (err4_2n <= 0.3) &
                (err5_2n <= 0.3) & (err8_2n <= 0.3))

    nfilt_5 = ((err2_2n <= 0.3) & (err3_2n <= 0.3) & (err4_2n <= 0.3) & (err5_2n <= 0.3) &
                (err8_2n <= 0.3))

    nfilt[nfilt_4] = 4
    nfilt[nfilt_5] = 5


    # Assign the galactic extinction for the five instrument/filter combinations
    # Read in extinction file as array (a vs u means acs/uvis)
    target_ext, uf2, uf3, af4, uf4, af5, uf5, af6, af8, uf8 = np.loadtxt(pydir + '/data/legus_galactic_extinction.tab', unpack=True, skiprows=2, dtype='str')

    # Pick extinctions for our target only
    id = (target_ext == target)

    imlist = glob.glob(target_dir + '/img/*_sci.fits')

    # Save extinction ordered by filter wavelength
    gal_ext = np.zeros(5)

    for image in imlist:

        # Get filter from filename
        filter = image.split('/')[-1].split('_')[2]

        # Get instrument from header
        instr = pyfits.getheader(image)['INSTRUME']
        instr = instr.lower()

        if (instr == 'wfc3') & (filter == 'f275w'):
            gal_ext[0] = float(uf2[id][0])

        if (instr == 'wfc3') & (filter == 'f336w'):
            gal_ext[1] = float(uf3[id][0])

        if (instr == 'acs') & (filter == 'f435w'):
            gal_ext[2] = float(af4[id][0])

        if (instr == 'wfc3') & (filter == 'f438w'):
            gal_ext[2] = float(uf4[id][0])

        if (instr == 'acs') & (filter == 'f555w'):
            gal_ext[3] = float(af5[id][0])

        if (instr == 'wfc3') & (filter == 'f555w'):
            gal_ext[3] = float(uf5[id][0])

        if (instr == 'acs') & (filter == 'f606w'):
            gal_ext[3] = float(af6[id][0])

        if (instr == 'acs' ) & (filter == 'f814w'):
            gal_ext[4] = float(af8[id][0])

        if (instr == 'wfc3') & (filter == 'f814w'):
            gal_ext[4] = float(uf8[id][0])


    # Make an array of cluster ids
    ids = np.arange(len(xc5_2n)) + 1

    # Convert xy coordinates into RA Dec of reference filter
    ref_image = [image for image in imlist if ref_filter in image][0]

    # Get header from reference image using pyfits
    header_ref = pyfits.getheader(ref_image)

    # Get wcs solution from reference image header
    wcs_ref = pywcs.WCS(header_ref)

    # Calculate RA and Dec for xy coordinates. 1 refers to origin of image in ds9.
    ra_2n, dec_2n = wcs_ref.wcs_pix2sky(xc5_2n, yc5_2n, 1)



#####  Make Average Aperture Correction Final Catalog  #####

    # Sort aperture correction file by filter wavelength
    filter, avg_apcor, err_avg_apcor, nclusters = np.loadtxt(target_dir + '/photometry/avg_aperture_correction.txt', unpack=True, dtype='str')

    ind = filter.argsort()
    avg_apcor_sort = avg_apcor[ind].astype(float)
    err_avg_apcor_sort = err_avg_apcor[ind].astype(float)

    # Make catalog of all sources that meet CI limit (including magnitude,
    # CI, and nfilt) and a catalog of xy coordinates of these sources

    # Apply aperture corrections and galactic extinction to photometry

    mag2_avg = mag2_2n + avg_apcor_sort[0] - gal_ext[0]
    mag3_avg = mag3_2n + avg_apcor_sort[1] - gal_ext[1]
    mag4_avg = mag4_2n + avg_apcor_sort[2] - gal_ext[2]
    mag5_avg = mag5_2n + avg_apcor_sort[3] - gal_ext[3]
    mag8_avg = mag8_2n + avg_apcor_sort[4] - gal_ext[4]

    err2_avg = np.sqrt(err2_2n**2 + err_avg_apcor_sort[0]**2)
    err3_avg = np.sqrt(err3_2n**2 + err_avg_apcor_sort[1]**2)
    err4_avg = np.sqrt(err4_2n**2 + err_avg_apcor_sort[2]**2)
    err5_avg = np.sqrt(err5_2n**2 + err_avg_apcor_sort[3]**2)
    err8_avg = np.sqrt(err8_2n**2 + err_avg_apcor_sort[4]**2)

    # If anything was set to 66.666 or 99.999 before, set them back to those values

    # These are sources that fall outside the FOV, both mag and merr are set to 66.666
    fov_m2 = (mag2_2n == 66.666)
    fov_m3 = (mag3_2n == 66.666)
    fov_m4 = (mag4_2n == 66.666)
    fov_m5 = (mag5_2n == 66.666)
    fov_m8 = (mag8_2n == 66.666)

    mag2_avg[fov_m2] = 66.666
    mag3_avg[fov_m3] = 66.666
    mag4_avg[fov_m4] = 66.666
    mag5_avg[fov_m5] = 66.666
    mag8_avg[fov_m8] = 66.666

    err2_avg[fov_m2] = 66.666
    err3_avg[fov_m3] = 66.666
    err4_avg[fov_m4] = 66.666
    err5_avg[fov_m5] = 66.666
    err8_avg[fov_m8] = 66.666

    # Sources with negative fluxes have both mag and merr set to 99.999
    neg_m2 = (mag2_2n == 99.999)
    neg_m3 = (mag3_2n == 99.999)
    neg_m4 = (mag4_2n == 99.999)
    neg_m5 = (mag5_2n == 99.999)
    neg_m8 = (mag8_2n == 99.999)

    mag2_avg[neg_m2] = 99.999
    mag3_avg[neg_m3] = 99.999
    mag4_avg[neg_m4] = 99.999
    mag5_avg[neg_m5] = 99.999
    mag8_avg[neg_m8] = 99.999

    # Sources with very small positive fluxes may have only merr set to 99.999
    tiny_m2 = (err2_2n == 99.999)
    tiny_m3 = (err3_2n == 99.999)
    tiny_m4 = (err4_2n == 99.999)
    tiny_m5 = (err5_2n == 99.999)
    tiny_m8 = (err8_2n == 99.999)

    err2_avg[tiny_m2] = 99.999
    err3_avg[tiny_m3] = 99.999
    err4_avg[tiny_m4] = 99.999
    err5_avg[tiny_m5] = 99.999
    err8_avg[tiny_m8] = 99.999

    # Save to compare to automatic_catalog_target.mag
    zip_table2 = zip(ids, xc5_2n, yc5_2n, ra_2n, dec_2n, mag2_avg, err2_avg, mag3_avg, err3_avg, mag4_avg, err4_avg, mag5_avg, err5_avg, mag8_avg, err8_avg, ci5_2n, nfilt)

    np.savetxt(target_dir + '/photometry/automatic_catalog_' + target + '_avgapcor_add.tab', zip_table2, fmt='%i  %.3f  %.3f  %.8f  %.8f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %i')


#####  Make final centered coordinate file for additional clusters #####

    xc_cen, yc_cen = np.loadtxt(target_dir + '/photometry/center_add.coo', unpack=True)


    # Make a regions file of final clusters
    outputfile = target_dir + '/photometry/catalog_ds9_clusters_final_add.reg'

    file = open(outputfile, 'w')
    file.write('global color=red width=3 font="helvetica 15 normal roman" highlite=1 \n')
    file.write('image\n')

    for i in range(len(xc_cen)):

        x = str(xc_cen[i])
        y = str(yc_cen[i])
        newline = 'circle('+ x + ',' + y +  ',7)\n'
        file.write(newline)

    file.close()



#####  Make CI-based Aperture Correction Final Catalog  #####

    apcorlist = glob.glob(target_dir + '/photometry/apcor_fromci_*_add.mag')

    # Find apcor files corresponding to filter names
    apcorfile_2 = [file for file in apcorlist if 'f275w' in file][0]
    apcorfile_3 = [file for file in apcorlist if 'f336w' in file][0]
    apcorfile_8 = [file for file in apcorlist if 'f814w' in file][0]

    # Special cases
    apcorfile_4 = [file for file in apcorlist if 'f435w' in file or 'f438w' in file][0]
    apcorfile_5 = [file for file in apcorlist if 'f555w' in file or 'f606w' in file][0]

    # Read in CI-based aperture corrections
    apcor2, err_apcor2 = np.loadtxt(apcorfile_2, unpack=True)
    apcor3, err_apcor3 = np.loadtxt(apcorfile_3, unpack=True)
    apcor4, err_apcor4 = np.loadtxt(apcorfile_4, unpack=True)
    apcor5, err_apcor5 = np.loadtxt(apcorfile_5, unpack=True)
    apcor8, err_apcor8 = np.loadtxt(apcorfile_8, unpack=True)

    # Must apply low photometric error selection to CI-based apcors too, otherwise array
    # lengths don't match!
    sel_2n = ((err2_sci <= 0.3) & (err3_sci <= 0.3)) | ((err3_sci <= 0.3) & (err4_sci <= 0.3)) | ((err4_sci <= 0.3) & (err5_sci <= 0.3)) | ((err5_sci <= 0.3) & (err8_sci <= 0.3))


    # Apply CI-based apcors to photometry

    apcor2_2n = apcor2[sel_2n]
    apcor3_2n = apcor3[sel_2n]
    apcor4_2n = apcor4[sel_2n]
    apcor5_2n = apcor5[sel_2n]
    apcor8_2n = apcor8[sel_2n]

    err_apcor2_2n = err_apcor2[sel_2n]
    err_apcor3_2n = err_apcor3[sel_2n]
    err_apcor4_2n = err_apcor4[sel_2n]
    err_apcor5_2n = err_apcor5[sel_2n]
    err_apcor8_2n = err_apcor8[sel_2n]

    mag2_cib = mag2_2n + apcor2_2n - gal_ext[0]
    mag3_cib = mag3_2n + apcor3_2n - gal_ext[1]
    mag4_cib = mag4_2n + apcor4_2n - gal_ext[2]
    mag5_cib = mag5_2n + apcor5_2n - gal_ext[3]
    mag8_cib = mag8_2n + apcor8_2n - gal_ext[4]

    err2_cib = np.sqrt(err2_2n**2 + err_apcor2_2n**2)
    err3_cib = np.sqrt(err3_2n**2 + err_apcor3_2n**2)
    err4_cib = np.sqrt(err4_2n**2 + err_apcor4_2n**2)
    err5_cib = np.sqrt(err5_2n**2 + err_apcor5_2n**2)
    err8_cib = np.sqrt(err8_2n**2 + err_apcor8_2n**2)

    # If magnitudes were 66.666 or 99.999 before, set them back to 66.666 or 99.999

    # These are sources that fall outside the FOV, both mag and merr are set to 66.666
    fov_m2 = (mag2_2n == 66.666)
    fov_m3 = (mag3_2n == 66.666)
    fov_m4 = (mag4_2n == 66.666)
    fov_m5 = (mag5_2n == 66.666)
    fov_m8 = (mag8_2n == 66.666)

    mag2_cib[fov_m2] = 66.666
    mag3_cib[fov_m3] = 66.666
    mag4_cib[fov_m4] = 66.666
    mag5_cib[fov_m5] = 66.666
    mag8_cib[fov_m8] = 66.666

    err2_cib[fov_m2] = 66.666
    err3_cib[fov_m3] = 66.666
    err4_cib[fov_m4] = 66.666
    err5_cib[fov_m5] = 66.666
    err8_cib[fov_m8] = 66.666

    # Do same for sources set to 99.999, but mag and err separately
    # This is because sometimes merr is 99.999, but mag is not
    neg_m2 = (mag2_2n == 99.999)
    neg_m3 = (mag3_2n == 99.999)
    neg_m4 = (mag4_2n == 99.999)
    neg_m5 = (mag5_2n == 99.999)
    neg_m8 = (mag8_2n == 99.999)

    mag2_cib[neg_m2] = 99.999
    mag3_cib[neg_m3] = 99.999
    mag4_cib[neg_m4] = 99.999
    mag5_cib[neg_m5] = 99.999
    mag8_cib[neg_m8] = 99.999

    tiny_m2 = (err2_2n == 99.999)
    tiny_m3 = (err3_2n == 99.999)
    tiny_m4 = (err4_2n == 99.999)
    tiny_m5 = (err5_2n == 99.999)
    tiny_m8 = (err8_2n == 99.999)

    err2_cib[tiny_m2] = 99.999
    err3_cib[tiny_m3] = 99.999
    err4_cib[tiny_m4] = 99.999
    err5_cib[tiny_m5] = 99.999
    err8_cib[tiny_m8] = 99.999

    # If apcorr is 99.999 but photometry is detected (less than 66.666), set to 44.444
    ap_m2 = (apcor2_2n == 99.999) & (mag2_2n < 66.666)
    ap_m3 = (apcor3_2n == 99.999) & (mag3_2n < 66.666)
    ap_m4 = (apcor4_2n == 99.999) & (mag4_2n < 66.666)
    ap_m5 = (apcor5_2n == 99.999) & (mag5_2n < 66.666)
    ap_m8 = (apcor8_2n == 99.999) & (mag8_2n < 66.666)

    mag2_cib[ap_m2] = 44.444
    mag3_cib[ap_m3] = 44.444
    mag4_cib[ap_m4] = 44.444
    mag5_cib[ap_m5] = 44.444
    mag8_cib[ap_m8] = 44.444

    err2_cib[ap_m2] = 44.444
    err3_cib[ap_m3] = 44.444
    err4_cib[ap_m4] = 44.444
    err5_cib[ap_m5] = 44.444
    err8_cib[ap_m8] = 44.444

    # Update number of filters
    nfilt_cib = np.copy(nfilt)

    nfilt_4 = (((mag2_cib < 44.444) & (mag4_cib < 44.444) & (mag5_cib < 44.444) & (mag8_cib < 44.444)) | ((mag3_cib < 44.444) & (mag4_cib < 44.444) & (mag5_cib < 44.444) & (mag8_cib < 44.444))) & ((nfilt_cib == 4) | (nfilt_cib == 5))

    nfilt_5 = ((mag2_cib < 44.444) & (mag3_cib < 44.444) & (mag4_cib < 44.444) & (mag5_cib < 44.444) & (mag8_cib < 44.444)) & (nfilt_cib == 5)

    nfilt_cib[nfilt_4] = 4
    nfilt_cib[nfilt_5] = 5
    nfilt_cib[~nfilt_5 & ~nfilt_4] = 0

    # Save to compare to automatic_catalog_target.mag
    zip_table3 = zip(ids, xc5_2n, yc5_2n, ra_2n, dec_2n, mag2_cib, err2_cib, mag3_cib, err3_cib, mag4_cib, err4_cib, mag5_cib, err5_cib, mag8_cib, err8_cib, ci5_2n, nfilt_cib)

    np.savetxt(target_dir + '/photometry/automatic_catalog_' + target + '_cibased_add.tab', zip_table3, fmt='%i  %.3f  %.3f  %.8f  %.8f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %i')


    # move final files to final directory
    source = target_dir + '/photometry/automatic_catalog_' + target + '_avgapcor_add.tab'
    destination = target_dir + '/final/automatic_catalog_' + target + '_avgapcor_add.tab'
    os.rename(source, destination)

    source = target_dir + '/photometry/automatic_catalog_' + target + '_cibased_add.tab'
    destination = target_dir + '/final/automatic_catalog_' + target + '_cibased_add.tab'
    os.rename(source, destination)

    source = target_dir + '/photometry/center_add.coo'
    destination = target_dir + '/final/cat_UBVI_' + target + '_center_add.coo'
    os.rename(source, destination)

    source = target_dir + '/photometry/catalog_ds9_clusters_final_add.reg'
    destination = target_dir + '/final/catalog_ds9_clusters_final_add.reg'
    os.rename(source, destination)


    print ''
    print 'Check final catalogs in /final directory.'
    print ''

    sys.exit('STEP 6 is finished')


sys.exit('fin')
