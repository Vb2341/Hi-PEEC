#!/usr/bin/env python
#-*- coding: utf-8 -*-

#------------------------------------------------------------------------------
#Title: Hi-PEEC Extraction and photometry
#Author:        Axel Runnholm
#Creation date: 2016-06-16
#Description:   This script does the source extraction and photometry for the
#               Hi-PEEC project.
#------------------------------------------------------------------------------


#import packages
#------------------------------------------------------------------------------
# compensation for main diffs between Python 2.7 and 3
from __future__ import division#,print_function

#import math, datahandling and plotting utils
import numpy as np
import math
import matplotlib.pyplot as plt

#sys utils
import os, glob
import time
import shutil
import sys
import string
import ast

#astronomy utils
import pyfits
from pyraf import iraf
import pywcs

#Import Hi-PEEC modules
sys.path.insert(0, './source/')
import filemanagement
#------------------------------------------------------------------------------

iraf.noao(_doprint=0)
iraf.digiphot(_doprint=0)
iraf.obsutil(_doprint=0)
iraf.daophot(_doprint=0)


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
iraf.phot.verbose = 'no'
iraf.phot.verify = 'no'
#------------------------------------------------------------------------------



def extraction(userinputs):
    #Set up required variables
    target_dir = userinputs['OUTDIR']
    seximage = userinputs['IMAGE']

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

    with open(outputfile, 'w') as file:
        file.write('global color=blue width=5 font="helvetica 15 normal roman" highlite=1 \n')
        file.write('image\n')

        for i in range(len(xx)):
            newline = 'circle(' + str(xx[i]) + ',' + str(yy[i]) +  ',7) \n'
            file.write(newline)
    print ''
    print 'Check catalog_ds9_sextractor.reg in the /s_extraction directory for'
    print 'the quality of source extraction.'
    print ''

    return target_dir + '/s_extraction/R2_wl_dpop_detarea.cat'

def photometry(userinputs, image, catalog, outputname, apertures):
    """Function for performing IRAF photometry
    @Params:
    userinputs  (dict)  - dictionary with results from the user input file.
    image       (STR)   - fitsfile we want to do photometry on
    catalog     (STR)   - input coordinates where we want to do photometry
    outputname  (STR)   - name of the file where we store the measured results
    apertures   (STR)   - what apertures to measure. Should be a string i.e. '1.0,3.0'
    """
    #set directory
    target_dir = userinputs['OUTDIR']

    #Update passed names  to be full paths if they are not

    if len(image.split('/'))==1:
        image = glob.glob(target_dir + '/img/' + image)
        if len(image)==0:
            sys.exit('Selected image does not exist')
        else:
            image = image[0]

    if len(catalog.split('/'))==1:
        catalog = target_dir + '/init/' + catalog

    if len(outputname.split('/'))==1:
        output = target_dir + '/photometry/' + outputname
    else:
        output = outputname
        outputname = outputname.split('/')[-1]


    #Load zeropoints
    inst_zp, filter_zp, zp_zp = np.loadtxt(target_dir + '/init/legus_zeropoints.tab', unpack=True, dtype='str')

    # Get filter from header
    filter = pyfits.getheader(image)['FILTER']


    # Set the necessary variables for photometry on the reference image
    exptime = pyfits.getheader(image)['EXPTIME']
    inst = pyfits.getheader(image)['INSTRUME']
    inst = inst.lower()

    match = (inst_zp == inst) & (filter_zp == filter.lower())
    zp = zp_zp[match]

    # zp is a string within an array, so need to turn into a float
    try:
        zp = float(zp[0])
        #If that cannot be done there was no match.
    except IndexError:
        sys.exit('No zeropoint was found for filter: {}'.format(filter))

    # Remove output file if it already exists
    if os.path.exists(output) == True:
        os.remove(output)

    # Run photometry
    iraf.datapars.epadu = exptime

    iraf.centerpars.calgorithm = 'centroid'

    iraf.fitskypars.annulus = userinputs['ANNULUS']
    iraf.fitskypars.dannulu = userinputs['D_ANNULUS']

    iraf.photpars.apertures = apertures
    iraf.photpars.zmag = zp

    iraf.phot(image, catalog, output)



    #Depending on the number of apertures used, different methods of saving the
    # results are required
    #--------------------------------------------------------------------------

    naper = len(apertures.split(','))

    #final output filename
    fullcat_mag_short = target_dir + '/photometry/short_' + outputname

    if naper > 1:
        # Removes all outputlines that do not contain the character '*'
        # ensures only phot results are kept
        cmd = 'grep "*" ' + output + ' > ' + fullcat_mag_short
        os.system(cmd)

        # Replace INDEFS:
        cmd = 'sed -i.bak "s/INDEF/99.999/g" ' + fullcat_mag_short
        os.system(cmd)

        dumpfile = target_dir+'/photometry/dumptest.mag'

        # Remove .bak files to prevent confusion
        bak_fullcat = fullcat_mag_short + '.bak'
        os.remove(bak_fullcat)


    else:
        #Dump results into a temp file
        temp = target_dir + '/photometry/phot_dump.mag'
        filemanagement.remove_if_exists(temp)
        iraf.txdump(output, 'XCENTER,YCENTER,FLUX,MAG,MERR,MSKY,ID', 'yes', Stdout = temp)

        # Set placeholders for sources outside of FOV and undetected sources
        # For outside of FOV, use 66.666 instead of INDEF
        # For undetected sources, use 99.999 instead of INDEF

        # Sources outside of FOV have exactly zero flux
        x, y, flux, mag, merr, msky, id = np.loadtxt(temp, unpack = True,
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

        zip_phot = zip(x, y, mag, merr, msky, id)

        np.savetxt(fullcat_mag_short, zip_phot,
                    fmt = '%.3f  %.3f  %.3f  %.3f  %.9f  %i')

    #--------------------------------------------------------------------------

    return fullcat_mag_short


def growth_curve(userinputs, catalog):

    #Load the photometry results from the catalog (that is returned by the phot
    #function)
    aper_st, flux_st = np.loadtxt(catalog, unpack=True, usecols=(0,3))

    #Growth curve is only done on the ref image so we get the filter from userinp.
    ref_filter = userinputs['REF_FILTER']

    ratio_st = np.empty(len(aper_st))

    #number of apertures
    naper = 20

    # Calculate the number of stars, make sure it is an integer
    nstar = int(len(aper_st)/naper)
    aper_ind = naper - 1

    for k in range(nstar):

        for i in range(naper):

            ratio_st[i + k*naper] = flux_st[i + k*naper]/flux_st[aper_ind + k*naper]



    # Find median ratio at each aperture between all the stars and all the clusters
    med_st = np.empty(naper)

    for i in range(naper):

        med_st[i] = np.median(ratio_st[i::naper])


    # Plot growth curves
    fig = plt.figure(figsize = (7,7))

    aper_x = np.arange(naper) + 1

    for i in range(nstar):

        ratio_y = ratio_st[i*naper:(i + 1)*naper]
        plt.plot(aper_x, ratio_y, 'y-')
        plt.annotate(str(i + 1), xy=(8.0, ratio_y[7]),
            horizontalalignment='left', verticalalignment='top', fontsize=6)


    plt.plot(aper_x, med_st, 'r-' , linewidth=4.0)
    plt.hlines(0.5, 0, 20, color='black', linewidth=2, zorder=10)
    plt.vlines(4, 0, 1.1, color='black', linewidth=2, linestyle='dashed', zorder=10)
    plt.vlines(5, 0, 1.1, color='black', linewidth=2, linestyle='dashed', zorder=10)
    plt.vlines(6, 0, 1.1, color='black', linewidth=2, linestyle='dashed', zorder=10)

    plt.ylabel('Normalized Flux ' + ref_filter.upper())
    plt.xlabel('Radius (pix)')
    plt.xlim(1,20)
    plt.minorticks_on()

    fig.savefig(userinputs['OUTDIR'] + '/plots/plot_growth_curve.pdf')

