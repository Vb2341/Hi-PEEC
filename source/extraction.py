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
import logging

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

def get_filter(image):
    """
    Returns the filter from a fits file by searching the header
    @params
    image (STR)     - path to image file

    @returns
    filter (STR)    - string with filter name
    """
    logging.debug('Retrieving filter for {}'.format(image))
    try:
        filter = pyfits.getheader(image)['FILTER']
    except KeyError:
        logging.debug('No FILTER key found, trying FILTER2')
        #The 814 image has the filter information under the keyword FILTER2:
        filter = pyfits.getheader(image)['FILTER2']

    return filter

def ACS_zeropoint(image):
    """
    Calculates the zeropoint for ACS images from the header info and returns it
    @params
    image (STR)     - path to image file

    @returns
    zeropoint (FLOAT)    - magnitude zeropoint for frame
    """
    logging.debug('Calculating zeropoint for {}'.format(image))

    PHOTFLAM = pyfits.getheader(image)['PHOTFLAM']
    PHOTPLAM = pyfits.getheader(image)['PHOTPLAM']

    ABMAG_ZEROPOINT=-2.5*np.log10(PHOTFLAM)-5*np.log10(PHOTPLAM)-2.408

    return ABMAG_ZEROPOINT

def extraction(userinputs):
    #Set up required variables
    target_dir = userinputs['OUTDIR']
    seximage = userinputs['IMAGE']
    logging.info('Running sextractor on {}'.format(userinputs['IMAGE']))

    print 'Executing SExtractor on user selected image : ', seximage

    # Verify that file exists
    if os.path.exists(userinputs['DATA'] + '/' + seximage) == False:
        print 'File ' + seximage + ' could not be found in ' + userinputs['DATA']
        logging.critical(' Could not find {}. Quitting'.format(seximage))
        logging.debug('Looking for {} but unable to locate'.format(userinputs['DATA'] + '/' + seximage))
        filemanagement.shutdown('Quitting now...',userinputs)

    # Run sextractor
    logging.info('Start sextractor')
    os.chdir(target_dir + '/s_extraction')
    logging.debug('Changed dir to {}'.format(os.getcwd()))
    command = 'sex ' + userinputs['DATA'] + '/' + seximage + '  -c R2_wl_aa.config'
    os.system(command)
    os.chdir(target_dir)
    logging.debug('Changed working directory back to {}'.format(target_dir))

    # Read in results and make regions file of objects sextracted
    logging.info('Read in Sextractor catalog')
    xx, yy = np.loadtxt(target_dir + '/s_extraction/R2_wl_dpop_detarea.cat', unpack=True,
                        skiprows=5, usecols=(0,1))

    outputfile = target_dir + '/s_extraction/catalog_ds9_sextractor.reg'

    logging.info('Writing region file from source extractor data')
    logging.debug('Sextractor file: {}'.format(target_dir + '/s_extraction/R2_wl_dpop_detarea.cat'))
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

def photometry(userinputs, image, catalog, outputname, apertures, annulus='', dannulus=''):
    """Function for performing IRAF photometry
    @Params:
    userinputs  (dict)  - dictionary with results from the user input file.
    image       (STR)   - fitsfile we want to do photometry on
    catalog     (STR)   - input coordinates where we want to do photometry
    outputname  (STR)   - name of the file where we store the measured results
    apertures   (STR)   - what apertures to measure. Should be a string i.e. '1.0,3.0'
    annulus     (FLOAT) - (optional) which skyannulus to use, if not set the one defined in
                          user inputs is used
    dannulus    (FLOAT) - (optional) which diameter to use for the sky annulus

    @Returns
    output      (STR)   - full path to the final catalog file
    """
    logging.info('Running photometry function on {}'.format(image))
    logging.info('Using {}px apertures'.format(apertures))

    #set directory
    target_dir = userinputs['OUTDIR']

    #Update passed names  to be full paths if they are not

    if len(image.split('/'))==1:
        logging.info('Looking for {} in {}.'.format(image,userinputs['DATA']))
        image = glob.glob(userinputs['DATA'] + '/' + image)
        if len(image)==0:
            logging.critical('No {} image found'.format(image))
            filemanagement.shutdown('Selected image does not exist',userinputs)
        else:
            image = image[0]
    logging.debug('Using image: {}'.format(image))

    if len(catalog.split('/'))==1:
        catalog = target_dir + '/init/' + catalog
        logging.debug('Input catalog: {}'.format(catalog))

    if len(outputname.split('/'))==1:
        output = target_dir + '/photometry/' + outputname
        logging.debug('Output name: {}'.format(output))
    else:
        output = outputname
        outputname = outputname.split('/')[-1]
        logging.debug('Output name: {}'.format(output))


    #Load zeropoints
    inst_zp, filter_zp, zp_zp = np.loadtxt(target_dir + '/init/Hi-PEEC_zeropoints.tab', unpack=True, dtype='str')

    # Get filter from header
    filter = get_filter(image)


    # Set the necessary variables for photometry on the reference image
    exptime = pyfits.getheader(image)['EXPTIME']
    logging.debug('Exposure time from header: {}'.format(exptime))
    inst = pyfits.getheader(image)['INSTRUME']
    logging.debug('Intrument from header: {}'.format(inst))
    inst = inst.lower()

    match = (inst_zp == inst) & (filter_zp == filter.lower())
    zp = zp_zp[match]

    # zp is a string within an array, so need to turn into a float
    try:
        zp = float(zp[0])
        #If that cannot be done there was no match.
    except IndexError:
        if inst == 'acs':
            logging.debug('Zeropoint not found in file, passing to ACS calculation')
            zp = ACS_zeropoint(image)
        else:
            logging.critical('No matching zeropoint found. Quitting.')
            logging.debug('No zeropoint match found for filter {} with instrument {}'\
                          .format(filter,inst))
            logging.debug('Available filters in zeropoint file : {} for instrument {}'\
                          .format(filter_zp, inst_zp))
            filemanagement.shutdown('No zeropoint was found for filter: {}'.format(filter),userinputs)

    logging.debug('Zeropoint from file: {}'.format(zp))
    # Remove output file if it already exists
    filemanagement.remove_if_exists(output)


    # Run photometry
    #--------------------------------------------------------------------------
    # Set up IRAF params:
    iraf.datapars.epadu = exptime

    iraf.centerpars.calgorithm = 'centroid'

    # Select the annulus depending on whether it is overwritten in the function call or not
    if annulus == '':
        iraf.fitskypars.annulus = userinputs['ANNULUS']
        logging.info('Using annulus from inputfile ({}px)'.format(userinputs['ANNULUS']))
    else:
        iraf.fitskypars.annulus = annulus
        logging.info('Using user specified annulus ({}px)'.format(annulus))
    if dannulus == '':
        iraf.fitskypars.dannulu = userinputs['D_ANNULUS']
        logging.info('Using annulus width from inputfile ({}px)'.format(userinputs['D_ANNULUS']))
    else:
        iraf.fitskypars.dannulu = dannulus
        logging.info('Using user specified annulus width ({}px)'.format(dannulus))

    iraf.photpars.apertures = apertures
    iraf.photpars.zmag = zp

    # Do phot
    iraf.phot(image, catalog, output)
    #--------------------------------------------------------------------------


    #Depending on the number of apertures used, different methods of saving the
    # results are required
    #--------------------------------------------------------------------------

    naper = len(apertures.split(','))
    logging.debug('Number of apertures used {}'.format(naper))

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
        logging.debug('Number of sources outside FOV: {}'.format(len(out_fov)))

        mag[out_fov] = 66.666
        merr[out_fov] = 66.666
        msky[out_fov] = 66.666

        # Undetected sources, those with negative flux or fluxes so small that mag err
        # is INDEF
        neg_flux = (flux < 0.)
        tiny_flux = (flux > 0.) & (merr == 'INDEF')

        mag[neg_flux] = 99.999
        merr[neg_flux] = 99.999
        msky[neg_flux] = 99.999

        merr[tiny_flux] = 99.999
        msky[tiny_flux] = 99.999

        logging.debug('Nr of undetected sources: {}'.format(len(tiny_flux)+len(neg_flux)))
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
    logging.info('Running growth curve analysis on {}'.format(catalog))
    # Load the photometry results from the catalog (that is returned by the phot
    # function)
    aper_st, flux_st = np.loadtxt(catalog, unpack=True, usecols=(0,3))

    #Growth curve is only done on the ref image so we get the filter from userinp.
    ref_filter = userinputs['REF_FILTER']

    ratio_st = np.empty(len(aper_st))

    #number of apertures
    naper = 20

    # Calculate the number of stars, make sure it is an integer
    nstar = int(len(aper_st)/naper)
    logging.info('Number of stars used: {}'.format(nstar))
    aper_ind = naper - 1

    for k in range(nstar):

        for i in range(naper):

            ratio_st[i + k*naper] = flux_st[i + k*naper]/flux_st[aper_ind + k*naper]


    # Find median ratio at each aperture between all the stars and all the clusters
    med_st = np.empty(naper)

    for i in range(naper):

        med_st[i] = np.median(ratio_st[i::naper])


    # Plot growth curves
    logging.info('Creating Growth curve plots')
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

