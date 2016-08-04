#!/usr/bin/env python
#-*- coding: utf-8 -*-

#------------------------------------------------------------------------------
#Title: Hi-PEEC Catalogue management
#Author:        Axel Runnholm
#Creation date: 2016-06-20
#Description:   This script contains the routines for estimating and applying
#               aperture corrections, extinction corrections, and catalogue
#               filtering.

#------------------------------------------------------------------------------


#import packages
#------------------------------------------------------------------------------
# compensation for main diffs between Python 2.7 and 3
from __future__ import division#,print_function

# import math, data handling and plotting utils
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd

# sys utils
import os, glob
import time
import shutil
import sys
import string
import ast
import logging
import warnings

# astronomy utils
import pyfits
from pyraf import iraf
import pywcs

# Import Hi-PEEC modules
sys.path.insert(0, './source/')
import filemanagement
import extraction
#-------------------------------------------------------------------------------

def apcorr_calc(userinputs):
    """Function for performing IRAF photometry
    @Params:
    userinputs  (dict)  - dictionary with results from the user input file.

    @Returns
    output      (STR)   - full path to the aperture corrections file
    """
    #---------------------------------------------------------------------------
    # Setup stuff
    #---------------------------------------------------------------------------
    # Read upper and lower limits for values to include in apcor average
    lowlim = userinputs['LOWLIM']
    uplim = userinputs['UPLIM']

    # Check if there is an additional single star file.
    try:
        # Split up the input
        add_stars = userinputs['ADD_STARS']
        add_stars_list = add_stars.split(',')
        if len(add_stars_list) == 1:
            add_stars_list = [add_stars]
        add_stars_list = add_stars_list

        add_filter = userinputs['ADD_FILTER']
        add_filter_list = add_filter.split(',')

        add_filter_list = [a.lower() for a in add_filter_list]

        logging.info('Using additional starfiles {} for filters {}'
                     .format(add_stars_list,add_filter_list))

        if len(add_filter_list) != len(add_stars_list):
            if len(add_stars_list) != 1:
                logging.info('Mismatch in number of additional filters and files.')
                filemanagement.shutdown('Number of filters and additional starfiles does not match.')
        additional_starfiles = True
    except KeyError:
        additional_starfiles = False

    #Set directory of the photometry as variable since it is target for the function
    phot_dir = userinputs['OUTDIR'] + '/photometry/'

    #Get the list of images
    imagelist = glob.glob(userinputs['DATA'] + '/*sci*.fits')

    #Clear old apcorr file
    apcorrfile = userinputs['OUTDIR'] + '/photometry/avg_aperture_correction.txt'
    filemanagement.remove_if_exists(apcorrfile)

    logging.info('Doing apcorr photometry')
    for image in imagelist:
        filter = extraction.get_filter(image)

        # Choose input parameters
        if additional_starfiles:
            if filter.lower() in add_filter_list:
                if len(add_stars_list) != 1:
                    starfile = [add_stars_list[a] for a in range(len(add_filter_list))\
                                if add_filter_list[a] == filter.lower()][0]
                else:
                    starfile = add_stars_list[0]
            else:
                starfile = userinputs['STARS']
        else:
            starfile = userinputs['STARS']

        print filter
        print starfile
        logging.info('Using {} as star input for apcorr'.format(starfile))

        #-----------------------------------------------------------------------
        # Do required photometry
        #-----------------------------------------------------------------------

        # Large (20px) aperture photometry
        photometry_file_large = phot_dir + 'large_apcorr_' + filter + '.mag'

        apcorr_cat_large = extraction.photometry(userinputs, image,
                            starfile, photometry_file_large,
                            '20.0', annulus=21.0, dannulus=1.0)

        # Small (user selected) aperture photometry
        photometry_file_small = phot_dir + 'small_apcorr_' + filter + '.mag'

        # calculate appropriate aperture to use
        ap = extraction.calc_aperture(userinputs,image)
        apcorr_cat_small = extraction.photometry(userinputs, image,
                            starfile, photometry_file_small,
                            ap)


        #-----------------------------------------------------------------------
        # Calculate apcorrs
        #-----------------------------------------------------------------------

        # Load the photometry
        x, y, mag_small = np.loadtxt(apcorr_cat_small, unpack=True, usecols=(0,1,2))
        mag_large = np.loadtxt(apcorr_cat_large, unpack=True, usecols=(2,))

        # Calculate aperture corrections
        apcor = mag_large - mag_small

        # Limit range of aperture corrections allowed to go into average
        lim = (apcor < uplim) & (apcor > lowlim)
        apcor_lim = apcor[lim]
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                apcor_avg = np.mean(apcor[lim])
                apcor_err = np.std(apcor_lim)/np.sqrt(len(apcor_lim))
                logging.info('Nr of stars in {} apcorr mean: {}'.format(filter,len(apcor[lim])))
            except RuntimeWarning:
                logging.warning('No stars in the apcorr star list for {} after filtering.'.format(filter))
                logging.debug('Check {} and the aperture correction requirements'\
                             .format(userinputs['STARS']))
                print ''
                print '\t WARNING:'
                print '\t No stars in the apcorr star list for {} after filtering. Check {} and the aperture correction requirements'\
                    .format(filter, userinputs['STARS'])
                apcor_avg = 0
                apcor_err = 0

        #Save these results to file
        with open(apcorrfile,'a') as file:
            avg = str(apcor_avg)
            err = str(apcor_err)
            lim = str(len(apcor_lim))
            file.write(filter + '\t' + avg + '\t' + err + '\t' + lim + '\n')
        print ''
        print '\t Filter: ' + filter
        print '\t Apcor: %.3f' % apcor_avg
        print '\t Apcor error: %.3f' % apcor_err
        logging.info('{}: Apcorr = {} +- {}'.format(filter, apcor_avg, apcor_err))

        # Make a plot for each filter
        fig = plt.figure(figsize = (7,7))

        # Draw histogram
        logging.info('Create aperture correction plots')
        num, bins, patches = plt.hist(apcor, bins=50, range=(-4.0,1.0),
                                      color='red', lw=2, histtype='step')

        plt.vlines(apcor_avg, 0, 1, color='blue', linewidth=5)
        plt.vlines(uplim, 0, 2, color='black', linestyle='--', linewidth=2)
        plt.vlines(lowlim, 0, 2, color='black', linestyle='--', linewidth=2)

        plt.xlabel('Aperture Correction')
        plt.ylabel('N')
        plt.minorticks_on()
        plt.axis([-4.0,1.0, 0, max(num) + 0.1*max(num)])

        fig.savefig(userinputs['OUTDIR'] + '/plots/apcor_' + filter + '.pdf')

    return apcorrfile



def BVI_detection(cat, filters, userinput):
    """
    Function for filtering out all sources not detected in the B,V and I filters

    Inputs:
        cat (pandas DATAFRAME) - dataframe containing the cluster data
        filters (LIST) - List of filters that are in the catalog
        userinput (DICTIONARY) - dictionary containing user input parameters
    Output:
        cat (PANDAS DATAFRAME) - dataframe containing the cluster data
    """
    if 'F438W' in filters:
        B = 'F438W'
    else:
        B = 'F435W'

    if 'F555W' in filters:
        V = 'F555W'
    else:
        V = 'F606W'

    I = 'F814W'
    maxerr = userinput['MAGERR']

    BVImags = ['Mag ' + s for s in filters if s in [B,V,I]]
    BVIerrs = ['Err ' + s for s in filters if s in [B,V,I]]
    logging.debug('selecting sources detected in {}'.format(BVImags))
    try:
        selection = ((cat[BVIerrs[0]]<=maxerr)
                 & (cat[BVIerrs[1]]<=maxerr)
                 & (cat[BVIerrs[2]]<=maxerr))
    except IndexError:
        logging.critical('F336W,F555W, and F814W are not all present. Unable to do selection')
        logging.debug('Filters present:{}'.format(filters))
        filemanagement.shutdown('F336W,F555W, and F814W are not all present. Unable to do selection. Shutting down',userinput)

    len_before = len(cat['X'])
    cat = cat.drop(cat[~selection].index)
    len_after = len(cat['X'])

    logging.info('Removed all sources not detected in BVI filters.\
    Nr of sources removed {}'.format(len_before-len_after))

    return cat



def remove_duplicates(cat, tolerance=0.5):
    """
    Function for removing all the duplicate sources from the catalog
    Inputs:
        cat (pandas DATAFRAME) - dataframe containing the cluster data
        tolerance - tolerance to use for establishing whether there is a
                    duplicate detection

    Output:
        cat (PANDAS DATAFRAME) - dataframe containing the cluster data
    """


    #Create boolean array
    for_removal = np.zeros(len(cat['X']), dtype = bool)

    prelength =len(cat['X'])

    x = cat['X'].as_matrix()
    y = cat['Y'].as_matrix()

    #Loop over sources to check if there are duplicates
    for k in np.arange(len(x - 1)):

        # Calculate distance between kth object and subsequent objects
        d = np.sqrt((x[k] - x[k + 1:])**2 + (y[k] - y[k + 1:])**2)

        # If any of the distances calculated is less than the tolerance, change removal status
        if (d < tolerance).any():
            for_removal[k] = True

    #Drop the sources found in the loop
    cat = cat.drop(cat[for_removal].index)
    postlength = len(cat['X'])

    nr_of_duplicates = prelength - postlength

    print '\t\t Duplicates removed: {}'.format(nr_of_duplicates)
    logging.info('{} duplicate sources removed'.format(nr_of_duplicates))

    return cat



def apply_corrs(userinput, cat):
    """
    Function for applying extinction and aperture corrections to the catalog
    Inputs:
        userinput (DICTIONARY) - dictionary containing user input parameters
        cat (pandas DATAFRAME) - dataframe containing the cluster data

    Output:
        cat (PANDAS DATAFRAME) - dataframe containing the cluster data
    """
    # Assign the galactic extinction for the five instrument/filter combinations
    # Read in extinction file as array (a vs u means acs/uvis)
    target_dir = userinput['OUTDIR']
    extfile = target_dir + '/init/Hi-PEEC_galactic_extinction.tab'
    extinctions = pd.read_table(extfile, header=[0,1],index_col=0, sep=r"\s*")

    # Pick extinctions for our target only
    target = userinput['TARGET'].lower()

    extinctions = extinctions.loc[['ngc1614']]

    imlist = glob.glob(userinput['DATA'] + '/*sci*.fits')


    # Load apcorrs
    logging.info('Loading apcorrs from file')
    apf, apc, ape = np.loadtxt(target_dir + '/photometry/avg_aperture_correction.txt',
                                  unpack=True, usecols=(0, 1, 2),dtype='str')
    apc = apc.astype(float)
    ape = ape.astype(float)

    for image in imlist:

        # Get filter from filename
        filter = extraction.get_filter(image)

        # Get instrument from header
        instr = pyfits.getheader(image)['INSTRUME']
        instr = instr.lower()

        # Select the correct extinction
        gal_ext = extinctions[filter.lower()][instr]

        # Select the appropriate mag column
        mag = cat['Mag ' + filter].as_matrix()
        err = cat['Err ' + filter].as_matrix()

        # Select the apcor and apcor error corresponding to the current filter
        ap_corr = apc[apf == filter]
        ap_err = ape[apf == filter]

        logging.debug('Corrections for {}: aperture:{}, extrinction: {}'.format(filter,ap_corr,float(gal_ext)))
        # Apply extinctions and apcorrs
        for a in range(len(mag)):
            if mag[a]!=99.999 and mag[a]!=66.666:
                mag[a] = mag[a] + ap_corr - gal_ext
                err[a] = np.sqrt(err[a]**2 + ap_err**2)

        # Insert the new mags into the catalog
        cat['Mag ' + filter] = mag
        cat['Err ' + filter] = err

    return cat

def insert_WCS(userinput, cat):
    """
    Function for inserting WCS coordinates into the catalog
    Inputs:
        cat (pandas DATAFRAME) - dataframe containing the cluster data
    Output:
        cat (PANDAS DATAFRAME) - dataframe containing the cluster data
    """
    target_dir = userinput['OUTDIR']
    imlist = glob.glob(userinput['DATA'] + '/*sci*.fits')

    # Convert xy coordinates into RA Dec of reference filter
    ref_image = [image for image in imlist if userinput['REF_FILTER'].lower() in image.lower()][0]

    # Get header from reference image using pyfits
    header_ref = pyfits.getheader(ref_image)

    # Get wcs solution from reference image header
    logging.info('Get and insert WCS coordinates for each source')
    wcs_ref = pywcs.WCS(header_ref)

    # Calculate RA and Dec for xy coordinates. 1 refers to origin of image in ds9.
    ra, dec = wcs_ref.wcs_pix2sky(cat['X'], cat['Y'], 1)
    cat.insert(2,'RA',ra)
    cat.insert(3,'DEC',dec)

    return cat
