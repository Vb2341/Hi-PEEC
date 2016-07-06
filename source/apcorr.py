#!/usr/bin/env python
#-*- coding: utf-8 -*-

#------------------------------------------------------------------------------
#Title: Hi-PEEC Aperture correction
#Author:        Axel Runnholm
#Creation date: 2016-06-20
#Description:   This script contains the routines for estimating and applying
#               aperture corrections.

#------------------------------------------------------------------------------


#import packages
#------------------------------------------------------------------------------
# compensation for main diffs between Python 2.7 and 3
from __future__ import division#,print_function

#import math, data handling and plotting utils
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
import extraction
#-------------------------------------------------------------------------------

def calculation(userinputs):
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

    # Set aperture sizes
    small_ap = 4
    big_ap = 20

    #Set directory of the photometry as variable since it is target for the function
    phot_dir = userinputs['OUTDIR'] + '/photometry/'

    #Get the list of images
    imagelist = glob.glob(userinputs['OUTDIR'] +'/img/*sci*.fits')

    #Clear old apcorr file
    apcorrfile = userinputs['OUTDIR'] + '/photometry/avg_aperture_correction.txt'
    filemanagement.remove_if_exists(apcorrfile)

    logging.info('Doing apcorr photometry')
    for image in imagelist:
        filter = extraction.get_filter(image)

        #-----------------------------------------------------------------------
        # Do required photometry
        #-----------------------------------------------------------------------

        # Large (20px) aperture photometry
        photometry_file_large = phot_dir + 'large_apcorr_' + filter + '.mag'

        apcorr_cat_large = extraction.photometry(userinputs, image,
                            userinputs['STARS'], photometry_file_large,
                            '20.0', annulus=21.0, dannulus=1.0)

        # Small (user selected) aperture photometry
        photometry_file_small = phot_dir + 'small_apcorr_' + filter + '.mag'

        apcorr_cat_small = extraction.photometry(userinputs, image,
                            userinputs['STARS'], photometry_file_small,
                            str(userinputs['AP_RAD']))


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
        try:
            apcor_avg = np.mean(apcor[lim])
            apcor_err = np.std(apcor_lim)/np.sqrt(len(apcor_lim))
        except RuntimeWarning:
            loggin.critical('No stars in the apcorr star list after filtering. Quitting')
            loggin.debug('Check {} and the aperture correction requirements'\
                         .format(userinputs['STARS']))
            sys.exit('No stars left after filtering. Check {} and the aperture \
                correction requirements'.format(userinputs['STARS']))

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
