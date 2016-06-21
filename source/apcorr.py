#!/usr/bin/env python
#-*- coding: utf-8 -*-

#------------------------------------------------------------------------------
#Title: Hi-PEEC Aperture correction
#Author:        Axel Runnholm
#Creation date: 2016-06-16
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
    phot_dir = userinputs['OUTDIR'] + 'photometry/

    #Get the list of images
    imagelist = glob.glob(target_dir+'/img/*sci*.fits')

    #Clear old apcorr file
    apcorrfile = target_dir + '/photometry/avg_aperture_correction.txt'
    filemanagement.remove_if_exists(apcorrfile)


    for image in imagelist:
        try:
            filter = pyfits.getheader(image)['FILTER']
        except KeyError:
            #The 814 image has the filter information under the keyword FILTER2:
            filter = pyfits.getheader(image)['FILTER2']

        #-----------------------------------------------------------------------
        # Do required photometry
        #-----------------------------------------------------------------------

        # 20 px aperture photometry
        photometry_file_20 = phot_dir + '20px_apcorr_' + filter + '.mag'

        apcorr_cat_20 = extraction.photometry(userinputs, image,
                            userinputs['STARS'], photometry_file_20,
                            '20.0', annulus=21.0, dannulus=1.0)

        # 4px aperture photometry 
        photometry_file_4 = phot_dir + '4px_apcorr_' + filter + '.mag'

        apcorr_cat_4 = extraction.photometry(userinputs, image,
                            userinputs['STARS'], photometry_file_20,
                            '4.0')


        #-----------------------------------------------------------------------
        # Calculate individual appcorrs
        #-----------------------------------------------------------------------
       
        # Load the photometry
        x, y, mag_4 = np.loadtxt(apcorr_cat_4, unpack=True, usecols=(0,1,2))
        mag_20 = np.loadtxt(apcorr_cat_20, unpack=True, usecols=(2,))

        # Calculate aperture corrections
        apcor = mag_20 - mag_4

        # Limit range of aperture corrections allowed to go into average
        lim = (apcor < uplim) & (apcor > lowlim)
        apcor_lim = apcor[lim]
        apcor_avg = np.mean(apcor[lim])
        apcor_err = np.std(apcor_lim)/np.sqrt(len(apcor_lim))

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
