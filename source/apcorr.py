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
#------------------------------------------------------------------------------

def calculation(userinputs):
    #------------------------------------------------------------------------------
    # Setup stuff
    #------------------------------------------------------------------------------

    small_ap = 4
    big_ap = 20

    #Set directory of the photometry as variable since it is target for the function
    phot_dir = userinputs['OUTDIR'] + 'photometry/

    #Get the list of images
    imagelist = glob.glob(target_dir+'/img/*sci*.fits')

    #Clear old apcorr file
    apcorrfile = target_dir + '/photometry/avg_aperture_correction.txt'
    filemanagement.remove_if_exists(apcorrfile)

    #------------------------------------------------------------------------------
    # Do required photometry
    #------------------------------------------------------------------------------


    for image in imagelist:
        try:
            filter = pyfits.getheader(image)['FILTER']
        except KeyError:
            #The 814 image has the filter information under the keyword FILTER2:
            filter = pyfits.getheader(image)['FILTER2']

        out_photometry_file = phot_dir + '20px_apcorr_' + filter + '.mag'

        extraction_catalog = extraction.photometry(userinputs, image, userinputs['STARS'], out_photometry_file, '20.0')





