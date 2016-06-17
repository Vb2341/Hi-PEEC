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

def photometry(userinputs):
    pass
