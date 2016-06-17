#!/usr/bin/env python
#-*- coding: utf-8 -*-

#------------------------------------------------------------------------------
#Title: Hi-PEEC Python Pipeline
#Author:        Axel Runnholm
#Creation date: 2016-06-16
#Description:   This is the cluster extraction pipeline for the Hi-PEEC project
#               (PI: Angela Adamo). For information on setup and running of the
#               file see the manual and the README included in this folder.
#
#               The code is based on the LEGUS cluster extraction pipeline
#               (see Calzetti(2015) and Adamo(in prep) )
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
import extraction
#------------------------------------------------------------------------------


#==============================================================================
#READ INPUT FILE
#==============================================================================
os.system('clear')
print 'Hi-PEEC CLUSTER EXTRACTION SOFTWARE'
print 'Version 0.1'
print 'Last changed: {}'.format(time.ctime(os.path.getmtime('Hi-PEEC_pipeline.py')))
print ''

# Location of python code
pydir = os.getcwd()

# Define input file
infile = 'Hi-PEEC_settings.input'

# Verify that input file exists
if os.path.exists(pydir + '/' + infile) == False:
    print ''
    print 'File ', infile, ' could not be found in ', pydir
    sys.exit('quitting now...')

#Read in file
raw_userinput = np.genfromtxt(infile, dtype=None)

#convert the input matrix to a dictionary for easy access to the inputs
userinput = dict(zip(raw_userinput[:, 0], raw_userinput[:, 1]))

#Convert from strings to the true datatypes.
for key in userinput:
    try:
        userinput[key] = ast.literal_eval(str(userinput[key]))
    except:
        userinput[key] = str(userinput[key])

#Set up directories
if userinput['OUTDIR']==False:
    target_dir = os.getcwd()
    userinput['OUTDIR'] = target_dir
else:
    target_dir = userinput['OUTDIR']


# Print contents of input file
print 'Inputs that will be used by the pipeline:'
print '_______________________________________________________________________'
print ''
for key in userinput:
    print '{}:  {}'.format(key, userinput[key])
print ''
print '_______________________________________________________________________'





#==============================================================================
#RUNNING THE PIPELINE
#==============================================================================

#Setting up folder structure at desired location
if userinput['SETUP']:
    filemanagement.setup(userinput, pydir)

#Running initial extraction
if userinput['EXTRACT']:
    extraction.extraction(userinput)

# Running initial photometry on the isolated stars
if userinput['DO_PHOT']:
    extraction.photometry(userinput, userinput['IMAGE'],
                          userinput['STARS'], 'stars.mag', '3.0' )
