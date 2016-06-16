#!/usr/bin/env python
#-*- coding: utf-8 -*-

#------------------------------------------------------------------------------
#Title: Hi-PEEC Python Pipeline
#Author: 		Axel Runnholm
#Creation date: 2016-06-16
#Description: 	This is the cluster extraction pipeline for the Hi-PEEC project
#				(PI: Angela Adamo). For information on setup and running of the
#				file see the manual and the README included in this folder.
#
#				The code is based on the LEGUS cluster extraction pipeline
#				(see Calzetti(2015) and Adamo(in prep) )
#------------------------------------------------------------------------------


#import packages
#------------------------------------------------------------------------------
# compensation for main diffs between Python 2.7 and 3
from __future__ import division,print_function

#import math, datahandling and plotting utils
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt

#sys utils
import os,glob
import time
import shutil
import sys
import string

#astronomy utils
import pyfits, pdb
from pyraf import iraf
import pywcs
#------------------------------------------------------------------------------


#==============================================================================
#READ INPUT FILE
#==============================================================================


# Location of target directory (should be current directory)
target_dir = os.getcwd()

# Location of /data directory and this code
pydir = target_dir 		#'/'.join(target_dir.split('/')[:-1])

# Open necessary packages in iraf
iraf.noao(_doprint=0)
iraf.digiphot(_doprint=0)
iraf.obsutil(_doprint=0)
iraf.daophot(_doprint=0)

# Define input file
infile = 'Hi-PEEC_settings.input'

# Verify that input file exists and load it in
if os.path.exists(target_dir + '/' + infile) == False:
    print ''
    print 'File ', infile, ' could not be found in ', target_dir
    sys.exit('quitting now...')

userinput = np.genfromtxt(infile,)

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
