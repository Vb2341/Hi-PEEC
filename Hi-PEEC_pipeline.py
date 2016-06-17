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
from __future__ import division#,print_function

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
import ast

#astronomy utils
import pyfits, pdb
from pyraf import iraf
import pywcs
#------------------------------------------------------------------------------


#==============================================================================
#READ INPUT FILE
#==============================================================================
os.system('clear')
print 'Hi-PEEC CLUSTER EXTRACTION SOFTWARE'
print 'Version 0.1'
print 'Last changed: {}'.format(time.ctime(os.path.getmtime('Hi-PEEC_pipeline.py')))
print ''
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

# Verify that input file exists
if os.path.exists(target_dir + '/' + infile) == False:
    print ''
    print 'File ', infile, ' could not be found in ', target_dir
    sys.exit('quitting now...')

#Read in file
raw_userinput = np.genfromtxt(infile, dtype=None)

#convert the input matrix to a dictionary for easy access to the inputs
userinput = dict(zip(raw_userinput[:,0], raw_userinput[:,1]))

#Convert from strings to the true datatypes.
for key in userinput:
	try:
		userinput[key] = ast.literal_eval(str(userinput[key]))
	except:
		userinput[key] = str(userinput[key])


# Print contents of input file
print 'Submitted inputs:'
print '________________________________________________________________________________________'
print ''
for key in userinput:
	print '{}:  {}'.format(key,userinput[key])
print ''
print '________________________________________________________________________________________'
