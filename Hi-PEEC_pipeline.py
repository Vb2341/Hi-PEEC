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
import pandas as pd
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
import apcorr
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

#------------------------------------------------------------------------------
#Setting up folder structure at desired location
#------------------------------------------------------------------------------
if userinput['SETUP']:
    filemanagement.setup(userinput, pydir)


#------------------------------------------------------------------------------
#Running initial extraction
#------------------------------------------------------------------------------
if userinput['EXTRACT']:
    print ''
    print 'Extracting sources from reference image:'
    print ''
    extraction_cat = extraction.extraction(userinput)


#------------------------------------------------------------------------------
#Create growth curve:
#------------------------------------------------------------------------------
if userinput['DO_GROWTH']:

    # Running initial photometry on the isolated stars & creating a growth curve
    print ''
    print 'Running initial photometry on the isolated stars'
    #create a string '1.0,2.0' etc up to 20
    growth_curve_apertures=','.join('{}'.format(float(i)) for i in range(1,21))

    #Do photometry for the growthcurve
    growth_catalog = extraction.photometry(userinput, userinput['IMAGE'],
                                userinput['STARS'], 'isolated_stars.mag',
                                growth_curve_apertures )
    #Create the growthcurve
    print ''
    print 'Creating growth curve'
    extraction.growth_curve(userinput, growth_catalog)


#------------------------------------------------------------------------------
#Do science photometry:
#------------------------------------------------------------------------------
if userinput['DO_PHOT']:

    #Do an initial centering photometry run
    print ''
    print 'Centering coordinates'
    extraction.photometry(userinput, userinput['IMAGE'],
                                   extraction_cat, 'centercoords.mag',
                                   str(userinput['AP_RAD']))

    center_catalog = target_dir + '/photometry/centercoords.mag'

    fullcat = target_dir + '/photometry/fullcat_ref_center.coo'
    filemanagement.remove_if_exists(fullcat)

    iraf.txdump(center_catalog,'XCENTER,YCENTER','yes',Stdout=fullcat)

    #Use this as input and do all image photometry


    print ''
    print 'Doing photometry on full image set'

    imagelist = glob.glob(target_dir+'/img/*sci*.fits')

    #print progressbar
    i=0
    l = len(imagelist)
    filemanagement.printProgress(i, l)

    for image in imagelist:
        try:
            filter = pyfits.getheader(image)['FILTER']
        except KeyError:
            #The 814 image has the filter information under the keyword FILTER2:
            filter = pyfits.getheader(image)['FILTER2']

        outputfile = target_dir + '/photometry/phot_'+filter+'.mag'
        extraction.photometry(userinput, image, fullcat, outputfile, str(userinput['AP_RAD']))
        i=i+1
        filemanagement.printProgress(i, l)


#------------------------------------------------------------------------------
# Calculate aperture corrections
#------------------------------------------------------------------------------

if userinput['APCORR']:
    print ''
    print 'Calculating aperture corrections'

    apcorr.calculation(userinput)


#------------------------------------------------------------------------------
# Create the final photometric catalogs
#------------------------------------------------------------------------------

# Get a list of the photometric catalogs & sort them by wavelength
phot_cats = glob.glob(target_dir + '/photometry/short_phot*')
phot_cats = sorted(phot_cats, key=lambda file: (os.path.basename(file)))

# Get a list of filters from the filenames
filters = [os.path.basename(i).split('_')[-1][0:5] for i in phot_cats]

# Get a list of images corresponding to the filters


final_cat =pd.DataFrame()
# Get the x y coordinates which are the same for all the filters.
x, y = np.loadtxt(phot_cats[filters==userinput['REF_FILTER']], unpack=True, usecols=(0,1))
final_cat['X'] = x
final_cat['Y'] = y

# Construct the photometric catalog
for a in range(len(phot_cats)):
    # Read in data, append to final catalog dataframe.
    mag, err = np.loadtxt(phot_cats[a], unpack=True, usecols=(2,3))
    mag_name = 'Mag ' + filters[a]
    err_name = 'Err ' + filters[a]
    final_cat[mag_name] = mag
    final_cat[err_name] = err

# Remove sources that are not detected in two contigous bands
