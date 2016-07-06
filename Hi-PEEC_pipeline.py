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
import matplotlib.pyplot as plt

#sys utils
import os, glob
import time
import shutil
import sys
import string
import ast
import datetime

import logging
#If there is a problem, set level=logging.DEBUG to get additional information in the log file
logging.basicConfig(filename='Hi-PEEC.log', filemode='w', level=logging.DEBUG,
                    format='%(levelname)s: %(filename)s:%(lineno)s:%(funcName)s() - %(asctime)s: %(message)s',
                    datefmt='%I:%M:%S')

#astronomy utils
import pyfits
from pyraf import iraf
import pywcs

#Import Hi-PEEC modules
sys.path.insert(0, './source/')
import filemanagement
import extraction
import catmanagement
#------------------------------------------------------------------------------


#==============================================================================
#READ INPUT FILE
#==============================================================================
os.system('clear')
print 'Hi-PEEC CLUSTER EXTRACTION SOFTWARE'
print 'Version 0.1'
print 'Last changed: {}'.format(time.ctime(os.path.getmtime('Hi-PEEC_pipeline.py')))
print ''
logging.info('Run started')


# Location of python code
pydir = os.getcwd()

# Define input file
infile = 'Hi-PEEC_settings.input'

# Verify that input file exists
if os.path.exists(pydir + '/' + infile) == False:
    print ''
    print 'File ', infile, ' could not be found in ', pydir
    logging.critical('Could not find input file, quit process')
    logging.debug('Looking for {} '.format(pydir + '/' + infile))
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
    logging.info('No output directory specified, using {}'.format(target_dir))
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
    logging.info('RUNNING SETUP')
    filemanagement.setup(userinput, pydir)


#------------------------------------------------------------------------------
#Running initial extraction
#------------------------------------------------------------------------------
if userinput['EXTRACT']:
    logging.info('RUNNING EXTRACTION')
    print ''
    print 'Extracting sources from reference image:'
    print ''
    extraction_cat = extraction.extraction(userinput)
else:
    target_dir + '/s_extraction/R2_wl_dpop_detarea.cat'


#------------------------------------------------------------------------------
#Create growth curve:
#------------------------------------------------------------------------------
if userinput['DO_GROWTH']:
    logging.info('RUNNING GROWTH CURVE ANALYSIS')

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
    logging.info('RUNNING SCIENCE PHOTOMETRY')

    #Do an initial centering photometry run
    print ''
    print 'Running science photometry step'
    print '\t Centering coordinates'
    extraction.photometry(userinput, userinput['IMAGE'],
                                   extraction_cat, 'centercoords.mag',
                                   str(userinput['AP_RAD']))

    center_catalog = target_dir + '/photometry/centercoords.mag'

    fullcat = target_dir + '/photometry/fullcat_ref_center.coo'
    filemanagement.remove_if_exists(fullcat)

    iraf.txdump(center_catalog,'XCENTER,YCENTER','yes',Stdout=fullcat)

    #Use this as input and do all image photometry


    print ''
    print '\tDoing photometry on full image set:'

    imagelist = glob.glob(target_dir+'/img/*sci*.fits')

    #print progressbar
    i=0
    l = len(imagelist)
    filemanagement.printProgress(i, l)

    for image in imagelist:
        filter = extraction.get_filter(image)

        outputfile = target_dir + '/photometry/phot_'+filter+'.mag'
        extraction.photometry(userinput, image, fullcat, outputfile, str(userinput['AP_RAD']))
        i=i+1
        filemanagement.printProgress(i, l)


#------------------------------------------------------------------------------
# Calculate aperture corrections
#------------------------------------------------------------------------------

if userinput['APCORR']:
    logging.info('RUNNING APCORR STEP')
    print ''
    print 'Calculating aperture corrections'

    catmanagement.apcorr_calc(userinput)


#------------------------------------------------------------------------------
# Create the final photometric catalogs
#------------------------------------------------------------------------------
logging.info('CREATING FINAL CATALOG')
print ''
print 'Creating final photometric catalog'

# Get a list of the photometric catalogs & sort them by wavelength
phot_cats = glob.glob(target_dir + '/photometry/short_phot*')
phot_cats = sorted(phot_cats, key=lambda file: (os.path.basename(file)))
logging.debug('List of photometric catalogs to include in final:{}'.format(phot_cats))

# Get a list of filters from the filenames
filters = [os.path.basename(i).split('_')[-1][0:5] for i in phot_cats]

# Create the catalog dataframe
cat = pd.DataFrame()

# Get the x y coordinates which are the same for all the filters.
x, y = np.loadtxt(phot_cats[filters==userinput['REF_FILTER']], unpack=True, usecols=(0,1))
cat['X'] = x
cat['Y'] = y

# Generate column labels
maglabels = ['Mag ' + s for s in filters]
errlabels = ['Err ' + s for s in filters]

# Construct the photometric catalog
logging.info('Assembling the catalog')
for a in range(len(phot_cats)):
    # Read in data, append to final catalog dataframe.
    mag, err = np.loadtxt(phot_cats[a], unpack=True, usecols=(2,3))
    cat[maglabels[a]] = mag
    cat[errlabels[a]] = err


# Remove sources that are not detected in BVI
#------------------------------------------------------------------------------
print '\t Selecting sources detected in B,V and I'
cat = catmanagement.BVI_detection(cat, filters,userinput)


# Remove duplicate sources
#------------------------------------------------------------------------------
print '\t Removing duplicate sources'
cat = catmanagement.remove_duplicates(cat)



# Apply Extinction and aperture corrections
#------------------------------------------------------------------------------
print '\t Apply extinction and aperture corrections'
cat = catmanagement.apply_corrs(userinput,cat)


# Insert WCS coordinates into catalog
cat = catmanagement.insert_WCS(cat)



# Save the final catalog to file
#------------------------------------------------------------------------------
print '\t Save the final catalog to file'
date = datetime.datetime.now().strftime ("%Y-%m-%d")

cat.reset_index(drop=True, inplace=True) #make sure ID numbers are reset to filtered cat.
cat.to_csv(target_dir+'/final_cat_'+userinput['TARGET'] + '_' + date + '.cat',
           sep='\t', float_format = '%.3f')
print ''
nr_of_clusters = len (cat['X'])
print '\t Number of clusters in the final catalogue: {}'.format(nr_of_clusters)
logging.info('Number of clusters in the final catalogue: {}'.format(nr_of_clusters))


#------------------------------------------------------------------------------
# FINAL CLEANUPS
#------------------------------------------------------------------------------

# Move the log file to the outdirectory
os.chdir(pydir)
cmd = 'mv Hi-PEEC.log ' + target_dir + '/Hi-PEEC.log'
os.system(cmd)
