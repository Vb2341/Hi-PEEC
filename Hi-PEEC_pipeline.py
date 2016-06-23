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
print ''
print 'Creating final photometric catalog'

# Get a list of the photometric catalogs & sort them by wavelength
phot_cats = glob.glob(target_dir + '/photometry/short_phot*')
phot_cats = sorted(phot_cats, key=lambda file: (os.path.basename(file)))

# Get a list of filters from the filenames
filters = [os.path.basename(i).split('_')[-1][0:5] for i in phot_cats]

# Get a list of images corresponding to the filters


cat = pd.DataFrame()

# Get the x y coordinates which are the same for all the filters.
x, y = np.loadtxt(phot_cats[filters==userinput['REF_FILTER']], unpack=True, usecols=(0,1))
cat['X'] = x
cat['Y'] = y

# Generate column labels
maglabels = ['Mag ' + s for s in filters]
errlabels = ['Err ' + s for s in filters]

# Construct the photometric catalog
for a in range(len(phot_cats)):
    # Read in data, append to final catalog dataframe.
    mag, err = np.loadtxt(phot_cats[a], unpack=True, usecols=(2,3))
    cat[maglabels[a]] = mag
    cat[errlabels[a]] = err


# Remove sources that are not detected in two contigous bands
#------------------------------------------------------------------------------
print '\t Removing sources not detected in 2 contiguous filters'
maxerr = userinput['MAGERR']
sel = ((cat[errlabels[0]] <= maxerr) & (cat[errlabels[1]] <= maxerr)) | \
         ((cat[errlabels[1]] <= maxerr) & (cat[errlabels[2]] <= maxerr)) | \
         ((cat[errlabels[2]] <= maxerr) & (cat[errlabels[3]] <= maxerr)) | \
         ((cat[errlabels[3]] <= maxerr) & (cat[errlabels[4]] <= maxerr))

# Remove all sources that do not match the above criterion
cat = cat.drop(cat[~sel].index)



# Assign number of filters
#------------------------------------------------------------------------------
print '\t Assigning number of filters'
nfilt_4 = ((cat[errlabels[0]] <= maxerr) & (cat[errlabels[1]] <= maxerr)
          & (cat[errlabels[2]] <= maxerr) &(cat[errlabels[3]] <= maxerr)) \
          | \
          ((cat[errlabels[1]] <= maxerr) & (cat[errlabels[2]] <= maxerr) &
          (cat[errlabels[3]] <= maxerr) & (cat[errlabels[4]] <= maxerr))

nfilt_5 = ((cat[errlabels[0]] <= maxerr) & (cat[errlabels[1]] <= maxerr)
          & (cat[errlabels[2]] <= maxerr) & (cat[errlabels[3]] <= maxerr)
          & (cat[errlabels[4]] <= maxerr))

cat['# filters'] = 2
cat['# filters'][cat[nfilt_4].index] = 4
cat['# filters'][cat[nfilt_5].index] = 5



# Remove duplicate sources
#------------------------------------------------------------------------------
print '\t Removing duplicate sources'

# Set tolerance to use
tolerance = 0.5

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



# Apply Extinction and aperture corrections
#------------------------------------------------------------------------------
print '\t Apply extinction and aperture corrections'
# Assign the galactic extinction for the five instrument/filter combinations
# Read in extinction file as array (a vs u means acs/uvis)
extfile = target_dir + '/init/Hi-PEEC_galactic_extinction.tab'
extinctions = pd.read_table(extfile, header=[0,1],index_col=0, sep=r"\s*")

# Pick extinctions for our target only
target = userinput['TARGET'].lower()

extinctions = extinctions.loc[['ngc1614']]

imlist = glob.glob(target_dir + '/img/*_sci.fits')


# Load apcorrs
apf, apc, ape = np.genfromtxt(target_dir + '/photometry/avg_aperture_correction.txt', unpack=True, usecols=(0, 1, 2))

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

    # Apply extinctions and apcorrs
    for a in range(len(mag)):
        if mag[a]!=99.999 and mag[a]!=66.666:
            mag[a] = mag[a] + ap_corr - gal_ext
            err[a] = np.sqrt(err[a]**2 + ap_err**2)

    # Insert the new mags into the catalog
    cat['Mag' + filter] = mag
    cat['Err ' + filter] = err


# Make an array of cluster ids
ids = np.arange(len(cat['X'])) + 1

# Convert xy coordinates into RA Dec of reference filter
ref_image = [image for image in imlist if userinput['REF_FILTER'] in image][0]

# Get header from reference image using pyfits
header_ref = pyfits.getheader(ref_image)

# Get wcs solution from reference image header
wcs_ref = pywcs.WCS(header_ref)

# Calculate RA and Dec for xy coordinates. 1 refers to origin of image in ds9.
ra, dec = wcs_ref.wcs_pix2sky(cat['X'], cat['Y'], 1)
cat.insert(2,'RA',ra)
cat.insert(3,'DEC',dec)

# Calculate apparent magnitude of absolute mag limit (-6)
distance = userinput['DISTANCE']

dist_mod = 5 * np.log10(distance) + 25.
m_apparent = dist_mod - 6.


# Save the final catalog to file
#------------------------------------------------------------------------------
print '\t Save the final catalog to file'
date = datetime.datetime.now().strftime ("%Y-%m-%d")
#roundcat = cat.map(lambda x: '%2.1f' % x)
cat = cat.reset_index()
cat.to_csv(target_dir+'/final_cat_'+userinput['TARGET'] + date + '.cat', sep='\t', float_format = '%.3f')
