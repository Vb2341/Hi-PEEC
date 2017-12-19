#!/usr/bin/env python
#-*- coding: utf-8 -*-

#------------------------------------------------------------------------------
#Title: Hi-PEEC Python Pipeline v1.4
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
import subprocess

from multiprocessing import Pool

import logging
#If there is a problem, set level=logging.DEBUG to get additional information in the log file
logging.basicConfig(filename='Hi-PEEC.log', filemode='w', level=logging.INFO,
                    format='%(levelname)s: %(filename)s:%(lineno)s:%(funcName)s() - %(asctime)s: %(message)s',
                    datefmt='%I:%M:%S')

#astronomy utils
from astropy.io import fits
from pyraf import iraf
from astropy import wcs
#Import Hi-PEEC modules
sys.path.insert(0, './source/')
import filemanagement
import extraction
import catmanagement
import linemask
#------------------------------------------------------------------------------


#==============================================================================
#READ INPUT FILE
#==============================================================================
os.system('clear')
print 'Hi-PEEC CLUSTER EXTRACTION SOFTWARE'
print 'Version 1.4'
print 'Last updated: {}'.format(time.ctime(os.path.getmtime('Hi-PEEC_pipeline.py')))
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
    sys.exit('Quitting now')

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
userinput['PYDIR'] = os.getcwd()

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
    filemanagement.setup(userinput)


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
    extraction_cat = target_dir + '/s_extraction/R2_wl_dpop_detarea.cat'


#------------------------------------------------------------------------------
#Removing edge detections
#------------------------------------------------------------------------------
if userinput['MASK_EDGES']:
    logging.info('Removing edge detections')
    print ''
    print 'Removing edge detections'
    print ''

    linemask.mask_edges(userinput, extraction_cat)

    # Create a new reg file
    # Save region file from extraction catalog

    xx, yy, fwhm, class_s, mag = np.loadtxt(extraction_cat, skiprows=5, unpack=True)
    extraction.create_regfile(userinput, xx, yy, 'edges_removed.reg')


#------------------------------------------------------------------------------
#Create growth curve:
#------------------------------------------------------------------------------
if userinput['DO_GROWTH']:
    logging.info('RUNNING GROWTH CURVE ANALYSIS')

    # Check if there is an additional single star file.
    try:
        add_stars = userinput['ADD_STARS']
        add_stars_list = add_stars.split(',')
        if len(add_stars_list) == 1:
            add_stars_list = [add_stars]

        add_filter = userinput['ADD_FILTER']
        add_filter_list = add_filter.split(',')
        if len(add_filter_list) == 1:
            add_filter_list = [add_filter]

        logging.info('Using additional starfiles {} for filters {}'
                     .format(add_stars_list,add_filter_list))

        if len(add_filter_list) != len(add_stars_list):
            if len(add_stars_list) != 1:
                logging.info('Mismatch in number of additional filters and files.')
                filemanagement.shutdown('Number of filters and additional starfiles does not match.')
        additional_starfiles = True
    except KeyError:
        additional_starfiles = False

    # Running initial photometry on the isolated stars & creating a growth curve
    print ''
    print 'Running initial photometry on the isolated stars'
    #create a string '1.0,2.0' etc up to 20
    growth_curve_apertures=','.join('{}'.format(float(i)) for i in range(1,21))

    #Do photometry for the growthcurve
    growth_catalog = extraction.photometry(userinput, userinput['IMAGE'],
                                userinput['STARS'], 'isolated_stars.mag',
                                growth_curve_apertures )
    #Create the growthcurve for ref filter
    print ''
    print 'Creating growth curve for {}'.format(userinput['REF_FILTER'])
    extraction.growth_curve(userinput, userinput['REF_FILTER'], growth_catalog)

    # Create a growth curve for any additional filters
    if additional_starfiles == True:
        for a in range(len(add_filter_list)):
            add_filter = add_filter_list[a]
            if len(add_stars_list) == 1:
                add_stars = add_stars_list[0]
            else:
                add_stars = add_stars_list[a]

            image_list = glob.glob(userinput['DATA'] + '/*' + add_filter + '*_dr?.fits')
            if not image_list:
                print 'No frame matching additional single star filter'
            else:
                image = image_list[0].split('/')[-1]

                # Do photometry for the growthcurve
                growth_catalog = extraction.photometry(userinput, image,
                                            add_stars, 'isolated_stars_{}.mag'.format(add_filter),
                                            growth_curve_apertures )
                #Create the growthcurve for ref filter
                print ''
                print 'Creating growth curve for {}'.format(add_filter)
                extraction.growth_curve(userinput, add_filter, growth_catalog)

#------------------------------------------------------------------------------
#Do science photometry:
#------------------------------------------------------------------------------
if userinput['DO_PHOT']:
    logging.info('RUNNING SCIENCE PHOTOMETRY')

    #Do an initial centering photometry run and compute concentration index
    print ''
    print 'Running science photometry step'
    print '\t Centering coordinates'
    extraction.photometry(userinput, userinput['IMAGE'],
                                   extraction_cat, 'centercoords.mag',
                                   '1,3',
                                   recenter=True)

    center_catalog = target_dir + '/photometry/centercoords.mag'

    fullcat = target_dir + '/photometry/fullcat_ref_center.coo'
    filemanagement.remove_if_exists(fullcat)

    iraf.txdump(center_catalog,'XCENTER,YCENTER,MAG[1-2]','yes',Stdout=fullcat)
    refx, refy, mag1, mag3 = np.genfromtxt(fullcat, missing_values='INDEF',
                                           filling_values=np.nan,unpack=True)
    concentration_index = mag1 - mag3
    header_ref = fits.getheader(userinput['DATA'] + '/' + userinput['IMAGE'],'SCI')
    wcs_ref = wcs.WCS(header_ref)
    ref_scale = fits.getval(userinput['DATA']+'/'+userinput['IMAGE'], 'D001SCAL', ext=0)
    # Calculate and Dec for xy coordinates. 1 refers to origin of image in ds9.
    ra, dec = wcs_ref.wcs_pix2world(refx, refy, 1)
    skycat = target_dir + '/photometry/fullcat_sky_center.coo'
    np.savetxt(skycat, np.array([ra,dec]).transpose())

    #Use this as input and do all image photometry

    print ''
    print '\tDoing photometry on full image set:'

    imagelist = glob.glob(userinput['DATA'] + '/*dr?.fits')

    #print progressbar
    i=0
    l = len(imagelist)
    # filemanagement.printProgress(i, l)

    # Wrapper for Parallelization
    def mproc_phot(image):
        filter = extraction.get_filter(image)

        outputfile = target_dir + '/photometry/phot_' + filter + '.mag'
        # calculate proper aperture to use
        ap = extraction.calc_aperture(userinput, image)

        #do photometry
        extraction.photometry(userinput, image, fullcat, outputfile, ap)

    # Parallelize it!
    # Pool(8).map(mproc_phot, imagelist)

    for image in imagelist:
        filter = extraction.get_filter(image)

        outputfile = target_dir + '/photometry/phot_' + filter + '.mag'
        # calculate proper aperture to use
        ap = extraction.calc_aperture(userinput, image)

        img_scale = fits.getval(image, 'D001SCAL', ext=0)
        ann = userinput['ANNULUS'] * ref_scale / img_scale
        dann = userinput['D_ANNULUS'] * ref_scale / img_scale

        # CHANGE COORDS HERE
        header_img = fits.getheader(image,'SCI')
        wcs_img = wcs.WCS(header_img)

        # Calculate RA and Dec for xy coordinates. 1 refers to origin of image in ds9.
        ix, iy = wcs_img.wcs_world2pix(ra, dec, 1)
        tmp_cat = target_dir + '/photometry/tmp.coo'
        np.savetxt(tmp_cat,np.array([ix,iy]).transpose())

        #do photometry
        extraction.photometry(userinput, image, tmp_cat, outputfile, ap, ann, dann)
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
if userinput['CREATE_CAT']:
    logging.info('CREATING FINAL CATALOG')
    print ''
    print 'Creating final photometric catalog'

    # Get a list of the photometric catalogs & sort them by wavelength
    phot_cats = glob.glob(target_dir + '/photometry/short_phot*')
    phot_cats = sorted(phot_cats, key=lambda file: (os.path.basename(file)))
    logging.debug('List of photometric catalogs to include in final:{}'.format(phot_cats))

    # Get a list of filters from the filenames
    filters = [os.path.basename(i).split('_')[-1].split('.')[0] for i in phot_cats]
    # Create the catalog dataframe
    cat = pd.DataFrame()
    print 'THE XY CATALOG IS {}'.format(phot_cats[filters.index(userinput['REF_FILTER'])])

    # Get XY coordinates for the reference filter.  NOT THE SAME FOR ALL FILTERS DUE TO DIFFERING PIXEL SCALES
    x, y = np.loadtxt(phot_cats[filters.index(userinput['REF_FILTER'])], unpack=True, usecols=(0,1))
    cat['X'] = x
    cat['Y'] = y
    cat['CI'] = concentration_index
    cat['CI'].replace(np.nan,99.999,inplace=True) # Replace CI nans with 99.999

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
    #------------------------------------------------------------------------------
    cat = catmanagement.insert_WCS(userinput,cat)



    # Save the final catalog to file
    #------------------------------------------------------------------------------
    print '\t Save the final catalog to file'
    date = datetime.datetime.now().strftime ("%Y-%m-%d")

    cat.reset_index(drop=True, inplace=True) #make sure ID numbers are reset to filtered cat.
    cat.to_csv(target_dir+'/final_cat_'+userinput['TARGET'] + '_' + date + '.cat',
               sep='\t', float_format = '%.8f')
    print ''
    nr_of_clusters = len (cat['X'])
    print '\t Number of clusters in the final catalogue: {}'.format(nr_of_clusters)
    logging.info('Number of clusters in the final catalogue: {}'.format(nr_of_clusters))

    # Create a region file from this catalog
    extraction.create_regfile(userinput, cat['X'], cat['Y'], 'final_catalog.reg')

#------------------------------------------------------------------------------
# MAKE CATALOGS FOR ADDED CLUSTERS
#------------------------------------------------------------------------------
if userinput['DO_CLUSTERS']:
    print ''
    print 'Running pipeline for manually added clusters'
    logging.info('RUNNING MANUALLY ADDED CLUSTER STEP')

    inputlist_path = target_dir + '/init/' + userinput['CLUSTERS']

    if os.path.exists(inputlist_path)==False:
        logging.critical('The specified cluster file {} does not exist.'\
                         .format(inputlist_path))
        logging.debug('The specified cluster file {} does not exist.'\
                         .format(inputlist_path))
        filemanagement.shutdown('The specified cluster file does not exist. Quitting.',userinput)

    # Select images
    ref_image = userinput['DATA'] + '/' + userinput['IMAGE']
    imagelist = glob.glob(userinput['DATA']+'/*dr?.fits')


    # Recenter the input coordinates
    #------------------------------------------------------------------------------
    print '\t Centering coordinates'
    extraction.photometry(userinput, userinput['IMAGE'],
                                   inputlist_path, 'manual_centers.mag',
                                   str(userinput['AP_RAD']))

    center_catalog = target_dir + '/photometry/manual_centers.mag'

    manualcat = target_dir + '/photometry/fullcat_ref_center.coo'
    filemanagement.remove_if_exists(manualcat)

    iraf.txdump(center_catalog,'XCENTER,YCENTER','yes',Stdout=manualcat)


    # Do photometry on the added cluster coordinates
    #------------------------------------------------------------------------------
    print '\t Doing photometry for manually added clusters on full image set:'


    #print progressbar
    i=0
    l = len(imagelist)
    filemanagement.printProgress(i, l)

    for image in imagelist:
        filter = extraction.get_filter(image)

        # Calculate appropriate aperture
        ap = extraction.calc_aperture(userinput,image)

        outputfile = target_dir + '/photometry/manual_'+filter+'.mag'
        extraction.photometry(userinput, image, manualcat, outputfile, ap)
        i=i+1
        filemanagement.printProgress(i, l)


    #------------------------------------------------------------------------------
    # Create the final photometric catalogs
    #------------------------------------------------------------------------------
    logging.info('CREATING MANUAL CATALOG')
    print ''
    print '\t Creating final photometric catalog for manually added clusters'

    # Get a list of the photometric catalogs & sort them by wavelength
    phot_cats = glob.glob(target_dir + '/photometry/short_manual*')
    phot_cats = sorted(phot_cats, key=lambda file: (os.path.basename(file)))
    logging.debug('List of photometric catalogs to include in manual cat:{}'.format(phot_cats))

    # Get a list of filters from the filenames
    filters = [os.path.basename(i).split('_')[-1].split('.')[0] for i in phot_cats]


    # Create the catalog dataframe
    mancat = pd.DataFrame()

    # Get the x y coordinates which are the same for all the filters.
    x, y = np.loadtxt(phot_cats[filters.index(userinput['REF_FILTER'])], unpack=True, usecols=(0,1))
    mancat['X'] = x
    mancat['Y'] = y

    # Generate column labels
    maglabels = ['Mag ' + s for s in filters]
    errlabels = ['Err ' + s for s in filters]

    # Construct the photometric catalog
    logging.info('Assembling the catalog')
    for a in range(len(phot_cats)):
        # Read in data, append to final catalog dataframe.
        mag, err = np.loadtxt(phot_cats[a], unpack=True, usecols=(2,3))
        mancat[maglabels[a]] = mag
        mancat[errlabels[a]] = err


    # Remove sources that are not detected in BVI
    #------------------------------------------------------------------------------
    print '\t Selecting sources detected in B,V and I'
    mancat = catmanagement.BVI_detection(mancat, filters,userinput)


    # Remove duplicate sources
    #------------------------------------------------------------------------------
    print '\t Removing duplicate sources'
    mancat = catmanagement.remove_duplicates(mancat)


    # Apply Extinction and aperture corrections
    #------------------------------------------------------------------------------
    print '\t Apply extinction and aperture corrections'
    mancat = catmanagement.apply_corrs(userinput,mancat)


    # Insert WCS coordinates into catalog
    #------------------------------------------------------------------------------
    mancat = catmanagement.insert_WCS(userinput,mancat)


    # Save the final catalog to file
    #------------------------------------------------------------------------------
    print '\t Save the final catalog to file'
    date = datetime.datetime.now().strftime ("%Y-%m-%d")

    mancat.reset_index(drop=True, inplace=True) #make sure ID numbers are reset to filtered cat.
    mancat.to_csv(target_dir+'/manual_cat_'+userinput['TARGET'] + '_' + date + '.cat',
               sep='\t', float_format = '%.3f')
    print ''
    nr_of_clusters = len (mancat['X'])
    print '\t Number of clusters in the final manual catalogue: {}'.format(nr_of_clusters)
    logging.info('Number of clusters in the final manual catalogue: {}'.format(nr_of_clusters))

#------------------------------------------------------------------------------
# FINAL CLEANUPS
#------------------------------------------------------------------------------


filemanagement.shutdown('All operations performed. Shutting down.',userinput)
