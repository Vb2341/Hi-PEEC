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

# Set science sky annulus and width here so that they can be changed
# without hunting for them
annulus_sci = 7.0
dannulus_sci = 1.0


# Set some IRAF phot parameters that will be common to all photometry runs
iraf.unlearn('datapars')
iraf.datapars.scale = 1.0
iraf.datapars.fwhmpsf = 2.0
iraf.datapars.sigma = 0.01
iraf.datapars.readnoise = 5.0
iraf.datapars.itime = 1.

iraf.unlearn('centerpars')
iraf.centerpars.cbox = 1
iraf.centerpars.cmaxiter = 3
iraf.centerpars.maxshift = 1

iraf.unlearn('fitskypars')
iraf.fitskypars.salgori = 'mode'

iraf.unlearn('photpars')

iraf.unlearn('phot')
iraf.phot.interactive = 'no'
iraf.phot.verbose = 'yes'
iraf.phot.verify = 'no'
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


def photometry(userinputs,catalog,outputname):

    #Load zeropoints
    inst_zp, filter_zp, zp_zp = np.loadtxt(pydir + '/data/legus_zeropoints.tab', unpack=True, dtype='str')

    #Get list of imagefiles
    imlist = glob.glob(target_dir + '/img/*_sci.fits')

    # Finds the full path of reference image from the filter
    ref_image = userinputs['OUTDIR'] + '/img/' + userinputs['IMAGE']


    # Set the necessary variables for photometry on the reference image
    ref_exptime = pyfits.getheader(ref_image)['EXPTIME']
    ref_inst = pyfits.getheader(ref_image)['INSTRUME']
    ref_inst = ref_inst.lower()

    match = (inst_zp == ref_inst) & (filter_zp == ref_filter)
    ref_zp = zp_zp[match]

    # ref_zp is a string within an array, so need to turn into a float
    ref_zp = float(ref_zp[0])
