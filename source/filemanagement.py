#!/usr/bin/env python
#-*- coding: utf-8 -*-

#------------------------------------------------------------------------------
#Title: Setup
#Author:        Axel Runnholm
#Creation date: 2016-06-17
#Description:   This script contains the functionality for setting up the
#               required directories and files needed for running as well as
#               the function responsible for cleaning the directory at the end
#               of the run.
#------------------------------------------------------------------------------


#import packages
#------------------------------------------------------------------------------
# compensation for main diffs between Python 2.7 and 3
from __future__ import division

#import math, datahandling and plotting utils
import numpy as np

#sys utils
import os,glob
import time
import shutil
import sys
import string
import ast

#astronomy utils
import pyfits

#------------------------------------------------------------------------------


def setup(userinputs,pydir):
    #setup directory paths
    if userinputs['OUTDIR']==False:
        target_dir = os.getcwd()
    else:
        target_dir = userinputs['OUTDIR']



    print 'Creating directories ...'
    #Create required directories
    if os.path.exists(target_dir + '/s_extraction') == False:
        os.makedirs(target_dir + '/s_extraction')

    if os.path.exists(target_dir + '/photometry') == False:
        os.makedirs(target_dir + '/photometry')

    if os.path.exists(target_dir + '/plots') == False:
        os.makedirs(target_dir + '/plots')

    if os.path.exists(target_dir + '/img') == False:
        os.makedirs(target_dir + '/img')


    #check that the init directory exists
    if os.path.exists(pydir + '/init') == False:
        print 'No init/ directory. Please create.'
        sys.exit()

    #Remove it if it exists to make sure that the config files are properly updated
    if os.path.exists(target_dir + '/init') == True:
        os.system('rm -r '+target_dir+'/init')

    #copy init directory to target directory
    source =  pydir+'/init'
    destination = target_dir+'/init'
    shutil.copytree(source, destination)


    # Only move images to /img directory if they aren't there already
    im_check = glob.glob(target_dir + '/img/*.fits')

    if len(im_check) == 0:
        #Move all the fits files
        print 'Moving fits files ...'

        imlist = glob.glob(userinputs['DATA'] + '/*.fits')

        for impath in imlist:
            source = impath
            image = impath.split('/')[-1]
            destination = target_dir + '/img/' + image
            os.rename(source, destination)

            # Update header of each image slightly
            pf = pyfits.open(destination, mode='update')
            pf[0].header['BUNIT']= 'ELECTRONS/S'
            pf.close()


    source =  pydir + '/init/output.param'
    destination = target_dir+ '/s_extraction/output.param'
    shutil.copyfile(source, destination)

    source =  pydir + '/init/R2_wl_aa.config'
    destination = target_dir + '/s_extraction/R2_wl_aa.config'
    shutil.copyfile(source, destination)

    source =  pydir + '/init/default.nnw'
    destination = target_dir + '/s_extraction/default.nnw'
    shutil.copyfile(source, destination)


def cleanup():
    pass
