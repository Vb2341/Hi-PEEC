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

def printProgress (iteration, total, prefix = '\t Progress', suffix = 'Complete', decimals = 2, barLength = 50):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    filledLength    = int(round(barLength * iteration / float(total)))
    percents        = round(100.00 * (iteration / float(total)), decimals)
    bar             = '#' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('%s [%s] %s%s %s\r' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        print("\n")

def remove_if_exists(filename):
    if os.path.exists(filename) == True:
        try:
            os.remove(filename)
            return 1
        except OSError:
            shutil.rmtree(filename)
            return 1
    else:
        return 0

def setup(userinputs,pydir):

    target_dir = userinputs['OUTDIR']

    print ''
    print 'Creating directories ...'
    print ''

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
    remove_if_exists(target_dir+'/init')

    print 'Copying files'
    #copy init directory to target directory
    source =  pydir+'/init'
    destination = target_dir+'/init'
    shutil.copytree(source, destination)



    # Only move images to /img directory if they aren't there already
    im_check = glob.glob(target_dir + '/img/*.fits')

    if len(im_check) == 0:
        #Move all the fits files
        print '\t Copying fits files ...'

        imlist = glob.glob(userinputs['DATA'] + '/*.fits')

        #print progressbar
        i=0
        l = len(imlist)
        printProgress(i, l, prefix = '\t Progress:', suffix = 'Complete', barLength = 50)

        for impath in imlist:
            #set path parameters
            source = impath
            image = impath.split('/')[-1]
            destination = target_dir + '/img/' + image

            #copy the file
            shutil.copy(source, destination)

            #Update the progressbar
            i = i+1
            printProgress(i, l, prefix = '\t Progress:', suffix = 'Complete', barLength = 50)

            # Update header of each image slightly
            pf = pyfits.open(destination, mode='update')
            pf[0].header['BUNIT']= 'ELECTRONS/S'
            pf.close()
        print('')

    #Copy the Source Extraction parameter files.
    source =  pydir + '/init/output.param'
    destination = target_dir+ '/s_extraction/output.param'
    shutil.copyfile(source, destination)

    source =  pydir + '/init/R2_wl_aa.config'
    destination = target_dir + '/s_extraction/R2_wl_aa.config'
    shutil.copyfile(source, destination)

    source =  pydir + '/init/default.nnw'
    destination = target_dir + '/s_extraction/default.nnw'
    shutil.copyfile(source, destination)

    #copy the file for star coordinates to photometry folder
    source = pydir + '/init/' + userinputs['STARS']
    if os.path.exists(source) == False:
        sys.exit("Coordinate file for isolated stars (keyword 'STARS')\
         was not found in the /init/ directory. Please add this file")
    else:
        destination = target_dir + '/photometry/' + userinputs['STARS']
        shutil.copyfile(source, destination)

def cleanup():
    pass

