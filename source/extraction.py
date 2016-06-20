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

iraf.noao(_doprint=0)
iraf.digiphot(_doprint=0)
iraf.obsutil(_doprint=0)
iraf.daophot(_doprint=0)


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
iraf.phot.verbose = 'no'
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


def photometry(userinputs, image, catalog, outputname, apertures):
    """Function for performing IRAF photometry
    @Params:
    userinputs  (dict)  - dictionary with results from the user input file.
    image       (STR)   - fitsfile we want to do photometry on
    catalog     (STR)   - input coordinates where we want to do photometry
    outputname  (STR)   - name of the file where we store the measured results
    apertures   (STR)   - what apertures to measure. Should be a string i.e. '1.0,3.0'
    """
    #set directory
    target_dir = userinputs['OUTDIR']

    #Update passed names  to be full paths
    image = glob.glob(target_dir + '/img/' + image)
    if len(image)==0:
        sys.exit('Selected image does not exist')
    else:
        image = image[0]

    catalog = target_dir + '/init/' + catalog
    output = target_dir + '/photometry/' + outputname


    #Load zeropoints
    inst_zp, filter_zp, zp_zp = np.loadtxt(target_dir + '/init/legus_zeropoints.tab', unpack=True, dtype='str')

    #get filtername from imagename
    filter = image.split('/')[-1].split('_')[0]


    # Set the necessary variables for photometry on the reference image
    exptime = pyfits.getheader(image)['EXPTIME']
    inst = pyfits.getheader(image)['INSTRUME']
    inst = inst.lower()

    match = (inst_zp == inst) & (filter_zp == filter.lower())
    zp = zp_zp[match]

    # zp is a string within an array, so need to turn into a float
    zp = float(zp[0])

    # Remove output file if it already exists
    if os.path.exists(output) == True:
        os.remove(output)

    # Run photometry
    iraf.datapars.epadu = exptime

    iraf.centerpars.calgorithm = 'centroid'

    iraf.fitskypars.annulus = userinputs['ANNULUS']
    iraf.fitskypars.dannulu = userinputs['D_ANNULUS']

    iraf.photpars.apertures = apertures
    iraf.photpars.zmag = zp

    iraf.phot(image, catalog, output)

    # Move phot results to new files, remove INDEFs
    #--------------------------------------------------------------------------
    fullcat_mag_short = target_dir + '/photometry/short_' + outputname

    # Removes all outputlines that do not contain the character '*'
    # ensures only phot results are kept
    cmd = 'grep "*" ' + output + ' > ' + fullcat_mag_short
    os.system(cmd)

    # Replace INDEFS:
    cmd = 'sed -i.bak "s/INDEF/99.999/g" ' + fullcat_mag_short
    os.system(cmd)

    # Remove .bak files to prevent confusion
    bak_fullcat = fullcat_mag_short + '.bak'
    os.remove(bak_fullcat)
    #--------------------------------------------------------------------------

    return fullcat_mag_short


def growth_curve(userinputs, catalog):
    aper_st, flux_st = np.loadtxt(catalog, unpack=True, usecols=(0,3))

    ratio_st = np.empty(len(aper_st))

    #number of apertures
    naper = 20

    # Calculate the number of stars, make sure it is an integer
    nstar = int(len(aper_st)/naper)
    aper_ind = naper - 1

    for k in range(nstar):

        for i in range(naper):

            ratio_st[i + k*naper] = flux_st[i + k*naper]/flux_st[aper_ind + k*naper]



    # Find median ratio at each aperture between all the stars and all the clusters
    med_st = np.empty(naper)

    for i in range(naper):

        med_st[i] = np.median(ratio_st[i::naper])


    # Plot growth curves
    fig = plt.figure(figsize = (7,7))

    aper_x = np.arange(naper) + 1

    for i in range(nstar):

        ratio_y = ratio_st[i*naper:(i + 1)*naper]
        plt.plot(aper_x, ratio_y, 'y-')
        plt.annotate(str(i + 1), xy=(8.0, ratio_y[7]),
            horizontalalignment='left', verticalalignment='top', fontsize=6)


    plt.plot(aper_x, med_st, 'r-' , linewidth=4.0)
    plt.hlines(0.5, 0, 20, color='black', linewidth=2, zorder=10)
    plt.vlines(4, 0, 1.1, color='black', linewidth=2, linestyle='dashed', zorder=10)
    plt.vlines(5, 0, 1.1, color='black', linewidth=2, linestyle='dashed', zorder=10)
    plt.vlines(6, 0, 1.1, color='black', linewidth=2, linestyle='dashed', zorder=10)

    plt.ylabel('Normalized Flux ' + ref_filter.upper())
    plt.xlabel('Radius (pix)')
    plt.xlim(1,20)
    plt.minorticks_on()

    fig.savefig(userinputs['OUTDIR'] + '/plots/plot_growth_curve.pdf')

