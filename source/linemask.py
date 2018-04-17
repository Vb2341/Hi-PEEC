#!/usr/bin/env python
#-*- coding: utf-8 -*-

#------------------------------------------------------------------------------
#Title: Linemask
#Authors:       Matthew Hayes & Axel Runnholm
#Creation date: 2016-07-13
#Description:   This script contains the routines for masking the image edges to
#               prevent spurious detections.

#------------------------------------------------------------------------------


#import packages
#------------------------------------------------------------------------------
# compensation for main diffs between Python 2.7 and 3
from __future__ import division#,print_function

import sys
import traceback
import scipy as S
import numpy as np
from astropy.io import fits
from astropy import wcs
import glob
import os
import subprocess
import shutil
import logging
from shutil import copyfile

sys.path.insert(0, './source/')
import filemanagement

#-------------------------------------------------------------------------------



def getlines(fn):
    """Function for getting lines from ds9 reg file
    Input:
        fn (str) - filepath to regfile
    Outputs:
        d (iterable) - iterable with lines from regfile
    """
    fh = open(fn, "r")
    d  = fh.readlines()
    fh.close()
    return d

def convert_to_decimal(RA,DEC):

    RAhms = RA.split(':')
    DEChms = DEC.split(':')
    h = float(DEChms[0])

    # RA
    ra = 15.*(float(RAhms[0]) +float(RAhms[1])/60. + float(RAhms[2])/3600.)

    # DEC
    if h < 0:
        dec = h - float(DEChms[1])/60. - float(DEChms[2])/3600.
    else:
        dec = h + float(DEChms[1])/60. + float(DEChms[2])/3600.

    return ra,dec

def get_ends(l,header):
    try:
        x1 = float(l[0])
        y1 = float(l[1])
        x2 = float(l[2])
        y2 = float(l[3])
    except ValueError:
        wcs_ref = wcs.WCS(header)

        ra1,dec1 = convert_to_decimal(l[0],l[1])
        ra2,dec2 = convert_to_decimal(l[2],l[3])

        x1,y1 = wcs_ref.all_world2pix(ra1, dec1, 1)
        x2,y2 = wcs_ref.all_world2pix(ra2, dec2, 1)

    if x1 == x2: x1 += 1.  # do not permit x1=x2, otherwise gradient undefined.
    return x1, y1, x2, y2

def get_nearedge(x1,y1,x2,y2, Nx, Ny):
    cx = (x1+x2)/2.
    cy = (y1+y2)/2.
    edge_arr = S.array([ cx, Nx-cx, cy, Ny-cy ])    # distance to [ left, right, bottom,  top ]
    if   edge_arr.argmin() == 0: return "left"
    elif edge_arr.argmin() == 1: return "right"
    elif edge_arr.argmin() == 2: return "bottom"
    elif edge_arr.argmin() == 3: return "top"
    else :
        print "something screwy finding the nearest edge!"
        sys.exit(1)

def get_fn(x1,y1,x2,y2):
    m = (y2-y1)/(x2-x1)
    c = y1 - m*x1              # c = y - mx
    return m,c

def create_mask(ref, lines):
    """Function for creating a boolean mask
    Inputs:
        ref (str) - path to reference image
        lines (str) - path to region file
    Outputs:
        mask (np.array) - matrix with 1 if we are inside the selected region and 0 if we are outside.
    """
    # Main function
    fnRef   = ref
    fnLines = lines

    dat  = [ l.replace("line(", "").split(")")[0].split(",")  for l in getlines(fnLines) if l.startswith("line") ]

    head   = fits.getheader(fnRef)
    head1   = fits.getheader(fnRef, 1)
    img    = fits.getdata(fnRef)
    mask   = S.ones_like(img)
    Ny, Nx = img.shape

    dx  = S.zeros_like(mask)
    dy  = S.zeros_like(mask)
    for i in range(Nx): dx[:,i] = i+1
    for i in range(Ny): dy[i]   = i+1

    for d in dat :
        x1, y1, x2, y2 = get_ends(d,head1)
        #print x1, y1, "  ->  ", x2, y2
        edge_near = get_nearedge(x1, y1, x2, y2, Nx, Ny)
        m, c = get_fn(x1, y1, x2, y2)
        #print "grad=", m, "    intercept=",c

        ineq = None
        if   edge_near == "top"    : ineq = "gt"
        elif edge_near == "bottom" : ineq = "lt"
        elif edge_near == "left":
            if   m>0: ineq = ineq = "gt"
            else    : eneq = ineq = "lt"
        elif edge_near == "right":
            if   m>0: ineq = ineq = "lt"
            else    : eneq = ineq = "gt"

        if ineq == "gt": ind = dy > (m*dx + c)
        else           : ind = dy < (m*dx + c)
        mask[ind]  = 0.

    fits.writeto('mask.fits', mask, header=head1, clobber=True)

    return mask

def remove_edgedetections(catalog, ref, lines):
    """Function for removing clusters from the catalog if they are not inside the region specified.
    Inputs:
        ref (str) - path to reference image
        lines (str) - path to region file
    Outputs:
        none
    Effects:
        removes detections outside defined region from cluster catalog
    """

    mask = create_mask(ref, lines)

    xx, yy, fwhm, class_s, mag = np.loadtxt(catalog, skiprows=5, unpack=True)


    # Read in and keep sources from catalog
    with open(catalog, 'r') as f:
        sourcelines = f.readlines()

    # Non destructively assemble new file
    newLines = []
    i = 0
    for line in sourcelines:
        if line.startswith('#'):
            newLines.append(line)

        else:
            # Use round and int to make sure that indices are proper integers
            # errors using this should be <= half pixel
            x = int(np.round(xx[i]))
            y = int(np.round(yy[i]))
            if mask[y, x] == 1:
                newLines.append('{}\t{}\t{}\t{}\t{}\n'.format(xx[i],
                                                              yy[i],
                                                              fwhm[i],
                                                              class_s[i],
                                                              mag[i]))
            # Increment linecount
            i += 1

    # If nothing bad occured in setting up the new lines: write back to file
    with open(catalog, 'w') as f:
        for l in newLines:
            f.write(l)


def select_regfile(files):
    if len(files) > 1:
        print 'There are multiple reg files in the init directory'
        logging.info('Multiple .reg files found')
        for a in range(len(files)):
            print '{}: {}'.format(a,files[a])

        while True:
            try:
                choice = int(raw_input('Choose file ({}-{}): '.format(0,len(files)-1)))
                try:
                    file = files[choice]
                    return file
                except IndexError:
                    print 'Not a valid choice. Please choose one of the indicated numbers'
            except ValueError:
                print 'Not a valid choice. Please choose one of the indicated numbers'
    else:
        file = files[0]
        return file

def mask_edges(userinput,extraction_cat):
    ref_image =userinput['DATA'] + '/' + userinput['IMAGE']
    regfilename = userinput['OUTDIR'] + '/init/*.reg'

    if not glob.glob(regfilename):
        print 'There is no .regfile in the init directory'
        inp = raw_input('Do you wish to create one? (y/n) ')
        if inp == 'y':
            os.chdir(userinput['PYDIR'] + '/init')
            cmd = 'ds9 ' + ref_image + ' -scale log -scale limits -0.1 100 -regions shape line'
            process  = subprocess.Popen(cmd, shell=True)
            process.wait()
            print 'Checking for saved image'
            os.chdir(userinput['OUTDIR'])
            if not glob.glob(userinput['PYDIR'] + '/init/*.reg'):
                print 'Still no .reg file detected. Skipping the edge-removal step.'
                logging.info('No regfile was detected. Edge-removal was skipped.')
                return 0
            else:

                files = glob.glob(userinput['PYDIR'] + '/init/*.reg')

                file = select_regfile(files)
                file = file.split('/')[-1]

                shutil.copyfile(userinput['PYDIR'] + '/init/' + file,userinput['OUTDIR'] + '/init/' + file)
                regfile = glob.glob(regfilename)[0]
        elif inp =='n':
            print 'Skipping edge removal step'
            return 0
        else:
            print 'Not a recognised input. Skipping edge-removal'
            return 0
    else:
        files = glob.glob(regfilename)
        regfile = select_regfile(files)

    try:
        print 'Removing sources outside mask.'
        print 'Using {} to mask sources'.format(regfile)
        remove_edgedetections(extraction_cat, ref_image, regfile)
        return 1
    except Exception as e:
        print 'Edge masking failed. Proceeding without it.'
        logging.warning('Edgemasking failed. No mask has been applied')
        logging.warning('Edgemask traceback: {}'.format(traceback.print_exc()))
        print traceback.print_exc()
        return 0
