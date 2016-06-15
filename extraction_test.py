#-*- coding: utf-8 -*-

#----------------------------------------------------------------------------------------
'''
#Title: extraction_test.py
#Author: Axel Runnholm
#Creation date: 2016-06-14
#Description: Script for running tests on Legus extraction pipeline.
'''
#----------------------------------------------------------------------------------------

#import packages
#----------------------------------------------------------------------------------------
from __future__ import division,print_function
from progress import printProgress

import numpy as np
import matplotlib.pyplot as plt

import shutil
import os
import glob
from time import gmtime, strftime

#----------------------------------------------------------------------------------------

#==============================================================================
#INPUT PARAMETERS
#==============================================================================
os.system('clear')
print('RUNNING TESTSCRIPT FOR LEGUS_CLUSTERS_EXTRACTION')
currtime = strftime("%H:%M:%S")
print('Starting time:{}'.format(currtime))

target = 'NGC1614'
print('')
print('Object currently selected: {}'.format(target))

edit = raw_input('Do you wish to change this? (y/n)(n)')
if (edit == 'y'):
    target = raw_input('New target:')

#Set directories
current_dir = os.getcwd()
pipeline_dir = '/home/user/Hi-PEEC/Pipeline/'
data_dir = '/home/user/Hi-PEEC/Data/{}/'.format(target)
test_dir = '/home/user/Hi-PEEC/Extraction_test/'


#Load and show input file parameters.
inputfile = 'legus_clusters_extraction.input'
input = np.loadtxt(pipeline_dir+inputfile, unpack=True, skiprows=0, dtype='str')

print('')
print('Settings read from inputfile:')
print('------------------------------------------------------------------------')

with open(pipeline_dir + inputfile) as f:
    print(f.read())


print('------------------------------------------------------------------------')
print('')


edit = raw_input('Do you wish to edit these settings? (y/n)(n)')
if (edit == 'y'):
    os.system('vim ' + pipeline_dir + inputfile)

#==============================================================================
#MAINSECTION
#==============================================================================

#Clear directory?
def resetdir():
    """Function for clearing old tests and recreating directory
    Removes old test directory and then recreates it, copying over all necessary files
    from Pipeline and Data directories
    """

    print('')
    print('Clearing and creating directories')

    if os.path.exists(test_dir) == False:
        os.makedirs(test_dir)
    else:
        os.system('rm -r '+test_dir)
        os.makedirs(test_dir)


    print('')
    print('Copying extraction files to the test directory')


    #copy init dir.
    source =  pipeline_dir+'init/'
    destination = test_dir+'init/'
    shutil.copytree(source, destination)

    #copy data dir.
    source =  pipeline_dir+'data/'
    destination = test_dir+'data/'
    shutil.copytree(source, destination)

    #copy python file
    source =  pipeline_dir+'legus_clusters_extraction.py'
    destination = test_dir+'legus_clusters_extraction.py'
    shutil.copyfile(source, destination)

    #copy inputfile
    source =  pipeline_dir+'legus_clusters_extraction.input'
    destination = test_dir+'legus_clusters_extraction.input'
    shutil.copyfile(source, destination)

    #copy IRAF login.cl
    source =  pipeline_dir+'login.cl'
    destination = test_dir+'login.cl'
    shutil.copyfile(source, destination)


    #copy imagefiles
    print('')
    print('Copying imagefiles to test directory')
    files = glob.glob(os.path.join(data_dir, '*.*'))
    i = 0
    l =len(files)

    printProgress(i, l, prefix = ' Progress:', suffix = 'Complete', barLength = 50)
    for filename in files:
        shutil.copy(filename, test_dir)
        i+=1
        printProgress(i, l, prefix = ' Progress:', suffix = 'Complete', barLength = 50)

if os.path.exists(test_dir) == False:
    resetdir()
else:
    print('')
    clear = raw_input('Do you wish to clear the test directory before running? (y/n)(n)')
    if (clear == 'y'):
        resetdir()



#==============================================================================
#RUN THE EXTRACTION PIPELINE
#==============================================================================
print('------------------------------------------------------------------------')
print('')
print('Running the LEGUS extraction pipeline')
print('')

os.chdir(test_dir)
os.system('python legus_clusters_extraction.py')




print('')
print('------------------------------------------------------------------------')
print('')
endtime = strftime("%H:%M:%S")
print('Ending time:{}'.format(endtime))
