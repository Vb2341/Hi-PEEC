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

#----------------------------------------------------------------------------------------

#==============================================================================
#INPUT PARAMETERS
#==============================================================================


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

#Clear directory
#Clear old test and recreate directory

print('')
print('Clearing and creating directories')

if os.path.exists(test_dir) == False:
    os.makedirs(test_dir)
else:
    os.system('rm -r '+test_dir)
    os.makedirs(test_dir)


print('')
print('Copying files to the test directory')


#copy init dir.
source =  pipeline_dir+'init/'
destination = test_dir+'init/'
shutil.copytree(source, destination)

#copy python file
source =  pipeline_dir+'legus_clusters_extraction.py'
destination = test_dir+'legus_clusters_extraction.py'
shutil.copyfile(source, destination)

#copy imagefiles
for filename in glob.glob(os.path.join(source_dir, '*.*')):
    shutil.copy(filename, dest_dir)



img = input[3]
extraction_image =data_dir+img
source =  extraction_image
destination = test_dir+img
shutil.copyfile(source, destination)

