#!/usr/bin/env python
#-*- coding: utf-8 -*-

#------------------------------------------------------------------------------
#Title: Hi-PEEC Aperture correction
#Author:        Axel Runnholm
#Creation date: 2016-06-16
#Description:   This script contains the routines for estimating and applying
#				aperture corrections.

#------------------------------------------------------------------------------


#import packages
#------------------------------------------------------------------------------
# compensation for main diffs between Python 2.7 and 3
from __future__ import division#,print_function

#import math, data handling and plotting utils
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

#Import Hi-PEEC modules
sys.path.insert(0, './source/')
import filemanagement
#------------------------------------------------------------------------------
