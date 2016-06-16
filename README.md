# README #


### What is this repository for? ###
This repo contains the code for cluster extraction in Hi-PEEC and primarily the scripts for running tests on the extraction.

### How do I get set up? ###

* In order to run the tests you need the following folder structure
```
#!python

/Data/
    /Object/
        -fitsfiles
/Pipeline/
    /init/
        default.nnw
        R2_wl_aa.config
        output.param
    /data/
        legus_galactic_extinction.tab
        legus_zeropoints.tab
    legus_clusters_extraction.py
    legus_clusters_extraction.input
    progress.py
    extraction_test.py
    
```
###How to run the test software###
1. Set up above folder structure

2. run 
```
python extraction.py
```
3. This will display a prompt where you are asked to select the object. This object name must match the filename ```/Object/``` above.

4. Next the program will display the contents of the legus_clusters_extraction.input file and ask if you want to change any settings. If the answer is yes the program will open vim allowing you to change the inputfile.

5. Next the program will check for the existence of a directory called ```Extraction_test```. If there is such a directory it will ask if you want to clear it before running the new test. If 'yes' or if there is no such directory the contents of the ```/Pipeline/``` dir and the fits files corresponding to the selected object will be copied over to the test dir.

6. The program now runs the extraction pipeline with the requested settings.