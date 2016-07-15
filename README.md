# README #


### What is this repository for? ###
This repo contains the code for cluster extraction in Hi-PEEC based on the LEGUS cluster extraction script.

### Changelog ###

#### v 1.2 ####
- Add edge masking
This version adds the ability to mask the edges of the image in order to prevent false detections from the enhanced noise levels. 

If the option in the settings file is set to True the pipeline looks for a .reg file in the init directory. If it does not find one it will issue a prompt to the user asking if he/she wants to create one. If the answer is yes ds9 is launched with the correct image and linedrawing selected. This then allows the user to draw lines around the region within which they want to detect sources. The user then has to save these regions ( name does not matter as long as it is .reg). The save format has to be ds9 in image coords. If another format is given the pipeline skips this step.

### How do I get set up? ###

####Using Git####
The easiest way to get the required folder structure is using [__Git__](https://git-scm.com/). To do this run the following in a terminal:

```
git clone https://C0gito@bitbucket.org/C0gito/hi-peec.git
```

which will give a copy of the entire source tree in the last stable version, the 'master'  branch. If you want the 
development version of the code which is not guaranteed to work you can then run 

```
git fetch && git checkout devel
```
which will then load the devel branch making your directory mirror that. 

#####Staying up to date:#####
Doing the above means it is easy to stay up to date with the changes in the code. You do this simply by running 
```
git pull
```
in the directory you placed the source files in. Note that if you have changed any of the input files that are tracked 
by git this might give you an error saying something like ```You have uncommited changes```. To get around this simply 
run 
```
git stash
```
and then try to pull from the repository again.


####Without Git####
The second way is to simply download the zip file that is contained in the master branch. This archive is kept up-to-date with the latest version released to master.

You download the file by going to ```Source``` in the left sidebar, clicking on the file in the filelist and then clicking ```Raw```


###How to run the pipeline###
1. Set up the required folder structure (see [Here](https://bitbucket.org/C0gito/hi-peec/wiki/config&setup))

3. Check the Hi-PEEC_settings.input file to make sure that the inputs are appropriate.

2. Run ```python Hi-PEEC_pipeline.py``` in a terminal

###Dependencies###
The pipeline has been thus far only been tested on a Linux Kubuntu system, but should work on any unix based OS.

Running the pipeline requires a working installation of Source Extractor and  python 2.7 with the following packages installed:


####math, datahandling and plotting utils####
These are the large and heavy packages 
```
numpy
pandas
matplotlib.pyplot
```

####sys utils####
Most of these should be included in a standard python installation 
```
os, glob
time
shutil
sys
string
ast
datetime
warnings
```

####astronomy utils####
These are normally included in scientific distributions e.g. Ureka
```
pyfits
pyraf
pywcs
```