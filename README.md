# README #


### What is this repository for? ###
This repo contains the code for cluster extraction in Hi-PEEC based on the LEGUS cluster extraction script.

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
The second way is to simply download the zip file that is contained in the master branch. Note that this is updated less
often.

You download the file by going to ```Source``` in the left sidebar, clicking on the file in the filelist and then clicking ```Raw```


###How to run the pipeline###
1. Set up the required folder structure (see [Here](https://bitbucket.org/C0gito/hi-peec/wiki/projectstructure))

3. Check the Hi-PEEC_settings.input file to make sure that the inputs are appropriate.

2. Run ```python Hi-PEEC_pipeline.py``` in a terminal

###Dependencies###
Running the pipeline requires a working python 2.7 installation with the following packages installed:

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
```

####astronomy utils####
These are normally included in scientific distributions e.g. Ureka
```
pyfits
pyraf
pywcs
```
>>>>>>> devel
