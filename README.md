[![Build Status](https://travis-ci.org/andycasey/smhr.svg?branch=master)](https://travis-ci.org/andycasey/smhr)

Spectroscopy Made Harder
------------------------
Gotta pay back that tech debt.


Authors
-------
 - Andrew R. Casey (Monash)
 - Alex Ji (Carnegie Observatories)

VERSION NOTE:
- The branch `refactor-scatterplot` has an updated and improved GUI (as of Jan 2020). These have not been merged into master yet but should be soon.
- The `master` branch is currently frozen to a version from about July 2019.
- v0.22 (formerly branch `better-errors`) is a frozen version that is the result of a big update on May 30, 2019. It is considered a stable version.
- v0.2 is a frozen development version, v0.21 is a slightly more recently frozen version. 
- v0.1 is the current stable version. Things are working and it is being used for papers.

If you are new to SMHR, you should use the branch `refactor-scatterplot`.
Note v0.1 and v0.2 files are not compatible, but there is a script to convert old save files into new save files.
There is not a way to convert files from the old SMH to new SMHR.


Note about this version
------------------------
 - This is a fork of the original SMHr code that has updates and additions for use by the *R*-Process Alliance.
 - Direct questions to Erika Holmbeck (RIT/Notre Dame) or Alex Ji (Carnegie Observatories).


Installation
------------

* Get anaconda
* Add conda-forge
```
conda config --add channels conda-forge
conda config --set channel_priority strict
```
* Install required libraries:
```
conda create --name smhr-py3 python=3.8 scipy numpy matplotlib=3.1.3 six astropy ipython
conda activate smhr-py3
conda install pyside2
conda install yaml
```
* Download and install this branch:
```
git clone https://github.com/andycasey/smhr.git 
cd smhr
git checkout -b py38-mpl313
git pull origin py38-mpl313
python setup.py develop
```
* Try running it:
```
cd smh/gui
ipython __main__.py #ipython is often needed for the menu bar it seems
```
* Install moog17scat (see below) and add it to your path.



MOOG
----
It is currently recommended that you use this version of MOOG: https://github.com/alexji/moog17scat

This version is not the most efficient, and it computes synthetic spectra only to about 0.003 accuracy. It is modified from the 2017 February version of MOOG from Chris Sneden's website (there has been an update made in Nov 2019 that Alex has not investigated yet). It includes Jennifer Sobeck's scattering routines (turned on and off with the flag `scat`, which is not true in the default MOOG 2017) and the fixes to the Barklem damping that were implemented in the 2014 MOOG refactoring.

Note that by default right now, we require you to have an executable called `MOOGSILENT` callable from your `$PATH` environment variable. Specifically, we use the version of MOOG that you get from `which MOOGSILENT`. We may change this in the future.
