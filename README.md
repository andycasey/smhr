[![Build Status](https://travis-ci.org/andycasey/smhr.svg?branch=master)](https://travis-ci.org/andycasey/smhr)

Spectroscopy Made Harder
------------------------
Gotta pay back that tech debt.


Authors
-------
 - Andrew R. Casey (Monash)
 - Alex Ji (University of Chicago)
 - Erika Holmbeck (Carnegie Observatories)


Installation
------------

* Get anaconda

* Add conda-forge
```
conda config --add channels conda-forge
```

* Install required libraries into the `smhr-py3` environment:
```
conda create --name smhr-py3 python=3.8 scipy numpy matplotlib=3.1.3 six astropy ipython python.app requests
conda activate smhr-py3
conda install pyside2
conda install yaml
```

* Download and install this branch:
```
git clone https://github.com/andycasey/smhr.git 
cd smhr
python setup.py develop
```
Due to data not being copied, it is best to do `develop` instead of install. Someday we will fix this.

* Try running it:
```
cd smh/gui
pythonw __main__.py #pythonw is installed with python.app and fixes menubar issues
```

Note: you can also open smhr with `python` or `ipython`, but the menu bar may not work.
It appears you can fix this by clicking outside SMHR then clicking back in. But using `pythonw` is better.
Details: https://stackoverflow.com/questions/48738805/mac-pyqt5-menubar-not-active-until-unfocusing-refocusing-the-app

* Install moog17scat (see below) and add it to your path.

MOOG
----
It is currently recommended that you use this version of MOOG: https://github.com/alexji/moog17scat

Follow the usual MOOG installation instructions. When you compile MOOG, make sure that you have not activated any anaconda environments, because it can mess up the gfortran flags.
Note that SMHR requires you to have an executable called `MOOGSILENT` callable from your `$PATH` environment variable. Specifically, it uses the version of MOOG that you get from `which MOOGSILENT`.

This version is modified from the 2017 February version of MOOG from Chris Sneden's website. It includes Jennifer Sobeck's scattering routines (turned on and off with the flag `scat`, which is not true in the default MOOG 2017) and the fixes to the Barklem damping that were implemented in the 2014 MOOG refactoring.
There is now a 2019 November version of MOOG, but it did not add anything different unless you use the HF molecule or work on combined spectra of globular clusters. It did also start forcing MOOG to read everything as loggf from linelists, rather than logging things if all the loggfs were positive. But in SMHR we add a fake line whenever this is detected, so it does impact anything here.

Note that Alex has recently (Nov 16, 2021) fixed a bug in moog17scat that existed since the beginning and resulted in continuum accuracy only at the 0.003 when scattering is on. He also fixed a bug in isotopes.
See the README for `moog17scat` if you have concerns.


VERSION HISTORY:
----------------
- The current master branch is python 3.
- Alex has ported SMHR to python 3 in branch `py38-mpl313`. It now uses pyside2 and updated libraries for matplotlib. It is also way easier to install, not relying on some obscure libraries that were no longer maintained.
- The branch `refactor-scatterplot` has an updated and improved GUI (as of Jan 2020). These have not been merged into master yet but should be soon.
- Until Feb 2022, the `master` branch was frozen to a version from about July 2019.
- v0.22 (formerly branch `better-errors`) is a frozen version that is the result of a big update on May 30, 2019. It is considered a stable version.
- v0.2 is a frozen development version, v0.21 is a slightly more recently frozen version. 
- v0.1 is the current stable version. Things are working and it is being used for papers.
Note v0.1 and v0.2 files are not compatible, but there is a script to convert old save files into new save files.
There is not a way to convert files from the old SMH to new SMHR.

TO DOs
------
-[ ] Fix GUI layout in Linux
