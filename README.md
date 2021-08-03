[![Build Status](https://travis-ci.org/andycasey/smhr.svg?branch=master)](https://travis-ci.org/andycasey/smhr)

Spectroscopy Made Harder
------------------------
Gotta pay back that tech debt.


Authors
-------
 - Andrew R. Casey (Monash)
 - Alex Ji (University of Chicago)


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
```
Note: if you previously used instructions where it said `conda config --set channel_priority strict` this makes installing on anaconda super slow; I would change this back to
`conda config --set channel_priority flexible`

* Install required libraries into the `smhr-py3` environment:
```
conda create --name smhr-py3 python=3.8 scipy numpy matplotlib=3.1.3 six astropy ipython python.app
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
pythonw __main__.py #pythonw is installed with python.app and fixes menubar issues
```
Note: if you use python or ipython on Big Sur, the menu bar may not work.
It appears you can fix this by clicking outside SMHR then clicking back in. But using pythonw is better.
Details: https://stackoverflow.com/questions/48738805/mac-pyqt5-menubar-not-active-until-unfocusing-refocusing-the-app
* Install moog17scat (see below) and add it to your path.

MOOG
----
It is currently recommended that you use this version of MOOG: https://github.com/alexji/moog17scat

This version is not the most efficient, and it computes synthetic spectra only to about 0.003 accuracy. It is modified from the 2017 February version of MOOG from Chris Sneden's website. It includes Jennifer Sobeck's scattering routines (turned on and off with the flag `scat`, which is not true in the default MOOG 2017) and the fixes to the Barklem damping that were implemented in the 2014 MOOG refactoring.

There is now a 2019 November version of MOOG, but it did not add anything different unless you use the HF molecule or work on combined spectra of globular clusters. It did also start forcing MOOG to read everything as loggf from linelists, rather than logging things if all the loggfs were positive. But in SMHR we add a fake line whenever this is detected, so it does impact anything here.

The 0.003 accuracy comes because this version of MOOG by default has a looser criterion for recomputing continuum opacity (compared to Jen's widely distributed version with scattering in 2011).
See the README for `moog17scat` if you have concerns.

Note that by default right now, we require you to have an executable called `MOOGSILENT` callable from your `$PATH` environment variable. Specifically, we use the version of MOOG that you get from `which MOOGSILENT`.

VERSION HISTORY:
----------------
- The current master branch is python 3.
- Alex has ported SMHR to python 3 in branch `py38-mpl313`. It now uses pyside2 and updated libraries for matplotlib. It is also way easier to install, not relying on some obscure libraries that were no longer maintained.
- The branch `refactor-scatterplot` has an updated and improved GUI (as of Jan 2020). These have not been merged into master yet but should be soon.
- Until Aug 2021, the `master` branch was frozen to a version from about July 2019.
- v0.22 (formerly branch `better-errors`) is a frozen version that is the result of a big update on May 30, 2019. It is considered a stable version.
- v0.2 is a frozen development version, v0.21 is a slightly more recently frozen version. 
- v0.1 is the current stable version. Things are working and it is being used for papers.

If you are new to SMHR, you should use the branch `refactor-scatterplot`.
Note v0.1 and v0.2 files are not compatible, but there is a script to convert old save files into new save files.
There is not a way to convert files from the old SMH to new SMHR.

TO DOs
------
-[ ] Fix GUI layout in Linux
