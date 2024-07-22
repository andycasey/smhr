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

* Get anaconda https://www.anaconda.com/

* Create a new environment and install required libraries:
For M1 macs, note that pyside2 has to be run in Rosetta. Thus, you can install it in this way:
```
conda create -c conda-forge/osx-64 --name smhr-py3 python=3.8 scipy numpy=1.22.4 matplotlib=3.1.3 six astropy ipython python.app requests pyside2=5.13.2 yaml
```
Currently (as of May 2022) anaconda on M1/ARM chips by default includes channels that search through `osx-arm64` and `noarch` but not `osx-64`.
Also, newer versions of pyside2 appear to have changed some syntax on dialog boxes. We will update this eventually but for now you can install the older pyside2 version.

For older Macs or other computers, this worked fine:
```
conda create -c conda-forge --name smhr-py3 python=3.8 scipy numpy=1.22.4 matplotlib=3.1.3 six astropy ipython python.app requests
conda activate smhr-py3
conda install -c conda-forge pyside2=5.13.2
conda install -c conda-forge yaml
```

* Download and install this branch:
```
git clone https://github.com/andycasey/smhr.git 
cd smhr
git checkout -b py38-mpl313
git pull origin py38-mpl313
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

* Some installation notes for Linux/Debian. It takes a very long time to install pyside2 (hours?) so be patient. Thanks to Shivani Shah and Terese Hansen for this information.
```
Install Python 3.7 from anaconda

create a new environment for smhr-py3:
> conda create --name smhr-py3 python=3.8 scipy numpy=1.22.4 matplotlib=3.1.3 six astropy ipython requests

Activate environment:
> conda activate smhr-py3

Install pyside2:
> conda install -c conda-forge pyside2=5.13.2


Install yaml
> conda install -c conda-forge yaml

Get smhr:
> git clone https://github.com/andycasey/smhr.git 
> cd smhr
> python setup.py develop

Start smhr:
> cd smh/gui
> ipython __main__.py
```

MOOG
----
It is currently recommended that you use this version of MOOG: https://github.com/alexji/moog17scat

Follow the usual MOOG installation instructions. When you compile MOOG, make sure that you have not activated any anaconda environments, because it can mess up the gfortran flags.
Note that SMHR requires you to have an executable called `MOOGSILENT` callable from your `$PATH` environment variable. Specifically, it uses the version of MOOG that you get from `which MOOGSILENT`.

This version is modified from the 2017 February version of MOOG from Chris Sneden's website. It includes Jennifer Sobeck's scattering routines (turned on and off with the flag `scat`, which is not true in the default MOOG 2017) and the fixes to the Barklem damping that were implemented in the 2014 MOOG refactoring.
There is now a 2019 November version of MOOG, but it did not add anything different unless you use the HF molecule or work on combined spectra of globular clusters. It did also start forcing MOOG to read everything as loggf from linelists, rather than logging things if all the loggfs were positive. But in SMHR we add a fake line whenever this is detected, so it does impact anything here.

Note that Alex has recently (Nov 16, 2021) fixed a bug in moog17scat that existed since the beginning and resulted in continuum accuracy only at the 0.003 when scattering is on. He also fixed a bug in isotopes.
See the README for `moog17scat` if you have concerns.
(Note May 2022: Alex has updated the master branch of moog17scat so this is done by default.)


VERSION HISTORY:
----------------
- The current master branch is python 3.
- March 5, 2024: starting from 4b7732ceaff1ba1bff9e5c36b891b2c0a8ab03a3, Alex has updated the mean abundances. Two important changes: (1) on the stellar parameters tab the average [Fe I,II/H] was previously reported as the median abundance. It now shows the mean abundance in parentheses. (2) on the review tab and in previous abundance summaries, I had mistakenly applied a weight based on the statistical uncertainty instead of using the straight mean abundance. It now uses no weights to calculate the mean.
- Alex has ported SMHR to python 3 in branch `py38-mpl313`. It now uses pyside2 and updated libraries for matplotlib. It is also way easier to install, not relying on some obscure libraries that were no longer maintained. Also added damping selection.
- The branch `refactor-scatterplot` has an updated and improved GUI (as of Jan 2020). These have not been merged into master yet but should be soon.
- Until Feb 2022, the `master` branch was frozen to a version from about July 2019.
- v0.22 (formerly branch `better-errors`) is a frozen version that is the result of a big update on May 30, 2019. It is considered a stable version.
- v0.2 is a frozen development version, v0.21 is a slightly more recently frozen version. 
- v0.1 is the current stable version. Things are working and it is being used for papers.

Note v0.1 and v0.2 files are not compatible, but there is a script to convert old save files into new save files.
There is not a way to convert files from the old SMH to new SMHR.
