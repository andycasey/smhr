[![Build Status](https://travis-ci.org/andycasey/smhr.svg?branch=master)](https://travis-ci.org/andycasey/smhr)

Spectroscopy Made Harder
------------------------
Gotta pay back that tech debt.


Authors
-------
 - Andrew R. Casey (Monash)
 - Alex Ji (Carnegie Observatories)

VERSION NOTE:
- The master branch is currently up to date and being actively developed.
- v0.22 (formerly branch `better-errors`) is a frozen version that is the result of a big update on May 30, 2019. It is considered a stable version.
- v0.2 is a frozen development version, v0.21 is a slightly more recently frozen version. 
- v0.1 is the current stable version. Things are working and it is being used for papers.

If you are new to SMHR, you should use the master branch or v0.22.
Note v0.1 and v0.2 files are not compatible, but there is a script to convert old save files into new save files.
There is not a way to convert files from the old SMH to new SMHR.

Installation
------------
This is one way that I (Alex) got things running from a fresh mac install and I'm putting it here as a record of some things I had to do.
Default anaconda does not officially support pyside, but I think conda-forge does. We may switch to that, or a different version of PySide, in the future. For now, this works...

- Download full anaconda with python 2.7. If you run into trouble later, specifically install an old anaconda version from [the anaconda archive](https://repo.continuum.io/archive/) around October 2016.
- `conda install pyside`: this will install pyside v1.1.1
- install qt4:
```
brew tap cartr/qt4
brew tap-pin cartr/qt4
brew install qt@4
```
- `conda install matplotlib=1.5.1` (fixing this is very painful because I need to update to qt5; I have started doing this but it will take a long time)
- `conda install numpy=1.11.3` (this is an issue with pyside)
- `conda install scipy=0.19.0` (this is an issue with newer versions of scipy, not sure yet why but it segfaults in `interpolate.griddata`)
- `conda install qt=4.8.7` (you may have to uninstall and downgrade some things for this to work; it should be safe to upgrade those later)
- Clone smhr
- Go into the smhr directory and `python setup.py develop`
- Go to `smhr/smh/gui` and open with `ipython __main__.py`. It should crash with a message about `libpyside`, saying something is not found. This is because a file is named wrong within anaconda.
  - Go to `~/anaconda/lib` (or the equivalent if you made a separate environment), and make a symlink or file copy for `libpyside-python2.7.1.1.dylib` from one that looks almost the same (e.g. `libpyside-python2.7.1.1.1.dylib`)
  - Try to open SMHR again, and it will crash. Do the same thing for `libshiboken` when you see this error message.
- If you have problems with `qt_menu.nib`, use `~/anaconda/bin/python.app __main__.py` instead of `ipython`. (I have fixed this on my laptop by copying it somewhere but I cannot find where.) Some possible places that could help:
```
~/anaconda/python.app/Contents/Resources/qt_menu.nib
~/anaconda/pkgs/launcher-1.0.0-1/launcherapp/Contents/Resources/qt_menu.nib
~/anaconda/pkgs/launcher-1.0.0-2/launcherapp/Contents/Resources/qt_menu.nib
~/anaconda/pkgs/python.app-1.2-py27_3/pythonapp/Contents/Resources/qt_menu.nib
~/anaconda/pkgs/python.app-1.2-py27_4/pythonapp/Contents/Resources/qt_menu.nib
~/anaconda/Launcher.app/Contents/Resources/qt_menu.nib
/usr/local/Cellar/qt/4.8.7_1/lib/QtGui.framework/Versions/4/Resources/qt_menu.nib
```
- There are sometimes some problems with segfaults. We hope this will go away when we move away from pyside. Sorry.
