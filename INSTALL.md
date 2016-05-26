
Installation
============

0.  Download and install [git](https://git-scm.com/downloads).

1.  Follow the 'Setting up Git' and 'Authenticating with GitHub from Git' guides [here](https://help.github.com/articles/set-up-git/).

2.  Install [Anaconda Python](https://www.continuum.io/downloads). For now, select Python 2.7 until some [PySide issues](https://github.com/andycasey/smhr/issues/14) are resolved.

3.  Clone SMH from the [GitHub repository](https://github.com/andycasey/smhr) using this terminal command:

    ``git clone git@github.com:andycasey/smhr.git``
    
    This will create a directory called `smhr`, where everything lives.
  
3.  Finally, you will need `MOOGSILENT` installed. If you don't have it already installed, first make sure you have a FORTRAN compiler installed (e.g., [gfortran](https://gcc.gnu.org/wiki/GFortran)). In OSX you can install `gfortran` and other command line tools with the terminal command `xcode-select --install`. Then to install `MOOGSILENT`:

    ``pip install moogsilent --user``
  
4.  Install `PySide` using the following terminal command:

    ``conda install pyside``

5.  Move to the new `smhr` directory and install the code:

    ````
    cd smhr
    python setup.py develop
    ````

That's it. Currently to run SMH, type the following from the `smhr` folder:

````python
ipython
run -i smh/gui/__main__.py
````

(Eventually this will be packaged into a `smh` command line tool or packaged into a clickable application.)
