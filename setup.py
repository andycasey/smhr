#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from setuptools import setup, find_packages
from codecs import open
from os import path, system
from re import compile as re_compile
from urllib import urlretrieve

# For convenience.
if sys.argv[-1] == "publish":
    system("python setup.py sdist upload")
    sys.exit()

def read(filename):
    kwds = {"encoding": "utf-8"} if sys.version_info[0] >= 3 else {}
    with open(filename, **kwds) as fp:
        contents = fp.read()
    return contents

# Get the version information.
here = path.abspath(path.dirname(__file__))
vre = re_compile("__version__ = \"(.*?)\"")
version = vre.findall(read(path.join(here, "smh", "__init__.py")))[0]

# External data.
if "--with-models" in map(str.lower, sys.argv):
    data_paths = [
        # Model photospheres:
        # Castelli & Kurucz (2004)
        ("https://zenodo.org/record/14964/files/castelli-kurucz-2004.pkl",
            "smh/photospheres/castelli-kurucz-2004.pkl"),
        # MARCS (2008)
        ("https://zenodo.org/record/14964/files/marcs-2011-standard.pkl",
            "smh/photospheres/marcs-2011-standard.pkl"),
        # Stagger-Grid <3D> (2013)
        ("https://zenodo.org/record/15077/files/stagger-2013-optical.pkl",
            "smh/photospheres/stagger-2013-optical.pkl"),
        ("https://zenodo.org/record/15077/files/stagger-2013-mass-density.pkl",
            "smh/photospheres/stagger-2013-mass-density.pkl"),
        ("https://zenodo.org/record/15077/files/stagger-2013-rosseland.pkl",
            "smh/photospheres/stagger-2013-rosseland.pkl"),
        ("https://zenodo.org/record/15077/files/stagger-2013-height.pkl",
            "smh/photospheres/stagger-2013-height.pkl"),
    ]
    for url, filename in data_paths:
        if path.exists(filename):
            print("Skipping {0} because file already exists".format(filename))
            continue
        print("Downloading {0} to {1}".format(url, filename))
        try:
            urlretrieve(url, filename)
        except IOError:
            raise("Error downloading file {} -- consider trying without the "
                "--with-models flag".format(url))
    sys.argv.remove("--with-models")


setup(
    name="smh",
    version=version,
    #author="",
    #author_email="",  # <-- Direct complaints to this address.
    description="Spectroscopy Made Harder",
    long_description=read(path.join(here, "README.md")),
    url="https://github.com/andycasey/smhr",
    license="MIT",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.4",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics"
    ],
    keywords="astronomy stellar spectroscopy spectra made easy hard",
    packages=find_packages(exclude=["documents", "tests"]),
    install_requires=[
        "numpy",
        "scipy>=0.14.0",
        "six",
        "pyside>=1.1.2",
        "astropy",
        "pyyaml"
        ],
    extras_require={
        "test": ["coverage"]
    },
    package_data={
        "": ["LICENSE"],
        "smh.gui": ["matplotlibrc"],
        "smh.photospheres": [
            "marcs-2011-standard.pkl",
            "castelli-kurucz-2004.pkl",
            "stagger-2013-optical.pkl",
            "stagger-2013-mass-density.pkl",
            "stagger-2013-rosseland.pkl",
            "stagger-2013-height.pkl"
        ],
        "smh.rt.moog": [
            "defaults.yaml",
            "abfind.in"
        ]
    },
    include_package_data=True,
    data_files=None,
    entry_points=None
)
