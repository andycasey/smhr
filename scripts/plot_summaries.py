from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import time
import yaml

import smh
from smh import Session
from smh.spectral_models import (ProfileFittingModel, SpectralSynthesisModel)
from smh.linelists import LineList

def fit_continuum(session, index):
    kwds = {}
    keys = ("function", "order", "low_sigma_clip", "high_sigma_clip",
            "knot_spacing", "max_iterations")
    for key in keys:
        kwds[key] = session.setting(("normalization", key))
    kwds["full_output"]=True
    # TODO should apply masks here
    order = session.input_spectra[index]
    try:
        normalized_spectrum, continuum, left, right \
            = order.fit_continuum(**kwds)
    except:
        print("Failed to fit {}".format(index))
    else:
        session.metadata["normalization"]["continuum"][index] = continuum
        session.metadata["normalization"]["normalization_kwargs"][index] = kwds

if __name__=="__main__":
    #prefix="/Users/alexji/MIKE_data/2016aug/70slit"
    #prefix="/Users/alexji/MIKE_data/2016aug/70slit"
    #prefix="/Users/alexji/MIKE_data/2016aug/35slit"
    #prefixblue = prefix+"/blue/Final-Products"
    #prefixred  = prefix+"/red/Final-Products"
    prefix="/Users/alexji/Dropbox/Observing/2016_8/reduced_data/rproc"
    prefixblue=prefix
    prefixred=prefix
    #outpath = "./n3"
    outpath = "/Users/alexji/Dropbox/Observing/2016_8/reduced_data/rproc/summaries"
    rv_fix = {} #{"retii-09":60.0}
    
    specnames = [x[:-15] for x in glob.glob(prefixblue+"/*blue*")]
    specnames = [os.path.basename(x) for x in specnames]
    print(specnames)

    for specname in specnames:
        bluespec = prefixblue+"/"+specname+"blue_multi.fits"
        redspec  = prefixred +"/"+specname+"red_multi.fits"
        session = Session([bluespec,redspec])
        #with open(smh.Session._default_settings_path, "rb") as fp:
        #    defaults = yaml.load(fp)
        #session.metadata.update(defaults)

        if specname in rv_fix:
            rv = rv_fix[specname]
        else:
            rv,rve = session.rv_measure()
        print("RV:",rv)
        session.rv_correct(rv)
        #session.normalize_input_spectra()
        for index,order in enumerate(session.input_spectra):
            fit_continuum(session,index)
        
        session.stitch_and_stack()
        #summary_kwds = session.setting(("summary_figure"))
        fig = session.make_summary_plot()
        fig.savefig(outpath+"/"+specname+"_summary.png")
        fig2 = session.make_ncap_summary_plot()
        fig2.savefig(outpath+"/"+specname+"_ncap_summary.png")
        session.save(outpath+"/"+specname+".smh",overwrite=True)
        plt.close(fig)
