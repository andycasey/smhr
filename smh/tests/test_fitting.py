"""
Based on Andy's pull request example
Currently just fits a bunch of individual lines with both synthesis and EW
TODO expand to do synthesis with multiple lines

TODO calculate abundances
TODO export tables
"""

import os

import smh.linelists
import smh.photospheres
import smh.radiative_transfer as rt
import smh.spectral_models as sm

import matplotlib.pyplot as plt

def show_result(region):
    fig, ax = plt.subplots()

    a, b, c = region._result
    ax.plot(c["model_x"], c["model_y"], c='r')
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    ax.plot(
        region.session.normalized_spectrum.dispersion,
        region.session.normalized_spectrum.flux,
        c='k', zorder=-1)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    l, u = c["model_yerr"]
    ax.fill_between(c["model_x"], c["model_y"] + l, c["model_y"] + u,
        facecolor="r", alpha=0.25, zorder=-2)

    return fig

if __name__=="__main__":
    datadir = os.path.dirname(os.path.abspath(__file__))+'/test_data'
    
    ## Test with "known" stellar parameters
    session = smh.Session([
            datadir+"/spectra/hd122563.fits"
            ])
    session.metadata['stellar_parameters']['effective_temperature'] = 4380
    session.metadata['stellar_parameters']['surface_gravity'] = 0.1
    session.metadata['stellar_parameters']['metallicity'] = -2.68
    session.metadata['stellar_parameters']['microturbulence'] = 2.65
    
    isotopes = {"H-C":{112:0.8,113:0.2},
                "C-N":{1214:0.8,1314:0.2},
                'Ba':{134:0.0,135:0.370,136:0.0,137:0.351,138:0.279}
                }
    session.metadata['isotopes'] = isotopes
    
    ## Start with a normalized, RV corrected spectrum
    spectrum = smh.specutils.Spectrum1D.read(datadir+"/spectra/hd122563.fits")
    session.normalized_spectrum = spectrum
    
    ## Primary line list (EW) measurements
    transitions = smh.linelists.LineList.read(datadir+"/linelists/complete.list")
    ew_measurements = []
    synth_measurements = []
    for transition in transitions:
        elem = transition['elem1']
        ew = sm.ProfileFittingModel(transition, session)
        synth = sm.SpectralSynthesisModel(transition, session, elem)
        ew_measurements.append(ew)
        synth_measurements.append(synth)
    
    ## Summarize
    ## EW Table
    
    ## Abundance Table
    
    ## Default plots
    
