from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import matplotlib.pyplot as plt
import time

import smh
from smh import Session
from smh.spectral_models import (ProfileFittingModel, SpectralSynthesisModel)
from smh.linelists import LineList

if __name__=="__main__":
    path = "/Users/alexji/lib/smhr/smh/tests/test_data/test2.smh"
    session = Session.load(path)
    # Already normalized and rv-corrected, but not stitched
    normalized_spectrum = session.stitch_and_stack()
    
    metadata = session.metadata
    line_list = metadata["line_list"]
    spectral_models = metadata["spectral_models"]
    
    """
    # Everything is fit already, but you can run fit again
    print("Fitting models..."); start = time.time()
    for model in spectral_models:
        try:
            res = model.fit()
        except Exception as e:
            print(model.elements, model.species, model._repr_wavelength, e)
    print("Time to fit all models: {:.1f}".format(time.time()-start))
    """

    # Abundances measured already, but let's do it here again to illustrate
    # Also some don't have uncertainties yet (BUG) so we have to do this
    print("Measuring abundances..."); start = time.time()
    session.measure_abundances()
    print("Time to measure abundances: {:.1f}".format(time.time()-start))

    ## Plot full spectrum
    fig,ax = plt.subplots()
    ax.plot(normalized_spectrum.dispersion, normalized_spectrum.flux,'k-',lw=.1)
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"Normalized flux")
    ax.set_ylim((0,1.2))

    ## Plot Ba4554 line
    fig,ax = plt.subplots()
    ii = np.logical_and(4570 >= normalized_spectrum.dispersion, 
                        normalized_spectrum.dispersion >= 4540)
    ax.plot(normalized_spectrum.dispersion[ii], normalized_spectrum.flux[ii],'k-',lw=.1)
    for model in spectral_models:
        if model.species[0] == 56.1 and float(model._repr_wavelength) < 4560 \
                and float(model._repr_wavelength) > 4550:
            break
    assert model.species[0] == 56.1
    fitted_parameters = model.metadata["fitted_result"][2]
    ax.plot(fitted_parameters["model_x"],fitted_parameters["model_y"], 'r-')

    ## Plot abundance versus wavelength, reduced equivalent width, excitation potential
    ## Gather the information from the spectral model list
    ## Only use acceptable models (i.e. checkmarked ones)
    is_acceptable = []
    species = []
    wavelengths = []
    excitation_potentials = []
    equivalent_widths = []
    abundances = []
    abundance_uncertainties = []
    for model in spectral_models:
        if isinstance(model, ProfileFittingModel):
            if not model.is_acceptable: continue

            transition = model.transitions[0]
            assert model.species[0] == transition["species"]
            species.append(model.species[0])

            wavelengths.append(transition["wavelength"])
            excitation_potentials.append(transition["expot"])

            try:
                fitted_result = model.metadata["fitted_result"]
            except KeyError:
                equivalent_widths.append(np.nan)
                abundances.append(np.nan)
                abundance_uncertainties.append(np.nan)
            else:
                equivalent_widths.append(fitted_result[2]["equivalent_width"][0]*1000.)
                # In some cases, a tiny EW causes MOOG to fail.
                try:
                    # model.abundances[0] gives the abundance, but tries to automatically
                    # run MOOG if it doesn't see an abundance. So it crashes.
                    abundances.append(fitted_result[2]["abundances"][0])
                    abundance_uncertainties.append(fitted_result[2]["abundance_uncertainties"][0])
                except (smh.radiative_transfer.RTError, IndexError, KeyError):
                    abundances.append(np.nan)
                    abundance_uncertainties.append(np.nan)
        else: #Synthesis; don't do for now
            pass

    species = np.array(species)
    wavelengths = np.array(wavelengths)
    excitation_potentials = np.array(excitation_potentials)
    equivalent_widths = np.array(equivalent_widths)
    abundances = np.array(abundances)
    abundance_uncertainties = np.array(abundance_uncertainties)

    fe_ii = np.array(species.astype(int) == 26)

    fig,axes = plt.subplots(3,1)
    ax = axes[0]
    ax.scatter(excitation_potentials[fe_ii], abundances[fe_ii], c='k')
    ax.set_xlabel(r"$\chi$")
    ax.set_ylabel(r"$A(Fe)$")
    
    ax = axes[1]
    ax.scatter(np.log10(equivalent_widths[fe_ii]/wavelengths[fe_ii]), abundances[fe_ii], c='k')
    ax.set_xlabel(r"$\log({\rm EW}/\lambda)$")
    ax.set_ylabel(r"$A(Fe)$")
    
    ax = axes[2]
    ax.scatter(wavelengths[fe_ii], abundances[fe_ii], c='k')
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$A(Fe)$")
    # fig.savefig("scatterplots.png", bbox_inches="tight")
    
    plt.show()
