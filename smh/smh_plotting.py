from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["make_summary_plot"]

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator, ScalarFormatter
import os
import logging

from . import specutils
from . import spectral_models

logger = logging.getLogger(__name__)

def make_summary_plot(summary_figure, normalized_spectrum,
                      figure=None):
    """
    :param summary_figure:
        A dict containing information about how to plot the summary figure
    :param normalized_spectrum:
        specutils.Spectrum1D
    :param figure:
        matplotlib figure object to make plot in.
        The function will adjust the figure object
        If not specified (None), creates a new figure and returns it.
    """
    all_axes_names = ["top_left", "top_right", "middle", "bottom_left", "bottom_right"]
    # Validate file names
    spectra_filenames = []
    Tlit = []
    Teff = []
    for spectrum_list in summary_figure["spectra_filenames"]:
        spectra_filenames.append(spectrum_list[0])
        Tlit.append(int(spectrum_list[1]))
        Teff.append(int(spectrum_list[2]))
    for key in all_axes_names:
        for spectrum_to_plot in summary_figure[key]["spectra_to_plot"]:
            assert spectrum_to_plot in spectra_filenames, spectrum_to_plot

    # Color palette generated from seaborn.color_palette("hls",6)
    hls_colors = [(0.86, 0.37119999999999997, 0.33999999999999997),
                  (0.82879999999999987, 0.86, 0.33999999999999997),
                  (0.33999999999999997, 0.86, 0.37119999999999997),
                  (0.33999999999999997, 0.82879999999999987, 0.86),
                  (0.37119999999999997, 0.33999999999999997, 0.86),
                  (0.86, 0.33999999999999997, 0.82879999999999987)]
    
    # Load spectra
    spectra_objects = {}
    spectra_colors = {}
    # HACK data path
    data_dir = os.path.dirname(os.path.abspath(__file__))+"/data/spectra"
    for i, spectrum_filename in enumerate(spectra_filenames):
        fname = "{}/{}.fits".format(data_dir, spectrum_filename)
        if os.path.exists(fname):
            spectra_objects[spectrum_filename] = specutils.Spectrum1D.read(fname)
            spectra_colors[spectrum_filename] = hls_colors[i % len(hls_colors)]
        else:
            logger.warn("Could not find spectrum data file {}".format(fname))
        
    # Make figure object and axes
    if figure==None: 
        figure = plt.figure(figsize=(10,8))
    figure.subplots_adjust(left=0.10, right=0.95)
    figure.patch.set_facecolor("w")
    gs = GridSpec(4, 2)
    gs.update(top=0.77)
    ax_top_left = figure.add_subplot(gs.new_subplotspec((0, 0), rowspan=2))
    ax_top_right = figure.add_subplot(gs.new_subplotspec((0, 1), rowspan=2))
    ax_middle = figure.add_subplot(gs.new_subplotspec((2, 0), colspan=2))
    ax_bot_left = figure.add_subplot(gs.new_subplotspec((3, 0)))
    ax_bot_right = figure.add_subplot(gs.new_subplotspec((3, 1)))
    all_axes = [ax_top_left, ax_top_right, ax_middle, ax_bot_left, ax_bot_right]
    
    # HACK
    def _label_Teff(ax, spectra_filenames, Tlit, Teff):
        ii = np.argsort(Tlit)
        spectra_filenames = list(np.array(spectra_filenames)[ii])
        Tlit = list(np.array(Tlit)[ii])
        Teff = list(np.array(Teff)[ii])

        temp_lit = [[r"$T_{\rm lit}$", 'k']]
        temp_eff = [[r"$T_{\rm eff}$", 'k']]

        for spec_name, _Tlit, _Teff in zip(spectra_filenames, Tlit, Teff):
            temp_lit.append( [str(_Tlit), spectra_colors[spec_name]] )
            temp_eff.append( [str(_Teff), spectra_colors[spec_name]] )
            
        print( temp_lit)
        print( temp_eff)

        coord = 100
        for i,(entrylit, entryeff) in enumerate(zip(temp_lit, temp_eff)):
            ax.annotate(entrylit[0], xy=(4858.33, 0), xytext=(40, coord),
                        textcoords='offset points', size=10, color=entrylit[1])
            ax.annotate(entryeff[0], xy=(4858.33, 0), xytext=(270, coord),
                        textcoords='offset points', size=10, color=entryeff[1])
            coord -= 15
        
    # Plot
    for key, ax in zip(all_axes_names, all_axes):
        ax_dict = summary_figure[key]
        this_star, = ax.plot(normalized_spectrum.dispersion, normalized_spectrum.flux, 'k')
        for spec_name in ax_dict["spectra_to_plot"]:
            spec = spectra_objects[spec_name]
            ax.plot(spec.dispersion, spec.flux, color=spectra_colors[spec_name])
        
        ax.set_ylabel(ax_dict["label"])
        ax.set_xlim(ax_dict["wavelength_range"])
        ax.set_ylim(ax_dict["ylim"])
    
        formatter = plt.matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax.xaxis.set_major_formatter(formatter)
        
        if ax_dict["label"]=="H-beta":
            _label_Teff(ax,spectra_filenames,Tlit,Teff)

    handles = [this_star]
    labels = ["This star"]
    for i, spec_name in enumerate(spectra_filenames):
        color = spectra_colors[spec_name]
        artist = plt.Line2D((0,1),(0,0), color=color)
        handles.append(artist)
        labels.append(spec_name)
    figure.legend(handles, labels, (.6, .8), prop = {'size': 12},
                  ncol = 2)
    return figure

def make_snr_plot(normalized_spectrum, figure=None):
    if figure==None: 
        figure = plt.figure(figsize=(10,8))
    figure.subplots_adjust(left=0.10, right=0.95)
    figure.patch.set_facecolor("w")
    
    ax = figure.add_subplot(111)
    ax.plot(normalized_spectrum.dispersion, np.sqrt(normalized_spectrum.ivar), lw=1, color='k')
    ax.set_xlabel("Wavelength (A)")
    ax.set_ylabel("SNR = sqrt(ivar)")
    return figure


def make_synthesis_plot(plotdata, default_err=0.1,
                        plot_resid=True, plot_data_err=True,
                        xlim=None, ylim=(0,1.2), r_ylim=(-.1,.1),
                        xmajlocator=1, ymajlocator=.1, r_ymajlocator=.05,
                        xminlocator=.1, yminlocator=.1, r_yminlocator=.01,
                        fig=None, ax=None, ax_residual=None, figsize=None):
    if isinstance(plotdata, spectral_models.SpectralSynthesisModel):
        plotdata = plotdata.export_plot_data(default_err)
    elem, orig_p, logeps_err, model, data = plotdata
    abund0 = orig_p["log_eps({})".format(elem)]
    
    # Get proper plotting axes
    if ax is None:
        assert ax_residual is None
        if plot_resid:
            fig = plt.figure(figsize=figsize)
            gs = GridSpec(2, 1, height_ratios=[3,1])
            ax = fig.add_subplot(gs[0])
            ax_residual = fig.add_subplot(gs[1])
        else:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
    else:
        fig = ax.get_figure()
        if plot_resid:
            assert ax_residual is not None
    
    # Set up limits and labeling
    if xlim is None:
        xlim = (min(model["wl"][0], data["wl"][0]), max(model["wl"][-1], data["wl"][-1]))
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.xaxis.set_major_locator(MultipleLocator(xmajlocator))
    ax.xaxis.set_minor_locator(MultipleLocator(xminlocator))
    ax.yaxis.set_major_locator(MultipleLocator(ymajlocator))
    ax.yaxis.set_minor_locator(MultipleLocator(yminlocator))
    ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    if plot_resid:
        ax_residual.plot(xlim, [0,0], 'k:', lw=1)
        ax_residual.set_xlim(xlim)
        ax_residual.set_ylim(r_ylim)
        ax_residual.xaxis.set_major_locator(MultipleLocator(xmajlocator))
        ax_residual.xaxis.set_minor_locator(MultipleLocator(xminlocator))
        ax_residual.yaxis.set_major_locator(MultipleLocator(r_ymajlocator))
        ax_residual.yaxis.set_minor_locator(MultipleLocator(r_yminlocator))
        ax_residual.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
            
    # Plot no elem model
    noelem_kwargs={'color':'b','linestyle':':','lw':1}
    noelem_label="No {}".format(elem)
    ax.plot(model["wl"], model["f_none"], label=noelem_label, **noelem_kwargs)
    if plot_resid:
        ax_residual.plot(data["wl"], data["r_none"], **noelem_kwargs)
    
    # Plot data
    data_kwargs={'fmt':'o-', 'color':'k', 'ecolor':'k', 'label':None}
    if plot_data_err:
        yerr = data["err"]
    else:
        yerr = None
    ax.errorbar(data["wl"], data["flux"], yerr=yerr, **data_kwargs)
    
    # Plot +/- error models
    error_kwargs={'color':'r','linestyle':':','lw':1}
    for suffix in ["-err", "+err"]:
        ax.plot(model["wl"], np.array(model["f_"+suffix]), label=None, **error_kwargs)
        if plot_resid:
            ax_residual.plot(data["wl"], data["r_"+suffix], label=None, **error_kwargs)
    
    # Plot fit model
    fit_kwargs={'color':'r','linestyle':'-','lw':1}
    fit_label=r"A({})={:5.2f}$\pm${:.2f}".format(elem, abund0, logeps_err)
    ax.plot(model["wl"], model["f_fit"], label=fit_label, **fit_kwargs)
    if plot_resid:
        ax_residual.plot(data["wl"], data["r_fit"], **fit_kwargs)
    
    return fig
