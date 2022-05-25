#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import sys
import os
from PySide2 import (QtCore, QtGui as QtGui2, QtWidgets as QtGui)
import time
from six import iteritems
import numpy as np

import matplotlib
from matplotlib.ticker import MaxNLocator, MultipleLocator

import smh
from smh import utils, LineList
from smh.gui import mpl, style_utils
from smh.spectral_models import (ProfileFittingModel, SpectralSynthesisModel)

import logging
logger = logging.getLogger(__name__)
logger.addHandler(smh.handler)

if sys.platform == "darwin":
    # See http://successfulsoftware.net/2013/10/23/fixing-qt-4-for-mac-os-x-10-9-mavericks/
    substitutes = [
        (".Lucida Grande UI", "Lucida Grande"),
        (".Helvetica Neue DeskInterface", "Helvetica Neue")
    ]
    for substitute in substitutes:
        QtGui2.QFont.insertSubstitution(*substitute)

DOUBLE_CLICK_INTERVAL = 0.1 # MAGIC HACK
PICKER_TOLERANCE = 10 # MAGIC HACK
_QFONT = QtGui2.QFont("Helvetica Neue", 10)
_ROWHEIGHT = 20

## These are valid attrs of a spectral model
_allattrs = ["wavelength","expot","species","elements","loggf",
             "equivalent_width","equivalent_width_uncertainty",
             "reduced_equivalent_width",
             "abundances", "abundances_to_solar", "abundance_uncertainties",
             "is_acceptable", "is_upper_limit", "user_flag",
             "use_for_stellar_parameter_inference", "use_for_stellar_composition_inference",
             "measurement_type","fwhm"]
_labels = ["Wavelength $\lambda$",u"Excitation potential","Species","Element","log gf",
           "Equivalent width", "Equivalent width error",
           r"$\log{EW}/\lambda$",
           u"log ε","[X/H]", u"σ(log ε)",
           "Acceptable", "Upper Limit", "User Flag",
           "Use for Spectroscopic Stellar Parameters", "Use for Stellar Abundances",
           "Measurement Type", "FWHM"]
_short_labels = [u"λ",u"χ","ID","El.","loggf",
                 "EW", u"σ(EW)",
                 "REW",
                 "A(X)","[X/H]",u"σ(X)",
                 "", "ul", "flag",
                 "spflag", "abundflag","type", "FWHM"]
_formats = [":.1f",":.2f",":.1f","",":.2f",
            ":.1f",":.1f",
            ":.2f",
            ":.2f",":.2f",":.2f",
            "","","","",
            "","",":.2f"]
_dtypes = [float, float, float, object, float,
           float, float,
           float,
           float, float, float,
           bool, bool, bool,
           bool, bool,
           object, float]
_formats = ["{"+fmt+"}" for fmt in _formats]
_attr2label = dict(zip(_allattrs,_labels))
_attr2slabel = dict(zip(_allattrs,_short_labels))
_attr2format = dict(zip(_allattrs,_formats))
_attr2dtype = dict(zip(_allattrs,_dtypes))

class SMHSpecDisplay(mpl.MPLWidget):
    """
    Refactored class to display spectrum and residual plot.
    This holds a single spectral model and draws it.
    
    INTERACTIVE PORTIONS
    * Click/shift-click to mask/antimask
    * Rightclick to zoom
    * Shift-rightclick to pan [TODO]
    
    METHODS
    reset():
        Clear internal variables (besides session)
    new_session(session)
        Clear graphics and prepare for new session
        Also called if you make changes to the spectrum
    update_after_selection(selected_models):
        Call after a selection is updated (to change 
    update_after_measurement_change(changed_model):
        Call after a measurement is changed in changed_model
    update_selected_model()
        Retrieve selected model from parent object.
        (Requires parent._get_selected_model() unless you specify get_selected_model kwarg)
    set_selected_model(model):
        Force selected model to be model
    update_spectrum_figure(redraw=False):
        The main workhorse of this class
        Call when you want to update the spectrum figure

    update_mask_region():
        TODO used to update just the mask

    """
    def __init__(self, parent, session=None, 
                 get_selected_model=None,
                 enable_zoom=True, enable_masks=False,
                 enable_model_modifications=False,
                 label_ymin=1.0, label_ymax=1.2,
                 callbacks_after_fit=[],
                 comparison_spectrum=None,
                 **kwargs):
        super(SMHSpecDisplay, self).__init__(parent=parent, session=session, 
                                             **kwargs)
        self.parent = parent
        self.get_selected_model = get_selected_model
        self.callbacks_after_fit = callbacks_after_fit
        self.session = session
        self.label_ymin = label_ymin
        self.label_ymax = label_ymax
        
        self.comparison_spectrum = comparison_spectrum

        self.setMinimumSize(QtCore.QSize(100,100))
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Expanding)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)

        gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[1,2])
        gs.update(top=.95, bottom=.05, hspace=0)

        self.ax_residual = self.figure.add_subplot(gs[0])
        self.ax_residual.axhline(0, c="#666666")
        self.ax_residual.xaxis.set_major_locator(MaxNLocator(10))
        self.ax_residual.xaxis.get_major_formatter().set_useOffset(False)
        #self.ax_residual.yaxis.set_major_locator(MaxNLocator(2))
        #self.ax_residual.set_xticklabels([])
        self.ax_residual.set_ylabel("Resid")
        
        self.ax_spectrum = self.figure.add_subplot(gs[1], sharex=self.ax_residual)
        self.ax_spectrum.axhline(1, c='k', ls=':')
        #self.ax_spectrum.xaxis.get_major_formatter().set_useOffset(False)
        #self.ax_spectrum.xaxis.set_major_locator(MaxNLocator(5))
        self.ax_spectrum.set_xlabel(u"Wavelength (Å)")
        self.ax_spectrum.set_ylabel(r"Normalized flux")
        self.ax_spectrum.set_ylim(0, 1.2)
        self.ax_spectrum.yaxis.set_major_locator(MultipleLocator(0.2))
        self.ax_spectrum.yaxis.set_minor_locator(MultipleLocator(0.02))
        #self.ax_spectrum.set_yticks([0, 0.5, 1])
        
        if enable_zoom:
            self.enable_interactive_zoom()
            self.mpl_connect("key_press_event", self.key_press_zoom)
        if enable_masks:
            ## TODO this is not currently compatible but would be nice if it was
            #self.enable_drag_to_mask([self.ax_residual, self.ax_spectrum])
            self.mpl_connect("button_press_event", self.spectrum_left_mouse_press)
            self.mpl_connect("button_release_event", self.spectrum_left_mouse_release)
            ## Connect shift and space keys
            self.mpl_connect("key_press_event", self.key_press_flags)
            self.mpl_connect("key_release_event", self.key_release_flags)
        if enable_model_modifications:
            self.mpl_connect("key_press_event", self.key_press_model)
        self.setFocusPolicy(QtCore.Qt.ClickFocus)

        # Internal MPL variables
        self._lines = {
            "comparison_spectrum": self.ax_spectrum.plot([np.nan], [np.nan],
                                                         c="c", alpha=.5, drawstyle="steps-mid")[0],
            "spectrum": self.ax_spectrum.plot(
                [np.nan], [np.nan], c="k", drawstyle="steps-mid")[0], #None,
            "spectrum_fill": None,
            "residual_fill": None,
            "transitions_center_main": self.ax_spectrum.axvline(
                np.nan, c="#666666", linestyle=":"),
            "transitions_center_residual": self.ax_residual.axvline(
                np.nan, c="#666666", linestyle=":"),
            "linelabels": self.ax_spectrum.vlines(
                np.nan, np.nan, np.nan, color="blue", lw=1),
            "weak_linelabels": self.ax_spectrum.vlines(
                np.nan, np.nan, np.nan, color="blue", linestyle=':', lw=1),
            "model_masks": [],
            "nearby_lines": [],
            "model_fit": self.ax_spectrum.plot([np.nan], [np.nan], c="r")[0],
            "model_residual": self.ax_residual.plot(
                [np.nan], [np.nan], c="k", drawstyle="steps-mid")[0],
            "interactive_mask": [
                self.ax_spectrum.axvspan(xmin=np.nan, xmax=np.nan, ymin=np.nan,
                                         ymax=np.nan, facecolor="r", edgecolor="none", alpha=0.25,
                                         zorder=-5),
                self.ax_residual.axvspan(xmin=np.nan, xmax=np.nan, ymin=np.nan,
                                         ymax=np.nan, facecolor="r", edgecolor="none", alpha=0.25,
                                         zorder=-5)
            ]
        }
        self.reset()

    def sizeHint(self):
        return QtCore.QSize(125,100)
    def minimumSizeHint(self):
        return QtCore.QSize(10,10)
    def add_callback_after_fit(self, callback):
        self.callbacks_after_fit.append(callback)
    def reset_callback_after_fit(self, callback):
        self.callbacks_after_fit = []
    def reset(self):
        """
        Clear all internal variables (except session)
        """
        #logger.debug("Resetting Spectrum Figure ({})".format(self))
        self.selected_model = None
        
        if self.session is not None:
            drawstyle = self.session.setting(["plot_styles","spectrum_drawstyle"],"steps-mid")
            self._lines["spectrum"].set_drawstyle(drawstyle)
            self._lines["comparison_spectrum"].set_drawstyle(drawstyle)
        for key in ["spectrum", "transitions_center_main", "transitions_center_residual",
                    "model_fit", "model_residual"]:
            self._lines[key].set_data([np.nan],[np.nan])
        self.label_lines(None)
    def new_session(self, session):
        self.session = session
        self.reset()
        self.draw() #update_spectrum_figure()
        return None

    def update_comparison_spectrum(self, new_comparison_spectrum):
        self.comparison_spectrum = new_comparison_spectrum
        self.update_spectrum_figure(redraw=True)
        return None

    def update_after_selection(self, selected_models):
        if self.session is None: return None
        ## Use the last selected model as current model
        self.set_selected_model(selected_models[-1])
    def update_after_measurement_change(self, changed_model):
        if self.session is None: return None
        ## Only do something if the changed model is the current model
        if changed_model == self.selected_model:
            raise NotImplementedError
        return None
    
    def update_selected_model(self):
        if self.session is None: return None
        if self.get_selected_model is not None:
            self.selected_model = self.get_selected_model()
        else:
            self.selected_model = self.parent._get_selected_model()
        self._set_xlimits(self._get_current_xlimits())
        self.reset_zoom_limits()
    def set_selected_model(self, model):
        if self.session is None: return None
        self.selected_model = model
        self._set_xlimits(self._get_current_xlimits())
        self.reset_zoom_limits()
    def _get_current_xlimits(self, scale=1.02):
        if self.session is None: return None
        if self.selected_model is None: return None
        transitions = self.selected_model.transitions
        window = self.selected_model.metadata["window"]
        limits = [
            transitions["wavelength"][0] - window*scale,
            transitions["wavelength"][-1]+ window*scale,
        ]
        return limits
    def _set_xlimits(self, limits):
        self.ax_spectrum.set_xlim(limits)
        self.ax_residual.set_xlim(limits)
    def update_spectrum_figure(self, redraw=False, reset_limits=True,
                               label_transitions=None, label_rv=None):
        #logger.debug("update spectrum figure ({}, {}, {})".format(self, redraw, reset_limits))
        if label_transitions is not None: logger.info("labelling {} transitions for {} (rv={})".format(
                len(label_transitions), np.array(np.unique(label_transitions["species"])), label_rv))
        if self.session is None: return None
        if reset_limits:
            self.update_selected_model()
        #logger.debug(" selected model: {}".format(self.selected_model))
        if self.selected_model is None: return None
        selected_model = self.selected_model
        
        limits = self._get_current_xlimits()
        
        ## Plot spectrum and error bars
        success = self._plot_normalized_spectrum(limits)
        if not success: return None
        
        self._plot_comparison_spectrum(limits)
        
        ## Plot indication of current lines
        self._plot_current_lines(selected_model)

        ## Plot masks
        success = self._plot_masks()
        
        ## Plot model
        success = self._plot_model()

        ## Plot labeled lines
        self.label_lines(label_transitions, rv=label_rv)

        if redraw: self.draw()
        
        return None
    def label_lines(self, transitions, label_elem=False,
                    ymin=None, ymax=None,
                    rv=None,
                    strong_loggf_min = -3.0,
                    strong_expot_max = 2.0):
        """
        Draw vertical lines at transition wavelengths
        """
        if ymin is None: ymin = self.label_ymin
        if ymax is None: ymax = self.label_ymax
        if label_elem: #TODO add text saying wl and elem
            raise NotImplementedError
        if rv is None: rv = 0.0

        collection = self._lines["linelabels"]
        collection2 = self._lines["weak_linelabels"]
        if transitions is None:
            collection.set_paths([])
            collection2.set_paths([])
            return None
        assert isinstance(transitions, LineList), type(transitions)
        ii_strong = np.logical_and(transitions["expot"] < strong_expot_max,
                                   transitions["loggf"] > strong_loggf_min)
        xs = np.array(transitions["wavelength"]) * (1 + rv/299792.458)
        paths = []
        for x in xs[ii_strong]:
            paths.append([[x, ymin], [x, ymax]])
        collection.set_paths(paths)
        paths = []
        for x in xs[~ii_strong]:
            paths.append([[x, ymin], [x, ymax]])
        collection2.set_paths(paths)
        return None
        
    def key_press_model(self, event):
        selected_model = self.selected_model
        if selected_model is None: return None
        logger.debug("key_press_model: {}".format(event.key))
        key = event.key.lower()
        #if event.key not in "auf": return None
        if event.key == "a":
            selected_model.is_acceptable = True
        elif event.key == "u":
            selected_model.is_acceptable = False
        elif event.key == "f":
            selected_model.user_flag = (~selected_model.user_flag)
        return None
        
    def key_press_zoom(self, event):
        if event.key not in "1234": return None
        if self.session is None: return None
        ylim = self.session.setting(["zoom_shortcuts",int(event.key)],
                                    default_return_value=[0.0,1.2])
        self.ax_spectrum.set_ylim(ylim)
        self.draw()
        return None
	    
    def spectrum_left_mouse_press(self, event):
        """
        Listener for if mouse button pressed in spectrum or residual axis
        
        Masking and mask removal 
        """
        if self.session is None: return None

        if event.button != 1: return None
        if event.inaxes not in (self.ax_residual, self.ax_spectrum):
            return None
        
        selected_model = self.selected_model
        if selected_model is None: return None
        
        ## Doubleclick: remove mask, and if so refit/redraw
        if event.dblclick:
            mask_removed = self._remove_mask(selected_model, event)
            if mask_removed and isinstance(selected_model, ProfileFittingModel):
                selected_model.fit()
                self.update_spectrum_figure(True,False)
                for callback in self.callbacks_after_fit:
                    callback()
            return None

        ## Normal click: start drawing mask
        # Clear all masks if shift key state is not same as antimask_flag
        # Also change the antimask state
        ## TODO this is very sensitive, only want to do this if you click and drag. But hard to do.
        if selected_model.metadata["antimask_flag"] != self.shift_key_pressed:
            selected_model.metadata["mask"] = []
            selected_model.metadata["antimask_flag"] = not selected_model.metadata["antimask_flag"]
            logger.debug("Switching antimask flag to {}".format(selected_model.metadata["antimask_flag"]))
            # HACK
            #self.update_fitting_options()

        # Single click.
        xmin, xmax, ymin, ymax = (event.xdata, np.nan, -1e8, +1e8)
        for patch in self._lines["interactive_mask"]:
            patch.set_xy([
                [xmin, ymin],
                [xmin, ymax],
                [xmax, ymax],
                [xmax, ymin],
                [xmin, ymin]
            ])
            patch.set_facecolor("g" if selected_model.metadata["antimask_flag"] else "r")

        # Set the signal and the time.
        self._interactive_mask_region_signal = (
            time.time(),
            self.mpl_connect(
                "motion_notify_event", self.update_mask_region)
        )
        
        return None
    
    def update_mask_region(self, event):
        """
        Update the visible selected masked region for the selected spectral
        model. This function is linked to a callback for when the mouse position
        moves.

        :param event:
            The matplotlib motion event to show the current mouse position.
        """
        if self.session is None: return None

        if event.xdata is None: return

        signal_time, signal_cid = self._interactive_mask_region_signal
        if time.time() - signal_time > DOUBLE_CLICK_INTERVAL:
            data = self._lines["interactive_mask"][0].get_xy()
            # Update xmax.
            data[2:4, 0] = event.xdata
            for patch in self._lines["interactive_mask"]:
                patch.set_xy(data)
            self.draw()
        return None

    def spectrum_left_mouse_release(self, event):
        if self.session is None: return None
        try:
            signal_time, signal_cid = self._interactive_mask_region_signal
        except (AttributeError, TypeError):
            return None
        xy = self._lines["interactive_mask"][0].get_xy()
        if event.xdata is None:
            # Out of axis; exclude based on the closest axis limit
            xdata = xy[2, 0]
        else:
            xdata = event.xdata

        # If the two mouse events were within some time interval,
        # then we should not add a mask because those signals were probably
        # part of a double-click event.
        if  time.time() - signal_time > DOUBLE_CLICK_INTERVAL \
        and np.abs(xy[0,0] - xdata) > 0:
            
            spectral_model = self.selected_model

            if spectral_model is None: 
                raise RuntimeError("""Must have a spectral model selected while making mask!
                                   Must have mouseover bug?""")

            # Add mask metadata.
            minx = min(xy[0,0],xy[2,0])
            maxx = max(xy[0,0],xy[2,0])
            spectral_model.metadata["mask"].append([minx,maxx])

            # Re-fit the spectral model and send to other widgets.
            if isinstance(spectral_model, ProfileFittingModel):
                spectral_model.fit()
                for callback in self.callbacks_after_fit:
                    callback()

        # Clean up interactive mask
        xy[:, 0] = np.nan
        for patch in self._lines["interactive_mask"]:
            patch.set_xy(xy)
        self.mpl_disconnect(signal_cid)
        del self._interactive_mask_region_signal
        
        self.update_spectrum_figure(True,False)
        return None

    def _remove_mask(self, selected_model, event):
        if self.session is None: return None
        # Called upon doubleclick
        logger.debug("Removing mask")
        if selected_model.metadata["antimask_flag"]:
            ## Cannot remove antimasks properly
            logger.info("Cannot remove antimasks right now")
            return False
        else:
            for i, (s, e) in enumerate(selected_model.metadata["mask"][::-1]):
                if e >= event.xdata >= s:
                    mask = selected_model.metadata["mask"]
                    index = len(mask) - 1 - i
                    del mask[index]
                    return True
        return False

    def _plot_normalized_spectrum(self, limits, extra_disp=10):
        if self.session is None: return False
        if not hasattr(self.session, "normalized_spectrum"): return False

        spectrum = self.session.normalized_spectrum
        
        try:
            # Fix the memory leak!
            #self.ax_spectrum.lines.remove(self._lines["spectrum"])
            self._lines["spectrum"].set_data([np.nan], [np.nan])
            self.ax_spectrum.collections.remove(self._lines["spectrum_fill"])
            self.ax_residual.collections.remove(self._lines["residual_fill"])
        except Exception as e:
            # TODO fail in a better way
            # print(e)
            pass

        plot_ii = np.logical_and(spectrum.dispersion > limits[0]-extra_disp,
                                 spectrum.dispersion < limits[1]+extra_disp)
        if np.sum(plot_ii)==0: # Can't plot, no points!
            return False
        
        # Draw the spectrum.
        self._lines["spectrum"].set_data(spectrum.dispersion[plot_ii], spectrum.flux[plot_ii])

        # Draw the error bars.
        sigma = 1.0/np.sqrt(spectrum.ivar[plot_ii])
        self._lines["spectrum_fill"] = \
        style_utils.fill_between_steps(self.ax_spectrum, spectrum.dispersion[plot_ii],
            spectrum.flux[plot_ii] - sigma, spectrum.flux[plot_ii] + sigma, 
            facecolor="#cccccc", edgecolor="#cccccc", alpha=1)

        # Draw the error bars.
        self._lines["residual_fill"] = \
        style_utils.fill_between_steps(self.ax_residual, spectrum.dispersion[plot_ii],
            -sigma, +sigma, facecolor="#CCCCCC", edgecolor="none", alpha=1)

        three_sigma = 3*np.median(sigma[np.isfinite(sigma)])
        if np.isfinite(three_sigma):
            self.ax_residual.set_ylim(-three_sigma, three_sigma)
        
        return True
    
    def _plot_comparison_spectrum(self, limits, extra_disp=10):
        if self.comparison_spectrum is None: 
            self._lines["comparison_spectrum"].set_data([np.nan], [np.nan])
            return False            
        spectrum = self.comparison_spectrum
        
        plot_ii = np.logical_and(spectrum.dispersion > limits[0]-extra_disp,
                                 spectrum.dispersion < limits[1]+extra_disp)
        if np.sum(plot_ii)==0: # Can't plot, no points!
            self._lines["comparison_spectrum"].set_data([np.nan], [np.nan])
            return False            
        
        self._lines["comparison_spectrum"].set_data(
            spectrum.dispersion[plot_ii], spectrum.flux[plot_ii])
        return True
    
    def _plot_current_lines(self, selected_model):
        if isinstance(selected_model, ProfileFittingModel):
            x = selected_model.wavelength
        else:
            # TODO how to deal with syntheses?
            # TODO need to know 
            x = np.nan
        self._lines["transitions_center_main"].set_data([x, x], [0, 1.2])
        self._lines["transitions_center_residual"].set_data([x, x], [0, 1.2])
        
    def _plot_masks(self):
        if self.session is None: return False
        selected_model = self.selected_model
        mask_color = "g" if "antimask_flag" in selected_model.metadata and \
            selected_model.metadata["antimask_flag"] else "r"
        for i, (start, end) in enumerate(selected_model.metadata["mask"]):
            try:
                patches = self._lines["model_masks"][i]

            except IndexError:
                self._lines["model_masks"].append([
                    self.ax_spectrum.axvspan(np.nan, np.nan,
                        facecolor=mask_color, edgecolor="none", alpha=0.25),
                    self.ax_residual.axvspan(np.nan, np.nan,
                        facecolor=mask_color, edgecolor="none", alpha=0.25)
                ])
                patches = self._lines["model_masks"][-1]

            for patch in patches:
                patch.set_xy([
                    [start, -1e8],
                    [start, +1e8],
                    [end,   +1e8],
                    [end,   -1e8],
                    [start, -1e8]
                ])
                patch.set_facecolor(mask_color)
                patch.set_visible(True)

        # Hide unnecessary ones.
        N = len(selected_model.metadata["mask"])
        for unused_patches in self._lines["model_masks"][N:]:
            for unused_patch in unused_patches:
                unused_patch.set_visible(False)
        
        return True
        
    def _plot_model(self):
        if self.session is None: return False
        # Hide previous model_errs
        try:
            # Note: leaks memory, but too hard to fix
            self._lines["model_yerr"].set_visible(False)
            del self._lines["model_yerr"]
        except KeyError:
            pass
        
        selected_model = self.selected_model
        try:
            (named_p_opt, cov, meta) = selected_model.metadata["fitted_result"]

            # Test for some requirements.
            if "plot_x" not in meta: #Profile
                _ = (meta["model_x"], meta["model_y"], meta["residual"])
                plotxkey = "model_x"
                plotykey = "model_y"
            else: #Synthesis
                _ = (meta["plot_x"], meta["plot_y"], meta["residual"])
                plotxkey = "plot_x"
                plotykey = "plot_y"

        except KeyError:
            meta = {}
            self._lines["model_fit"].set_data([np.nan], [np.nan])
            self._lines["model_residual"].set_data([np.nan], [np.nan])

        else:
            assert len(meta[plotxkey]) == len(meta[plotykey])
            assert len(meta["model_x"]) == len(meta["residual"])
            assert len(meta["model_x"]) == len(meta["model_yerr"])

            self._lines["model_fit"].set_data(meta[plotxkey], meta[plotykey])
            self._lines["model_fit"].set_linestyle("-" if self.selected_model.is_acceptable else "--")
            self._lines["model_fit"].set_color("r" if self.selected_model.is_acceptable else "b")
            self._lines["model_residual"].set_data(meta["model_x"], meta["residual"])

            # Model yerr.
            if np.any(np.isfinite(meta["model_yerr"])):
                self._lines["model_yerr"] = self.ax_spectrum.fill_between(
                    meta["model_x"],
                    meta["model_y"] + meta["model_yerr"],
                    meta["model_y"] - meta["model_yerr"],
                    facecolor="r" if self.selected_model.is_acceptable else "b",
                    edgecolor="none", alpha=0.5)

            # Model masks due to nearby lines.
            if "nearby_lines" in meta:
                for i, (_, (start, end)) in enumerate(meta["nearby_lines"]):
                    try:
                        patches = self._lines["nearby_lines"][i]
                
                    except IndexError:
                        self._lines["nearby_lines"].append([
                            self.ax_spectrum.axvspan(np.nan, np.nan,
                                facecolor="b", edgecolor="none", alpha=0.25),
                            self.ax_residual.axvspan(np.nan, np.nan,
                                facecolor="b", edgecolor="none", alpha=0.25)
                        ])
                        patches = self._lines["nearby_lines"][-1]

                    for patch in patches:                            
                        patch.set_xy([
                            [start, -1e8],
                            [start, +1e8],
                            [end,   +1e8],
                            [end,   -1e8],
                            [start, -1e8]
                        ])
                        patch.set_visible(True)
                    
        # Hide any masked model regions due to nearby lines.
        N = len(meta.get("nearby_lines", []))
        for unused_patches in self._lines["nearby_lines"][N:]:
            for unused_patch in unused_patches:
                unused_patch.set_visible(False)

        return True
        
class SMHScatterplot(mpl.MPLWidget):
    """
    Scatterplot of things in a spectral model table.
    Gets all information using self.tablemodel.data(), so does not see the session.
    Interfaces with other things only by using the tableview to select or the tablemodel to set.
    
    SMHScatterplot.__init__:
        parent = parent of the widget
        xattr: name of x variable in the tablemodel to plot
        yattr: name of y variable in the tablemodel to plot
        tableview: the tableview that interacts with this plot [also used to access tablemodel]
        enable_zoom: allow zoom with rightclick
        enable_pick: allow point selection with leftclick
        enable_keyboard_shortcuts: enable ability to do things with selected points
        e_xattr, e_yattr: name of the error on x and y variable to plot, required if error_styles is not None
    The following things are used for setting styles and subselections.
    They are all lists, and must be the same length or None.
    Can put None into a list to skip that type of plot for that filter
        filters: list of additional spectral_model filters for substyling
        point_styles: list of plt.scatterplot keyword dictionaries
        error_styles: list of plt.errorbar keyword dictionaries (e.g. ecolor, elinewidth)
        linefit_styles: list of plt.plot keyword dictionaries (e.g. color, linestyle)
        linemean_styles: list of plt.plot keyword dictionaries (e.g. color, linestyle)
                         [note we actually plot median instead of mean, sorry for misnomer]

    Once you have set up __init__ properly, that sets what and how you want to plot.
    Then to update the figure, things you need to call:
    - linkToTable(tableview)
    - update_scatterplot()
    - update_selected_points()
    """
    allattrs = _allattrs
    labels = _labels
    attr2label = _attr2label
    def __init__(self, parent, xattr, yattr,
                 tableview=None,
                 enable_zoom=True, enable_pick=True,
                 enable_keyboard_shortcuts=True,
                 do_not_select_unacceptable=False,
                 ## These are settings for multiple things
                 exattr=None, eyattr=None,
                 filters=[lambda x: True],
                 point_styles=[{"s":30,"facecolor":"k","edgecolor":"k","alpha":0.5}],
                 error_styles=None,
                 linefit_styles=None,
                 linemean_styles=None,
                 **kwargs):
        assert xattr in self.allattrs, xattr
        assert yattr in self.allattrs, yattr
        self.xattr = xattr
        self.yattr = yattr
        if error_styles is not None:
            assert (e_xattr is not None) or (e_yattr is not None), "Must specify e_xattr and/or e_yattr for error_styles"
            assert (exattr is None) or (exattr in self.allattrs), exattr
            assert (eyattr is None) or (eyattr in self.allattrs), eyattr
        self.exattr = exattr
        self.eyattr = eyattr

        super(SMHScatterplot, self).__init__(parent=parent,
                                             **kwargs)
        
        self.parent = parent
        self.linkToTable(tableview)
        
        self.ax = self.figure.add_subplot(1,1,1)
        self.ax.xaxis.get_major_formatter().set_useOffset(False)
        self.ax.yaxis.set_major_locator(MaxNLocator(5))
        self.ax.yaxis.set_major_locator(MaxNLocator(4))
        self.ax.set_xlabel(self.attr2label[xattr])
        self.ax.set_ylabel(self.attr2label[yattr])

        ## Verify that all the style inputs are the right length
        assert len(filters) == len(point_styles)
        if error_styles is None:
            error_styles = [None for f in filters]
        else:
            assert len(filters)==len(error_styles)
            raise NotImplementedError("Need to implement error bar graphics updating!")
        if linefit_styles is None:
            linefit_styles = [None for f in filters]
        else:
            assert len(filters)==len(linefit_styles)
        if linemean_styles is None:
            linemean_styles = [None for f in filters]
        else:
            assert len(filters)==len(linemean_styles)
        
        ## Create graphic objects
        point_objs = []
        error_objs = []
        linefit_objs = []
        linemean_objs = []
        for filt, point_kw, error_kw, linefit_kw, linemean_kw in zip(
                filters, point_styles, error_styles, linefit_styles, linemean_styles):
            if point_kw is None: point_objs.append(None)
            else:
                point_objs.append(
                    self.ax.scatter([], [], picker=PICKER_TOLERANCE,
                                    **point_kw))
            
            if error_kw is None: error_objs.append(None)
            else:
                error_objs.append(
                    self.ax.errorbar(np.nan * np.ones(2), np.nan * np.ones(2),
                                     yerr=np.nan * np.ones((2, 2)),
                                     fmt=None,zorder=-10,
                                     **error_kw))
            
            if linefit_kw is None: linefit_objs.append(None)
            else:
                linefit_objs.append(
                    self.ax.plot([np.nan], [np.nan], **linefit_kw)[0])
            
            if linemean_kw is None: linemean_objs.append(None)
            else:
                linemean_objs.append(
                    self.ax.axhline(np.nan, **linemean_kw))
        
        ## Save graphic objects
        self._selected_points = self.ax.scatter([], [],
                                                edgecolor="b", facecolor="none", 
                                                s=150, linewidth=3, zorder=2)
        self._filters = filters
        self._points = point_objs
        self._errors = error_objs
        self._linefits = linefit_objs
        self._linemeans = linemean_objs
        self._graphics = list(zip(self._filters, self._points, self._errors, self._linefits, self._linemeans))
        
        ## Connect Interactivity
        if enable_zoom:
            self.enable_interactive_zoom()
        if enable_pick:
            self.canvas.callbacks.connect("pick_event", self.figure_mouse_pick)
        if enable_keyboard_shortcuts:
            self.mpl_connect("key_press_event", self.key_press_event)
        self.do_not_select_unacceptable = do_not_select_unacceptable
        self.setFocusPolicy(QtCore.Qt.ClickFocus)

        self.reset()
        self.update_scatterplot()
        self.update_selected_points(True)

    def sizeHint(self):
        return QtCore.QSize(125,100)
    def minimumSizeHint(self):
        return QtCore.QSize(10,10)
    def reset(self):
        self._selected_points.set_offsets(np.array([np.nan, np.nan]).T)
        for filt, point, error, linefit, linemean in self._graphics:
            if point is not None: point.set_offsets(np.array([np.nan, np.nan]).T)
            if error is not None: pass # TODO!!!
            if linefit is not None: linefit.set_data([np.nan],[np.nan])
            if linemean is not None: linemean.set_data([0,1],[np.nan,np.nan])
    def linkToTable(self, tableview):
        """
        view for selection; model for data
        """
        if tableview is None:
            self.tableview = None
            self.tablemodel= None
            self.xcol = None
            self.ycol = None
            return
        tablemodel = tableview.model()
        assert isinstance(tablemodel, MeasurementTableModelProxy) or \
               isinstance(tablemodel, MeasurementTableModelBase), \
               type(tablemodel)
        self.tableview = tableview
        self.tablemodel= tablemodel
        assert self.xattr in self.tablemodel.attrs, (self.xattr, self.tablemodel.attrs)
        assert self.yattr in self.tablemodel.attrs, (self.yattr, self.tablemodel.attrs)
        self.xcol = self.tablemodel.attrs.index(self.xattr)
        self.ycol = self.tablemodel.attrs.index(self.yattr)
        if self.exattr is None:
            self.excol = None
        else:
            assert self.exattr in self.tablemodel.attrs, (self.exattr, self.tablemodel.attrs)
            self.excol = self.tablemodel.attrs.index(self.exattr)
        if self.eyattr is None:
            self.eycol = None
        else:
            assert self.eyattr in self.tablemodel.attrs, (self.eyattr, self.tablemodel.attrs)
            self.eycol = self.tablemodel.attrs.index(self.eyattr)
        #logger.debug("Linked {} to {}/{}".format(self, self.tableview, self.tablemodel))
        #logger.debug("{}->{}, {}->{}".format(self.xattr, self.xcol, self.yattr, self.ycol))
        #if (self.exattr is not None) or (self.eyattr is not None):
        #    logger.debug("Err col: {}->{}, {}->{}".format(self.exattr, self.excol, self.eyattr, self.eycol))
    def _load_value_from_table(self, index):
        val = self.tablemodel.data(index, QtCore.Qt.DisplayRole)
        try:
            val = float(val)
        except ValueError as e:
            if val != "": logger.debug(e)
            val = np.nan
        return val
    def _get_spectral_models_from_rows(self, rows):
        ### TODO
        try:
            spectral_models = self.tablemodel().get_models_from_rows(rows)
            return spectral_models
        except:
            return None
    def update_scatterplot(self, redraw=False):
        if self.tableview is None or self.tablemodel is None: return None
        xs, ys, exs, eys = [], [], [], []
        Nrows = self.tablemodel.rowCount()
        if Nrows==0: return None
        spectral_models = self.tablemodel.get_models_from_rows(np.arange(Nrows))
        for i in range(Nrows):
            ix = self._ix(i, self.xcol)
            x = self._load_value_from_table(ix)
            ix = self._ix(i, self.ycol)
            y = self._load_value_from_table(ix)
            if self.excol is None: ex = np.nan
            else:
                ix = self._ix(i, self.excol)
                ex = self._load_value_from_table(ix)
            if self.eycol is None: ey = np.nan
            else:
                ix = self._ix(i, self.eycol)
                ey = self._load_value_from_table(ix)
            xs.append(x)
            ys.append(y)
            exs.append(ex)
            eys.append(ey)
        xs = np.array(xs); ys = np.array(ys); exs = np.array(exs); eys = np.array(eys)
        valids = []
        for filt in self._filters:
            valids.append([filt(sm) for sm in spectral_models])
        valids = np.atleast_2d(np.array(valids, dtype=bool))
        for ifilt,(filt, point, error, linefit, linemean) in enumerate(self._graphics):
            valid = valids[ifilt,:]
            nonzero = valid.sum() > 0
            x, y = xs[valid], ys[valid]
            if point is not None:
                if nonzero: point.set_offsets(np.array([x,y]).T)
                else: point.set_offsets(np.array([np.nan,np.nan]).T)
            if error is not None:
                ## TODO not doing anything with error bars right now
                pass
            if nonzero and ((linefit is not None) or (linemean is not None)):
                ## TODO: Figure out how best to save and return info about the lines
                ## For now, just refitting whenever needed
                try:
                    m, b, medy, stdy, stdm, N = utils.fit_line(x, y, None)
                except ValueError as e:
                    return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                #xlim = np.array(self.ax.get_xlim())
                xlim = np.array([x.min(), x.max()])
                if (linefit is not None) and nonzero:
                    linefit.set_data(xlim, m*xlim + b)
                else:
                    linefit.set_data([np.nan], [np.nan])
                if (linemean is not None) and nonzero:
                    linemean.set_data([0,1], [medy, medy])
                else:
                    linemean.set_data([0,1], [np.nan, np.nan])
        style_utils.relim_axes(self.ax)
        self.reset_zoom_limits()
        if redraw: self.draw()
        return None
    def update_selected_points(self, redraw=False):
        if self.tableview is None or self.tablemodel is None: return None
        #logger.debug("update_selected_points ({}, {})".format(self, redraw))
        xs = []
        ys = []
        rows = np.array([row.row() for row in self.tableview.selectionModel().selectedRows()])
        if self.do_not_select_unacceptable and len(rows)>0:
            not_acceptable = np.logical_not(self.tablemodel.get_data_column("is_acceptable"))
            is_upper_limit = self.tablemodel.get_data_column("is_upper_limit")
            skip_plot = not_acceptable | is_upper_limit
        for i in rows:
            if self.do_not_select_unacceptable and skip_plot[i]:
                continue
            ix = self._ix(i, self.xcol)
            x = self._load_value_from_table(ix)
            ix = self._ix(i, self.ycol)
            y = self._load_value_from_table(ix)
            xs.append(x)
            ys.append(y)
        self._selected_points.set_offsets(np.array([xs,ys]).T)
        if redraw: self.draw()
        return None

    def _ix(self,row,col):
        if self.tablemodel is None: return None
        return self.tablemodel.createIndex(row,col)

    def figure_mouse_pick(self, event):
        if event.mouseevent.button != 1: return None
        ## this is fast but broken when multiple points are on the scatterplot.
        #self.tableview.selectRow(event.ind[0])
        xscale = np.ptp(self.ax.get_xlim())
        yscale = np.ptp(self.ax.get_ylim())
        xall = self.tablemodel.get_data_column(self.xcol)
        yall = self.tablemodel.get_data_column(self.ycol)
        xpick = event.mouseevent.xdata
        ypick = event.mouseevent.ydata
        dist = np.sqrt((xall-xpick)**2 + (yall-ypick)**2)
        ix = np.nanargmin(dist)
        #logger.debug("picked row {} with {} {}".format(ix, xpick, ypick))
        self.tableview.selectRow(ix)
        return None
    
    def key_press_event(self, event):
        if event.key not in "uUaA": return None
        if self.tableview is None or self.tablemodel is None: return None
        if event.key in "uU":
            self.mark_selected_models_as_unacceptable()
        elif event.key in "aA":
            self.mark_selected_models_as_acceptable()
    def mark_selected_models_as_unacceptable(self):
        col = self.tablemodel.attrs.index("is_acceptable")
        rows = self.tableview.selectionModel().selectedRows()
        for row in rows:
            i = row.row()
            ix = self._ix(i, col)
            self.tablemodel.setData(ix, False)
            self.tableview.update_row(i)
    def mark_selected_models_as_acceptable(self):
        col = self.tablemodel.attrs.index("is_acceptable")
        rows = self.tableview.selectionModel().selectedRows()
        for row in rows:
            i = row.row()
            ix = self._ix(i, col)
            self.tablemodel.setData(ix, True)
            self.tableview.update_row(i)
    
class BaseTableView(QtGui.QTableView):
    """ Basic sizing and options for display table """
    def __init__(self, parent, *args):
        super(BaseTableView, self).__init__(parent, *args)
        self.parent = parent
        self.setSortingEnabled(False)
        self.verticalHeader().setDefaultSectionSize(_ROWHEIGHT)
        self.horizontalHeader().setStretchLastSection(True)
        self.horizontalHeader().setSectionResizeMode(QtGui.QHeaderView.Stretch)
        self.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        
    def sizeHint(self):
        return QtCore.QSize(125,100)
    def minimumSizeHint(self):
        return QtCore.QSize(125,100)
class MeasurementTableView(BaseTableView):
    def __init__(self, parent, session,
                 callbacks_after_menu=[], display_fitting_options=True):
        super(MeasurementTableView, self).__init__(parent)
        self.session = session
        self.callbacks_after_menu = callbacks_after_menu
        self.display_fitting_options = display_fitting_options
    def update_session(self, session):
        self.session = session
        self.model().reindex()
    def update_row(self,row):
        """ Used for proxy models to efficiently update data"""
        self.rowMoved(row, row, row)
        return None
    def menu_finished(self):
        """ Use to refresh GUI after changing fitting options """
        for callback in self.callbacks_after_menu:
            callback()
    #############################
    # Stuff for rightclick menu #
    #############################
    def contextMenuEvent(self, event):
        indices = self.selectionModel().selectedRows()
        rows = [index.row() for index in indices]
        N = len(rows)
        
        ### Create rightclick menu
        menu = QtGui.QMenu(self)
        fit_selected_models = menu.addAction(
            "Fit selected model{}..".format(["", "s"][N != 1]))
        measure_selected_models = menu.addAction(
            "Measure selected model{}..".format(["", "s"][N != 1]))

        menu.addSeparator()

        mark_as_acceptable = menu.addAction("Mark as acceptable")
        mark_as_unacceptable = menu.addAction("Mark as unacceptable")

        set_user_flag = menu.addAction("Set user flag")
        unset_user_flag = menu.addAction("Unset user flag")

        if self.display_fitting_options:
            menu.addSeparator()
            
            set_fitting_window = menu.addAction("Set fitting window..")
            continuum_menu = menu.addMenu("Set continuum")
            set_no_continuum = continuum_menu.addAction("No continuum",
                checkable=True)
            continuum_menu.addSeparator()
            set_continuum_order = [continuum_menu.addAction(
                "Order {}".format(i), checkable=True) for i in range(0, 10)]
    
            menu.addSeparator()
    
            menu_profile_type = menu.addMenu("Set profile type")
    
            set_gaussian = menu_profile_type.addAction("Gaussian")
            set_lorentzian = menu_profile_type.addAction("Lorentzian")
            set_voigt = menu_profile_type.addAction("Voigt")
    
            menu.addSeparator()
    
            enable_central_weighting = menu.addAction("Enable central weighting")
            disable_central_weighting = menu.addAction("Disable central weighting")
    
            menu.addSeparator()
    
            set_detection_sigma = menu.addAction("Set detection sigma..")
            set_detection_pixels = menu.addAction("Set detection pixels..")
    
            menu.addSeparator()
    
            set_rv_tolerance = menu.addAction("Set RV tolerance..")
            set_wl_tolerance = menu.addAction("Set WL tolerance..")

        if N == 0:
            menu.setEnabled(False)
            
        action = menu.exec_(self.mapToGlobal(event.pos()))
    
        if action == fit_selected_models:
            self.fit_selected_models()
        elif action == measure_selected_models:
            self.measure_selected_models()
        elif action in (mark_as_acceptable, mark_as_unacceptable):
            toggle = (action == mark_as_acceptable)
            self.set_flag("is_acceptable", toggle)
        elif action in (set_user_flag, unset_user_flag):
            toggle = (action == set_user_flag)
            self.set_flag("user_flag", toggle)
        if self.display_fitting_options:
            if action == set_fitting_window:
                self.set_fitting_window()
            elif action == set_no_continuum:
                self.set_continuum_order(-1)
            elif action == set_continuum_order:
                order = set_continuum_order.index(action)
                self.set_continuum_order(order)
            elif action in (set_gaussian, set_lorentzian, set_voigt):
                kind = {
                    set_gaussian: "gaussian",
                    set_lorentzian: "lorentzian",
                    set_voigt: "voigt"
                }[action]
                self.set_profile(kind)
            elif action in (enable_central_weighting, disable_central_weighting):
                toggle = (action == enable_central_weighting)
                self.set_central_weighting(toggle)
            elif action == set_detection_sigma:
                self.set_detection_sigma()
            elif action == set_detection_pixels:
                self.set_detection_pixels()
            elif action == set_rv_tolerance:
                self.set_rv_tolerance()
            elif action == self.set_wl_tolerance:
                self.set_wl_tolerance()
        
        self.menu_finished()
        return None

    def get_selected_models(self, getrows=False):
        indices = self.selectionModel().selectedRows()
        rows = [index.row() for index in indices]
        spectral_models = self.model().get_models_from_rows(rows)
        if getrows:
            return rows, spectral_models
        return spectral_models

    def fit_selected_models(self):
        rows, spectral_models = self.get_selected_models(getrows=True)
        for row, model in zip(rows, spectral_models):
            try:
                model.fit()
            except:
                logger.exception("Error fitting row {}".format(row))
            self.update_row(row)
        return None
    
    def measure_selected_models(self):
        rows, spectral_models = self.get_selected_models(getrows=True)
        self.session.measure_abundances(spectral_models)
        for row in rows:
            self.update_row(row)
        return None
    
    def set_flag(self, flag_field, toggle):
        assert flag_field in ["is_acceptable", "is_upper_limit", "user_flag"], flag_field
        value = 2 if toggle else 0
        rows, spectral_models = self.get_selected_models(getrows=True)
        data_model = self.model()
        col = self.model().attrs.index(flag_field)
        for row in rows:
            index = data_model.createIndex(row, col)
            data_model.setData(index, value)
        return None
    
    # E. Holmbeck added extra keypress for convenience.
    # Press the leftarrow to deselect and right arrow to select a line
    def keyPressEvent (self, eventQKeyEvent):
        key = eventQKeyEvent.key()
        if key == QtCore.Qt.Key_Left:
            self.set_flag("is_acceptable", False)
        elif key == QtCore.Qt.Key_Right:
            self.set_flag("is_acceptable", True)
        else:
            super(MeasurementTableView, self).keyPressEvent(eventQKeyEvent)
    
      
    def set_fitting_option_value(self, key, value,
                                 valid_for_profile=False,
                                 valid_for_synth=False):
        rows, spectral_models = self.get_selected_models(getrows=True)
        num_fit = 0
        num_unacceptable = 0
        num_profile_models = 0
        num_synthesis_models = 0
        num_error = 0
        for row, spectral_model in zip(rows, spectral_models):
            run_fit = False
            if not spectral_model.is_acceptable: 
                num_unacceptable += 1
                continue
            if valid_for_profile and isinstance(spectral_model,ProfileFittingModel):
                num_profile_models += 1
                spectral_model.metadata[key] = value
                run_fit = True
            if valid_for_synth and isinstance(spectral_model,SpectralSynthesisModel):
                num_synthesis_models += 1
                spectral_model.metadata[key] = value
                run_fit = True
                
            if run_fit and "fitted_result" in spectral_model.metadata:
                num_fit += 1
                try:
                    spectral_model.fit()
                except:
                    num_error += 1
                    logger.exception("Error fitting row {} after modifying {} to {}".format(row, key, value))
                self.update_row(row)
        logger.info("Changed {0}={1}, fit {2} out of {3} models ({4} profile, {5} synth, {6} unacceptable, {7} fit fail)".format(\
                key, value, num_fit, len(spectral_models), num_profile_models, num_synthesis_models, num_unacceptable, num_error))
        return None
    def set_fitting_window(self):
        window, is_ok = QtGui.QInputDialog.getDouble(
            None, "Set fitting window", u"Fitting window (Å):", 
            value=5, minValue=0.1, maxValue=1000)
        if not is_ok: return None
        self.set_fitting_option_value("window", window,
                                      valid_for_profile=True,
                                      valid_for_synth=True)
        return None
    def set_continuum_order(self, order):
        self.set_fitting_option_value("continuum_order", order,
                                      valid_for_profile=True,
                                      valid_for_synth=True)
        return None
    def set_profile(self, kind):
        self.set_fitting_option_value("profile", kind,
                                      valid_for_profile=True,
                                      valid_for_synth=False)
        return None
    def set_central_weighting(self, toggle):
        self.set_fitting_option_value("central_weighting", toggle,
                                      valid_for_profile=True,
                                      valid_for_synth=False)
        return None
    def set_detection_sigma(self):
        detection_sigma, is_ok = QtGui.QInputDialog.getDouble(
            None, "Set detection sigma", u"Detection sigma:", 
            value=0.5, minValue=0.1, maxValue=1000)
        if not is_ok: return None
        self.set_fitting_option_value("detection_sigma", detection_sigma,
                                      valid_for_profile=True,
                                      valid_for_synth=False)
    def set_detection_pixels(self):
        detection_pixels, is_ok = QtGui.QInputDialog.getInt(
            None, "Set detection pixel", u"Detection pixels:", 
            value=3, minValue=1, maxValue=1000)
        if not is_ok: return None
        self.set_fitting_option_value("detection_pixels", detection_pixels,
                                      valid_for_profile=True,
                                      valid_for_synth=False)
    def set_rv_tolerance(self):
        velocity_tolerance, is_ok = QtGui.QInputDialog.getDouble(
            None, "Set velocity tolerance", u"Velocity tolerance:", 
            value=5, minValue=0.01, maxValue=100)
        if not is_ok: return None
        # TODO cannot turn it back to None right now
        self.set_fitting_option_value("velocity_tolerance", velocity_tolerance,
                                      valid_for_profile=True,
                                      valid_for_synth=True)
    def set_wl_tolerance(self):
        wavelength_tolerance, is_ok = QtGui.QInputDialog.getDouble(
            None, "Set wavelength tolerance", u"Wavelength tolerance:", 
            value=0.1, minValue=0.01, maxValue=1)
        if not is_ok: return None
        # TODO cannot turn it back to None right now
        self.set_fitting_option_value("wavelength_tolerance", wavelength_tolerance,
                                      valid_for_profile=True,
                                      valid_for_synth=False)

def create_measurement_table_with_buttons(parent, filtermodel, session, **kwargs):
    """
    kwargs go to MeasurementTableView (esp. callbacks_after_menu, display_fitting_options)
    """
    vbox = QtGui.QVBoxLayout()
    tableview = MeasurementTableView(parent, session, **kwargs)
    tableview.setModel(filtermodel)
    vbox.addWidget(tableview)
    
    hbox = QtGui.QHBoxLayout()

    btn_filter = QtGui.QPushButton(parent)
    btn_filter.setText("Hide unacceptable")
    show_or_hide_unacceptable = lambda: filtermodel.show_or_hide_unacceptable(btn_filter)
    btn_filter.clicked.connect(show_or_hide_unacceptable)
    
    btn_refresh = QtGui.QPushButton(parent)
    btn_refresh.setText("Refresh")

    hbox.addWidget(btn_filter)
    hbox.addWidget(btn_refresh)
    
    vbox.addLayout(hbox)
    return vbox, tableview, btn_filter, btn_refresh
    
class MeasurementTableDelegate(QtGui.QItemDelegate):
    ## TODO this doesn't work
    ## It doesn't paint checkboxes or get the font right anymore
    COLOR = "#FFFF00" #yellow
    def __init__(self, parent, view, *args):
        super(MeasurementTableDelegate, self).__init__(parent, *args)
        self.view = view
    def paint(self, painter, option, index):
        painter.save()
        painter.setPen(QtGui2.QPen(QtCore.Qt.NoPen))
        if option.state & QtGui.QStyle.State_Selected:
            painter.setBrush(QtGui2.QBrush(
                self.parent().palette().highlight().color()))
        else:
            model = self.view.model()
            if model is not None and "user_flag" in model.attrs:
                row = index.row()
                col = model.attrs.index("user_flag")
                state = model.data(model.createIndex(row,col))
                if state == QtCore.Qt.Checked:
                    painter.setBrush(QtGui2.QBrush(QtGui2.QColor(self.COLOR)))
                else:
                    painter.setBrush(QtGui2.QBrush(QtCore.Qt.white))
            else:
                painter.setBrush(QtGui2.QBrush(QtCore.Qt.white))
        painter.drawRect(option.rect)
        painter.setPen(QtGui2.QPen(QtCore.Qt.black))
        painter.drawText(option.rect, QtCore.Qt.AlignLeft|QtCore.Qt.AlignCenter, index.data())
        painter.restore()
        

class MeasurementTableModelProxy(QtCore.QSortFilterProxyModel):
    """
    Proxy model allowing for filtering (and eventually sorting) of the full MeasurementTableModelBase
    Based on the old SpectralModelsFilterProxyModel
    """
    def __init__(self, parent=None, views_to_update=[], callbacks_after_setData=[]):
        """
        Views to update must implement update_row(proxy_index).
        """
        super(MeasurementTableModelProxy, self).__init__(parent)
        self.filter_functions = {}
        self.views_to_update = views_to_update
        self.callbacks_after_setData = callbacks_after_setData
        return None
    @property
    def attrs(self):
        try:
            return self.sourceModel().attrs
        except:
            return []
    def show_or_hide_unacceptable(self, btn):
        hide = btn.text().startswith("Hide")
        if hide:
            self.add_filter_function(
                "is_acceptable", lambda model: model.is_acceptable)
        else:
            self.delete_filter_function("is_acceptable")
        text = "{} unacceptable".format(("Hide","Show")[hide])
        btn.setText(text)
        return None
    def add_view_to_update(self, view):
        self.views_to_update.append(view)
    def reset_views_to_update(self):
        self.views_to_update = []
    def add_callback_after_setData(self, callback):
        self.callbacks_after_setData.append(callback)
    def reset_callbacks_after_setData(self):
        self.callbacks_after_setData = []
    def get_data_column(self, column, rows=None):
        """ Function to quickly go under the hood and access one column """
        if rows is None:
            rows = self.lookup_indices
        else:
            rows = self.lookup_indices[np.array(rows)]
        data = self.sourceModel().get_data_column(column, rows=rows)
        return data
    def setData(self, proxy_index, value, role=QtCore.Qt.DisplayRole):
        """
        Only allow checking/unchecking of is_acceptable and user_flag
        """
        attr = self.sourceModel().attrs[proxy_index.column()]
        if attr not in ["is_acceptable", "is_upper_limit", "user_flag"]:
            return False
        else:
            proxy_row = proxy_index.row()
            data_row = self.mapToSource(proxy_index).row()
            model = self.sourceModel().spectral_models[data_row]
            # value appears to be 0 or 2. Set it to True or False
            value = (value != 0)
            setattr(model, attr, value)
            
            for view in self.views_to_update:
                view.update_row(proxy_row)
            for callback in self.callbacks_after_setData:
                callback()
            return value
    def add_filter_function(self, name, filter_function):
        self.filter_functions[name] = filter_function
        self.invalidateFilter()
        self.reindex()
        return None
    def delete_filter_function(self, name):
        try:
            del self.filter_functions[name]
            self.invalidateFilter()
            self.reindex()
        except KeyError:
            raise
        else:
            return None
    def delete_all_filter_functions(self):
        self.filter_functions = {}
        self.invalidateFilter()
        self.reindex()
        return None
    def reset(self, *args):
        #super(MeasurementTableModelProxy, self).reset(*args)
        self.beginResetModel()
        self.reindex()
        self.endResetModel()
        return None
    def reindex(self):
        try: 
            self.sourceModel().spectral_models
        except AttributeError:
            return None

        lookup_indices = []
        for i, model in enumerate(self.sourceModel().spectral_models):
            for name, filter_function in self.filter_functions.items():
                if not filter_function(model):
                    break
            else:
                # No problems with any filter functions.
                lookup_indices.append(i)

        self.lookup_indices = np.array(lookup_indices)
        return None
    def filterAcceptsRow(self, row, parent):
        # Check if we need to update the filter indices for mapping.
        model = self.sourceModel().spectral_models[row]
        for filter_name, filter_function in self.filter_functions.items():
            if not filter_function(model): break
        else:
            # No problems.
            return True

        # We broke out of the for loop.
        return False
    def mapFromSource(self, data_index):
        if not data_index.isValid():
            return data_index
        #self.reindex()
        return self.createIndex(
            np.where(self.lookup_indices == data_index.row())[0],
            data_index.column())
    def mapToSource(self, proxy_index):
        if not proxy_index.isValid():
            return proxy_index
        try:
            #self.reindex()
            return self.createIndex(self.lookup_indices[proxy_index.row()],
                proxy_index.column())
        except AttributeError:
            return proxy_index
        except IndexError as e:
            print("INDEX ERROR IN TABLE: probably you loaded a new file while 'hide unacceptable' was activated")
            print(e)
    def get_models_from_rows(self, rows):
        actual_rows = [self.lookup_indices[row] for row in rows]
        return self.sourceModel().get_models_from_rows(actual_rows)

class MeasurementTableModelBase(QtCore.QAbstractTableModel):
    """
    A table model that accesses the properties of a list of spectral models
    (calling them "measurements" here to avoid ambiguity)
    Based on the old SpectralModelsTableModel
    """
    allattrs = _allattrs
    attr2slabel = _attr2slabel
    attr2format = _attr2format
    
    def __init__(self, parent, session, columns, *args):
        """
        An abstract table model for accessing spectral models from the session.
        
        :param parent:
        :param session:
        :param columns:
            List of attributes in order of columns in the table.
        """
        
        super(MeasurementTableModelBase, self).__init__(parent, *args)
        self.verify_columns(columns)
        self.attrs = columns
        self.header = [self.attr2slabel[attr] for attr in self.attrs]
        
        # Normally you should never do this, but here I know "better". See:
        #http://stackoverflow.com/questions/867938/qabstractitemmodel-parent-why
        self.parent = parent 
        self.session = session
        return None
    
    def new_session(self, session):
        self.beginResetModel()
        self.session = session
        self.verify_columns(self.attrs)
        self.endResetModel()
        return None

    def verify_columns(self, columns):
        for col in columns:
            assert col in self.allattrs, col
        assert len(columns) == len(np.unique(columns))
        return True

    def get_data_column(self, column, rows=None):
        """ Function to quickly go under the hood and access one column """
        if isinstance(column, int):
            attr = self.attrs[column]
        else:
            attr = column
        models = self.spectral_models
        if rows is None:
            rows = np.arange(len(models))
        getter = lambda ix: getattr(models[ix], attr, np.nan)
        data = np.array([np.ravel(getter(r)) for r in rows], dtype=_attr2dtype[attr])
        return data

    def get_models_from_rows(self, rows):
        models_to_return = []
        for row in rows:
            models_to_return.append(self.spectral_models[row])
        return models_to_return

    def data(self, index, role):
        """
        Retrieve data from spectral models and turn into string
        """
        if not index.isValid():
            return None
        if role==QtCore.Qt.FontRole:
            return _QFONT
        attr = self.attrs[index.column()]
        spectral_model = self.spectral_models[index.row()]
        
        value = getattr(spectral_model, attr, None)
        if value is None: return ""
        
        ## Deal with checkboxies
        if attr in ["is_acceptable","is_upper_limit","user_flag",
                    "use_for_stellar_parameter_inference",
                    "use_for_stellar_composition_inference"] \
        and role in (QtCore.Qt.DisplayRole, QtCore.Qt.CheckStateRole):
            if role == QtCore.Qt.CheckStateRole:
                return QtCore.Qt.Checked if value else QtCore.Qt.Unchecked
            else:
                return None
        
        if role != QtCore.Qt.DisplayRole: return None
        
        # any specific hacks to display attrs are put in here
        # TODO the spectral model itself should decide how its own information is displayed.
        # Make a _repr_attr for every attr?
        if attr == "wavelength":
            try:
                return str(float(spectral_model.wavelength))
            except:
                return spectral_model._repr_wavelength
        if attr == "elements":
            return spectral_model._repr_element

        fmt = self.attr2format[attr]
        if isinstance(value, (list, np.ndarray)):
            try: 
                if isinstance(value[0], (list, np.ndarray)): # list of lists for syntheses species
                    mystrs = [[fmt.format(v) for v in vlist] for vlist in value]
                    mystrs = [item for sublist in mystrs for item in sublist]
                else:
                    mystrs = [fmt.format(v) for v in value]
            except ValueError as e:
                logger.debug("A: {} has fmt {} and failing on {}".format(attr, fmt, value))
                mystr = str(value)
            else:
                mystr = ";".join(mystrs)
        else:
            try:
                mystr = fmt.format(value)
            except ValueError as e:
                logger.debug("B: {} has fmt {} and failing on {}".format(attr, fmt, value))
                mystr = str(value)
        return mystr

    @property
    def spectral_models(self):
        if self.session is None: return []
        return self.session.metadata.get("spectral_models", [])
    
    def rowCount(self, parent=None):
        """ Return the number of rows in the table. """
        return len(self.spectral_models)
    
    def columnCount(self, parent=None):
        """ Return the number of columns in the table. """
        return len(self.header)

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal \
        and role == QtCore.Qt.DisplayRole:
            return self.header[col]
        return None

    def setData(self, index, value, role=QtCore.Qt.DisplayRole):
        """
        Only allow checking/unchecking of is_acceptable and user_flag
        """
        attr = self.attrs[index.column()]
        if attr not in ["is_acceptable", "is_upper_limit", "user_flag"]:
            return False
        else:
            row = index.row()
            model = self.spectral_models[row]
            # value appears to be 0 or 2. Set it to True or False
            value = (value != 0)
            setattr(model, attr, value)
            ## TODO this emit is SUPER slow with proxy models.
            ## You should overwrite setData for the proxy model
            ## and explicitly connect it to the view.
            self.dataChanged.emit(self.createIndex(row, 0),
                                  self.createIndex(row, 
                                  self.columnCount(None)))
            return value

    def flags(self, index):
        if not index.isValid(): return
        return  QtCore.Qt.ItemIsSelectable|\
                QtCore.Qt.ItemIsEnabled|\
                QtCore.Qt.ItemIsUserCheckable
    
    
class MeasurementSummaryTableModel(QtCore.QAbstractTableModel):
    header = ["El.", "Species", "N", "A(X)", "σ(X)", "[X/H]", "[X/Fe]"]
    def __init__(self, parent, session, *args):
        super(MeasurementSummaryTableModel, self).__init__(parent, *args)
        self.parent = parent 
        self.new_session(session)
        return None
    def summarize(self):
        if self.session is None:
            self.summary = {}
        else:
            self.summary = self.session.summarize_spectral_models(what_fe=self.what_fe)
    def new_session(self, session):
        """
        Reset the table based on the new session
        """
        self.beginResetModel()
        self.session = session
        if session is None:
            self.what_fe = 1
        else:
            self.what_fe = session.setting("what_fe", 1)
        self.summarize()
        self.all_species = np.sort(list(self.summary.keys()))
        self.endResetModel()
    def update_summary(self, species=None):
        """
        Call this when you update any measurement
        """
        if species is None:
            # Reset the whole model
            self.new_session(self.session)
            return None
        row = np.where(species == self.all_species)[0]
        if len(row) != 1:
            logger.debug("Invalid species: {} ({})".format(species, self.all_species))
            logger.debug("Resetting whole model")
            self.new_session(self.session)
            return None
        # Summarize everything, but only emit that you changed one row
        row = row[0]
        self.summarize()
        self.dataChanged.emit(self.createIndex(row, 0),
                              self.createIndex(row, self.columnCount()))
        return None
    def rowCount(self, parent=None):
        return len(self.all_species)
    def columnCount(self, parent=None):
        return len(self.header)
    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal \
        and role == QtCore.Qt.DisplayRole:
            return self.header[col]
        return None
    def flags(self, index):
        if not index.isValid(): return
        return  QtCore.Qt.ItemIsSelectable|\
                QtCore.Qt.ItemIsEnabled
    def data(self, index, role):
        if not index.isValid(): return None
        if role==QtCore.Qt.FontRole: return _QFONT
        if role != QtCore.Qt.DisplayRole: return None

        row = index.row()
        species = self.all_species[row]
        vals = self.summary[species]
        num_models, logeps, stdev, stderr, XH, XFe = vals
        
        col = index.column()
        if col == 0:
            return utils.species_to_element(species)
        elif col == 1:
            return str(species)
        elif col == 2:
            return str(num_models)
        elif col == 3:
            return "{:.2f}".format(logeps)
        elif col == 4:
            return "{:.2f}".format(stdev)
        elif col == 5:
            return "{:.2f}".format(XH)
        elif col == 6:
            return "{:.2f}".format(XFe)
        raise ValueError("row={} col={} species={}".format(row, col, species))
