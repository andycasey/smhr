#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The stellar parameters tab in Spectroscopy Made Hard """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["ChemicalAbundancesTab"]

import matplotlib.gridspec
import numpy as np
import sys
from PySide import QtCore, QtGui
import time

import mpl, style_utils
from matplotlib.ticker import MaxNLocator
#from smh.photospheres import available as available_photospheres
import smh.radiative_transfer as rt
from smh.spectral_models import (ProfileFittingModel, SpectralSynthesisModel)
from abund_tree import AbundTreeView, AbundTreeModel, AbundTreeMeasurementItem, AbundTreeElementSummaryItem
from linelist_manager import TransitionsDialog

import logging
logger = logging.getLogger(__name__)

if sys.platform == "darwin":
        
    # See http://successfulsoftware.net/2013/10/23/fixing-qt-4-for-mac-os-x-10-9-mavericks/
    substitutes = [
        (".Lucida Grande UI", "Lucida Grande"),
        (".Helvetica Neue DeskInterface", "Helvetica Neue")
    ]
    for substitute in substitutes:
        QtGui.QFont.insertSubstitution(*substitute)


DOUBLE_CLICK_INTERVAL = 0.1 # MAGIC HACK
PICKER_TOLERANCE = 10 # MAGIC HACK


class ChemicalAbundancesTab(QtGui.QWidget):
    
    def __init__(self, parent):
        super(ChemicalAbundancesTab, self).__init__(parent)
        self.parent = parent
        self.parent_splitter = QtGui.QSplitter(self)
        self.parent_layout = QtGui.QHBoxLayout(self)
        self.parent_splitter.setContentsMargins(20, 20, 20, 0)
        self.parent_layout.addWidget(self.parent_splitter)
        
        self.spectral_models = []
        
        ################
        # LEFT HAND SIDE
        ################
        lhs_widget = QtGui.QWidget(self)
        lhs_layout = QtGui.QVBoxLayout()
        
        # Abund tree
        self.abundtree = AbundTreeView(self)
        self.abundtree.setModel(AbundTreeModel(self))
        self.abundtree.span_cols()
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, 
            QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.abundtree.sizePolicy().hasHeightForWidth())
        self.abundtree.setSizePolicy(sp)
        lhs_layout.addWidget(self.abundtree)

        # Buttons
        hbox = QtGui.QHBoxLayout()
        self.btn_fit_all = QtGui.QPushButton(self)
        self.btn_fit_all.setText("Fit and Measure all")
        hbox.addWidget(self.btn_fit_all)
        lhs_layout.addLayout(hbox)

        hbox = QtGui.QHBoxLayout()
        self.btn_refresh = QtGui.QPushButton(self)
        self.btn_refresh.setText("Refresh table")
        self.btn_replot  = QtGui.QPushButton(self)
        self.btn_replot.setText("Refresh plots")
        hbox.addWidget(self.btn_refresh)
        hbox.addWidget(self.btn_replot)
        lhs_layout.addLayout(hbox)

        hbox = QtGui.QHBoxLayout()
        self.btn_save_to_session = QtGui.QPushButton(self)
        self.btn_save_to_session.setText("Save to session")
        hbox.addWidget(self.btn_save_to_session)
        lhs_layout.addLayout(hbox)
        
        # TODO Model options
        self.model_options = QtGui.QWidget(self)
        lhs_layout.addWidget(self.model_options)

        lhs_widget.setLayout(lhs_layout)
        self.parent_splitter.addWidget(lhs_widget)

        #############################
        # RIGHT HAND SIDE: MPL WIDGET
        #############################
        rhs_layout = QtGui.QVBoxLayout()
        self.figure = mpl.MPLWidget(None, tight_layout=True, autofocus=True)
        self.figure.setMinimumSize(QtCore.QSize(800, 300))
        
        gs_top = matplotlib.gridspec.GridSpec(3,1,height_ratios=[1,2,1])
        gs_top.update(top=.95,bottom=.05,hspace=0)
        gs_bot = matplotlib.gridspec.GridSpec(3,1,height_ratios=[1,2,1])
        gs_bot.update(top=.95,bottom=.05,hspace=.3)
        
        self._points = None
        self._trend_lines = None
        
        self.ax_residual = self.figure.figure.add_subplot(gs_top[0])
        self.ax_residual.axhline(0, c="#666666")
        self.ax_residual.xaxis.set_major_locator(MaxNLocator(5))
        self.ax_residual.yaxis.set_major_locator(MaxNLocator(2))
        self.ax_residual.set_xticklabels([])
        self.ax_residual.set_ylabel("Residual")
        
        self.ax_spectrum = self.figure.figure.add_subplot(gs_top[1])
        self.ax_spectrum.xaxis.get_major_formatter().set_useOffset(False)
        self.ax_spectrum.xaxis.set_major_locator(MaxNLocator(5))
        self.ax_spectrum.set_xlabel(u"Wavelength (Å)")
        self.ax_spectrum.set_ylabel(r"Normalized flux")
        
        self.ax_line_strength = self.figure.figure.add_subplot(gs_bot[2])
        self.ax_line_strength.xaxis.get_major_formatter().set_useOffset(False)
        self.ax_line_strength.yaxis.set_major_locator(MaxNLocator(5))
        self.ax_line_strength.yaxis.set_major_locator(MaxNLocator(4))
        self.ax_line_strength.set_xlabel(r"$\log({\rm EW}/\lambda)$")
        self.ax_line_strength.set_ylabel("A(X)")
        
        # Some empty figure objects that we will use later.
        self._lines = {
            "selected_point": [
                self.ax_line_strength.scatter([], [],
                    edgecolor="b", facecolor="none", s=150, linewidth=3, zorder=2)
            ],
            "spectrum": None,
            "transitions_center_main": self.ax_spectrum.axvline(
                np.nan, c="#666666", linestyle=":"),
            "transitions_center_residual": self.ax_residual.axvline(
                np.nan, c="#666666", linestyle=":"),
            "model_masks": [],
            "nearby_lines": [],
            "model_fit": self.ax_spectrum.plot([], [], c="r")[0],
            "model_residual": self.ax_residual.plot([], [], c="k")[0],
            "interactive_mask": [
                self.ax_spectrum.axvspan(xmin=np.nan, xmax=np.nan, ymin=np.nan,
                    ymax=np.nan, facecolor="r", edgecolor="none", alpha=0.25,
                    zorder=-5),
                self.ax_residual.axvspan(xmin=np.nan, xmax=np.nan, ymin=np.nan,
                    ymax=np.nan, facecolor="r", edgecolor="none", alpha=0.25,
                    zorder=-5)
            ]
        }
        
        self.parent_splitter.addWidget(self.figure)

        # Connect selection model
        _ = self.abundtree.selectionModel()
        _.selectionChanged.connect(self.selected_model_changed)

        # Connect buttons
        self.btn_fit_all.clicked.connect(self.fit_all)
        self.btn_refresh.clicked.connect(self.refresh_table)
        self.btn_replot.clicked.connect(self.refresh_plots)
        self.btn_save_to_session.clicked.connect(self.save_to_session)

        # Connect matplotlib.
        self.figure.mpl_connect("button_press_event", self.figure_mouse_press)
        self.figure.mpl_connect("button_release_event", self.figure_mouse_release)
        self.figure.figure.canvas.callbacks.connect(
            "pick_event", self.figure_mouse_pick)
        
        self.currently_plotted_species = None
        self.populate_widgets()

    def populate_widgets(self):
        """
        Refresh widgets from session
        """
        self.refresh_table()
        return None

    def refresh_table(self):
        if self.parent.session is None: return None
        self._check_for_spectral_models()
        self.updated_spectral_models() # TODO duplicated
        logger.debug("Resetting tree model")
        self.abundtree.model().reset()
        return None

    def refresh_plots(self):
        self.update_spectrum_figure(False)
        self.update_line_strength_figure(True)
        return None

    def fit_all(self):
        self._check_for_spectral_models()
        self.updated_spectral_models()
        # TODO order by complexity
        logger.debug("Looping through spectral models...")
        for m in self.spectral_models:
            if isinstance(m, SpectralSynthesisModel):
                # TODO 
                logger.debug("Skipping syntheses",m)
                continue
            try:
                res = m.fit()
            except (ValueError, RuntimeError) as e:
                logger.debug("Fitting error",m)
                logger.debug(e)
                continue
            # TODO
            #if not m.is_acceptable:
            #    logger.debug("Skipping",m)
            #    continue
            try:
                ab = m.abundances
            except rt.RTError as e:
                logger.debug("Abundance error",m)
                logger.debug(e)
                continue
        self.abundtree.model().reset()
        self.selected_model_changed()

    def save_to_session(self):
        self.parent.session.metadata["spectral_models"] = self.spectral_models
        logger.debug("ChemicalAbundanceTab: Overwrote session spectral_models!")
        # TODO trigger relevant stuff in other tabs

    def updated_spectral_models(self):
        self.spectral_models = []
        for model in self.parent.session.metadata["spectral_models"]:
            if model.use_for_stellar_composition_inference:
                self.spectral_models.append(model)
        return None
    def _check_for_spectral_models(self):
        for sm in self.parent.session.metadata.get("spectral_models", []):
            if sm.use_for_stellar_composition_inference: break
        else:
            reply = QtGui.QMessageBox.information(self,
                "No spectral models found",
                "No spectral models are currently associated with the "
                "determination of chemical abundances.\n\n"
                "Click 'OK' to load the transitions manager.")
            if reply == QtGui.QMessageBox.Ok:
                # Load line list manager.
                dialog = TransitionsDialog(self.parent.session)
                dialog.exec_()

                # Do we even have any spectral models now?
                for sm in self.parent.session.metadata.get("spectral_models", []):
                    if sm.use_for_stellar_composition_inference: break
                else:
                    return False
                self.updated_spectral_models()
            else:
                return False
        return True


    def figure_mouse_pick(self, event):
        """
        Trigger for when the mouse is used to select an item in the figure.

        :param event:
            The matplotlib event.
        """
        if self.currently_plotted_species is None: return None
        logger.debug("Mouse picked {} from {}".format(event.ind,self.currently_plotted_species))
        model = self.abundtree.model()
        ii = model.all_species == self.currently_plotted_species
        assert np.sum(ii) == 1, "{} {}".format(self.currently_plotted_species, model.all_species)
        summary_index = np.where(ii)[0]
        summary = model.summaries[summary_index]
        item = summary.subnodes[event.ind[0]]
        index = model.createIndex(event.ind[0],0,item)
        self.abundtree.setCurrentIndex(index)
        self.selected_model_changed()
        return None


    def figure_mouse_press(self, event):
        """
        Trigger for when the mouse button is pressed in the figure.

        :param event:
            The matplotlib event.
        """
        logger.debug("Mouse pressed"+str(event))

        if event.inaxes in (self.ax_residual, self.ax_spectrum):
            self.spectrum_axis_mouse_press(event)
        return None


    def figure_mouse_release(self, event):
        """
        Trigger for when the mouse button is released in the figure.

        :param event:
            The matplotlib event.
        """
        logger.debug("Mouse released"+str(event))

        if event.inaxes in (self.ax_residual, self.ax_spectrum):
            self.spectrum_axis_mouse_release(event)
        return None

    def spectrum_axis_mouse_press(self, event):
        """
        The mouse button was pressed in the spectrum axis.

        :param event:
            The matplotlib event.
        """
        logger.debug("Spectrum pressed"+str(event))

        if event.dblclick:

            # Double click.
            spectral_model, index = self._get_selected_model(True)
            if spectral_model is None:
                return None #TODO is this right?
            for i, (s, e) in enumerate(spectral_model.metadata["mask"][::-1]):
                if e >= event.xdata >= s:
                    # Remove a mask
                    mask = spectral_model.metadata["mask"]
                    index = len(mask) - 1 - i
                    del mask[index]

                    # Re-fit the current spectral_model.
                    spectral_model.fit()

                    # Update the table view for this row.
                    tree_model = self.abundtree.model()
                    # It will be pointing to an item
                    item = index.internalPointer() 
                    tree_model.dataChanged.emit(
                        tree_model.createIndex(0, 0, item),
                        tree_model.createIndex(0, item.columnCount(), item))

                    # Update the view of the current model.
                    self.update_spectrum_figure(True)
                    break

            else:
                # No match with a masked region. 
                # TODO: Add a point that will be used for the continuum?
                # For the moment just refit the model.
                spectral_model.fit()

                # Update the table view for this row.
                tree_model = self.abundtree.model()
                # It will be pointing to an item
                item = index.internalPointer() 
                tree_model.dataChanged.emit(
                    tree_model.createIndex(0, 0, item),
                    tree_model.createIndex(0, item.columnCount(), item))

                # Update the view of the current model.
                self.update_spectrum_figure(True)
                return None

        else:
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

            # Set the signal and the time.
            self._interactive_mask_region_signal = (
                time.time(),
                self.figure.mpl_connect(
                    "motion_notify_event", self.update_mask_region)
            )

        return None


    def update_mask_region(self, event):
        """
        Update the visible selected masked region for the selected spectral
        model. This function is linked to a callback for when the mouse position
        moves.

        :para event:
            The matplotlib motion event to show the current mouse position.
        """

        if event.xdata is None: return

        signal_time, signal_cid = self._interactive_mask_region_signal
        if time.time() - signal_time > DOUBLE_CLICK_INTERVAL:

            data = self._lines["interactive_mask"][0].get_xy()

            # Update xmax.
            data[2:4, 0] = event.xdata
            for patch in self._lines["interactive_mask"]:
                patch.set_xy(data)

            self.figure.draw()

        return None



    def spectrum_axis_mouse_release(self, event):
        """
        Mouse button was released from the spectrum axis.

        :param event:
            The matplotlib event.
        """
        logger.debug("Spectrum released"+str(event))

        try:
            signal_time, signal_cid = self._interactive_mask_region_signal

        except AttributeError:
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
            
            # Get current spectral model.
            spectral_model, index = self._get_selected_model(True)
            if spectral_model is None: 
                raise RuntimeError("""Must have a spectral model selected while making mask!
                                   Must have mouseover bug?""")

            # Add mask metadata.
            spectral_model.metadata["mask"].append([xy[0,0], xy[2, 0]])

            # Re-fit the spectral model.
            spectral_model.fit()

            # Update the table view for this row.
            tree_model = self.abundtree.model()
            # It will be pointing to an item
            item = index.internalPointer() 
            tree_model.dataChanged.emit(
                tree_model.createIndex(0, 0, item),
                tree_model.createIndex(0, item.columnCount(), item))

            # Update the view of the spectral model.
            self.update_spectrum_figure()

        xy[:, 0] = np.nan
        for patch in self._lines["interactive_mask"]:
            patch.set_xy(xy)

        self.figure.mpl_disconnect(signal_cid)
        self.figure.draw()
        del self._interactive_mask_region_signal
        return None

    def _get_selected_model(self, full_output=False):
        index = self.abundtree.selectionModel().currentIndex()
        item = index.internalPointer()
        if isinstance(item,AbundTreeMeasurementItem):
            model = self.spectral_models[item.sm_ix]
            return (model, index) if full_output else model
        else:
            return (None, None) if full_output else None

    def selected_model_changed(self):
        try:
            selected_model = self._get_selected_model()
        except IndexError:
            for collection in self._lines["selected_point"]:
                collection.set_offsets(np.array([np.nan, np.nan]).T)
            self.figure.draw()
            return None
        if selected_model is None: 
            for collection in self._lines["selected_point"]:
                collection.set_offsets(np.array([np.nan, np.nan]).T)
            self.figure.draw()
            return None
        
        logger.debug("Changing selected model: "+str(selected_model))
        # TODO currently assumes all models are ProfileFittingModel
        assert isinstance(selected_model, ProfileFittingModel)
        try:
            metadata = selected_model.metadata["fitted_result"][-1]
            abundances = metadata["abundances"]
            equivalent_width = metadata["equivalent_width"][0]
        except (IndexError, KeyError):
            abundances = [np.nan]
        
        abundance = abundances[0]
        if not np.isfinite(abundance):
            excitation_potential, rew = (np.nan, np.nan)
        else:
            transitions = selected_model.transitions
            assert len(transitions) == 1
            excitation_potential = transitions["expot"][0]
            rew = np.log10(equivalent_width/transitions["wavelength"][0])
            
        point_strength = self._lines["selected_point"][0]
        point_strength.set_offsets(np.array([rew, abundance]).T)
        
        self.update_spectrum_figure(False)
        self.update_line_strength_figure(True)
        
        return None

    def update_spectrum_figure(self, refresh=False):
        """
        TODO refactor
        Currently copied straight from stellar_parameters.py with minor changes
        """
        if self._lines["spectrum"] is None \
        and hasattr(self.parent, "session") \
        and hasattr(self.parent.session, "normalized_spectrum"):
            # Draw the spectrum.
            spectrum = self.parent.session.normalized_spectrum
            self._lines["spectrum"] = self.ax_spectrum.plot(spectrum.dispersion,
                spectrum.flux, c="k", drawstyle="steps-mid")

            sigma = 1.0/np.sqrt(spectrum.ivar)
            style_utils.fill_between_steps(self.ax_spectrum, spectrum.dispersion,
                spectrum.flux - sigma, spectrum.flux + sigma, 
                facecolor="#cccccc", edgecolor="#cccccc", alpha=1)

            style_utils.fill_between_steps(self.ax_residual, spectrum.dispersion,
                -sigma, +sigma, facecolor="#CCCCCC", edgecolor="none", alpha=1)

            self.ax_spectrum.set_xlim(
                spectrum.dispersion[0], spectrum.dispersion[-1])
            self.ax_residual.set_xlim(self.ax_spectrum.get_xlim())
            self.ax_spectrum.set_ylim(0, 1.2)
            self.ax_spectrum.set_yticks([0, 0.5, 1])
            self.ax_residual.set_ylim(-0.05, 0.05)

            if refresh: self.figure.draw()
        
        selected_model = self._get_selected_model()
        transitions = selected_model.transitions
        window = selected_model.metadata["window"]
        limits = [
            transitions["wavelength"][0] - window,
            transitions["wavelength"][-1] + window,
        ]

        # Zoom to region.
        self.ax_spectrum.set_xlim(limits)
        self.ax_residual.set_xlim(limits)
            
        # If this is a profile fitting line, show where the centroid is.
        x = transitions["wavelength"][0] \
            if isinstance(selected_model, ProfileFittingModel) else np.nan
        self._lines["transitions_center_main"].set_data([x, x], [0, 1.2])
        self._lines["transitions_center_residual"].set_data([x, x], [0, 1.2])
        # Model masks specified by the user.
        # (These should be shown regardless of whether there is a fit or not.)
        for i, (start, end) in enumerate(selected_model.metadata["mask"]):
            try:
                patches = self._lines["model_masks"][i]

            except IndexError:
                self._lines["model_masks"].append([
                    self.ax_spectrum.axvspan(np.nan, np.nan,
                        facecolor="r", edgecolor="none", alpha=0.25),
                    self.ax_residual.axvspan(np.nan, np.nan,
                        facecolor="r", edgecolor="none", alpha=0.25)
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
                patch.set_visible(True)

        # Hide unnecessary ones.
        N = len(selected_model.metadata["mask"])
        for unused_patches in self._lines["model_masks"][N:]:
            for unused_patch in unused_patches:
                unused_patch.set_visible(False)

        # Hide previous model_errs
        try:
            self._lines["model_yerr"].set_visible(False)
            del self._lines["model_yerr"]
            # TODO: This is wrong. It doesn't actually delete them so if
            #       you ran this forever then you would get a real bad 
            #       memory leak in Python. But for now, re-calculating
            #       the PolyCollection is in the too hard basket.

        except KeyError:
            None

        # Things to show if there is a fitted result.
        try:
            (named_p_opt, cov, meta) = selected_model.metadata["fitted_result"]

            # Test for some requirements.
            _ = (meta["model_x"], meta["model_y"], meta["residual"])

        except KeyError:
            meta = {}
            self._lines["model_fit"].set_data([], [])
            self._lines["model_residual"].set_data([], [])

        else:
            assert len(meta["model_x"]) == len(meta["model_y"])
            assert len(meta["model_x"]) == len(meta["residual"])
            assert len(meta["model_x"]) == len(meta["model_yerr"])

            self._lines["model_fit"].set_data(meta["model_x"], meta["model_y"])
            self._lines["model_residual"].set_data(meta["model_x"], 
                meta["residual"])

            # Model yerr.
            if np.any(np.isfinite(meta["model_yerr"])):
                self._lines["model_yerr"] = self.ax_spectrum.fill_between(
                    meta["model_x"],
                    meta["model_y"] + meta["model_yerr"],
                    meta["model_y"] - meta["model_yerr"],
                    facecolor="r", edgecolor="none", alpha=0.5)

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

        if refresh: self.figure.draw()

        return None

    def update_line_strength_figure(self, refresh=False):
        selected_model, index = self._get_selected_model(True)
        if selected_model is None:
            # TODO clear plot?
            return None

        item = index.internalPointer() # AbundTreeMeasurementItem
        summary = item.parent

        # If the species is already plotted, don't replot
        assert isinstance(selected_model, ProfileFittingModel)
        ## TODO this doesn't update the plot when selecting/deselecting
        #if selected_model.transitions["species"][0] == self.currently_plotted_species:
        #    if refresh: self.figure.draw()
        #    return None
        
        rew_list = []
        abund_list = []
        for node in summary.subnodes:
            m = self.spectral_models[node.sm_ix]
            if not m.is_acceptable: 
                # Append these to keep indices straight between [de]selected points and the plot
                rew_list.append(np.nan)
                abund_list.append(np.nan)
                continue
            # TODO SpectralSynthesisModel
            assert isinstance(m, ProfileFittingModel), m
            try:
                metadata = m.metadata["fitted_result"][-1]
                abundances = metadata["abundances"]
                equivalent_width = metadata["equivalent_width"][0]
            except (IndexError, KeyError):
                abundances = [np.nan]
            abundance = abundances[0]
            if not np.isfinite(abundance):
                rew = np.nan
            else:
                rew = np.log10(equivalent_width/m.transitions["wavelength"][0])
            rew_list.append(rew)
            abund_list.append(abundance)
        if self._points is None:
            self._points = [self.ax_line_strength.scatter([], [], s=30, \
                                 facecolor="k", edgecolor="k", picker=PICKER_TOLERANCE, \
                                 alpha=0.5)]
        collections = self._points
        collections[0].set_offsets(np.array([rew_list,abund_list]).T)
        style_utils.relim_axes(self.ax_line_strength)
        
        # TODO this will work for now but is a hack
        self.currently_plotted_species = m.transitions["species"][0]
        
        # TODO trend lines
        
        if refresh: self.figure.draw()
        return None
