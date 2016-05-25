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
from smh.spectral_models import (ProfileFittingModel, SpectralSynthesisModel)
from abund_tree import AbundTreeView, AbundTreeModel
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
        self.btn_refresh = QtGui.QPushButton(self)
        self.btn_refresh.setText("Refresh table")
        hbox.addWidget(self.btn_refresh)
        lhs_layout.addLayout(hbox)
        
        # Model options
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
        
        self._points = {}
        self._trend_lines = {}
        
        self.ax_residual = self.figure.figure.add_subplot(gs_top[0])
        self.ax_residual.axhline(0, c="#666666")
        self.ax_residual.xaxis.set_major_locator(MaxNLocator(5))
        self.ax_residual.yaxis.set_major_locator(MaxNLocator(2))
        self.ax_residual.set_xticklabels([])
        
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
        self.ax_line_strength.set_ylabel("[X/H]") # TODO: X/M
        
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

        # Connect buttons
        self.btn_refresh.clicked.connect(self.refresh_table)

        # Connect matplotlib.
        self.figure.mpl_connect("button_press_event", self.figure_mouse_press)
        self.figure.mpl_connect("button_release_event", self.figure_mouse_release)
        self.figure.figure.canvas.callbacks.connect(
            "pick_event", self.figure_mouse_pick)
        
    def populate_widgets(self):
        """
        Refresh widgets from session
        """
        return None

    def refresh_table(self):
        self._check_for_spectral_models()
        self.updated_spectral_models() # TODO duplicated
        print("Resetting model")
        self.abundtree.model().reset()

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

        # Select the row(s). That will trigger the rest.
        #for index in self._spectral_model_indices[event.ind]:
        #    self.table_view.selectRow(index)
        return None


    def figure_mouse_press(self, event):
        """
        Trigger for when the mouse button is pressed in the figure.

        :param event:
            The matplotlib event.
        """

        if event.inaxes in (self.ax_residual, self.ax_spectrum):
            self.spectrum_axis_mouse_press(event)
        return None


    def figure_mouse_release(self, event):
        """
        Trigger for when the mouse button is released in the figure.

        :param event:
            The matplotlib event.
        """

        if event.inaxes in (self.ax_residual, self.ax_spectrum):
            self.spectrum_axis_mouse_release(event)
        return None

    def spectrum_axis_mouse_press(self, event):
        """
        The mouse button was pressed in the spectrum axis.

        :param event:
            The matplotlib event.
        """

        if event.dblclick:

            # Double click.
            spectral_model, index = self._get_selected_model(True)
            for i, (s, e) in enumerate(spectral_model.metadata["mask"][::-1]):
                if e >= event.xdata >= s:

                    mask = spectral_model.metadata["mask"]
                    index = len(mask) - 1 - i
                    del mask[index]

                    # Re-fit the current spectral_model.
                    spectral_model.fit()

                    # Update the table view for this row.
                    table_model = self.table_view.model()
                    table_model.dataChanged.emit(
                        table_model.createIndex(index, 0),
                        table_model.createIndex(
                            index, table_model.columnCount(0)))

                    # Update the view of the current model.
                    self.update_spectrum_figure()
                    break

            else:
                # No match with a masked region. 

                # TODO: Add a point that will be used for the continuum?

                # For the moment just refit the model.
                spectral_model.fit()

                # Update the table view for this row.
                table_model = self.table_view.model()
                table_model.dataChanged.emit(
                    table_model.createIndex(index, 0),
                    table_model.createIndex(
                        index, table_model.columnCount(0)))

                # Update the view of the current model.
                self.update_spectrum_figure()
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
                time(),
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
        if time() - signal_time > DOUBLE_CLICK_INTERVAL:

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
        if  time() - signal_time > DOUBLE_CLICK_INTERVAL \
        and np.abs(xy[0,0] - xdata) > 0:
            
            # Get current spectral model.
            spectral_model, index = self._get_selected_model(True)

            # Add mask metadata.
            spectral_model.metadata["mask"].append([xy[0,0], xy[2, 0]])

            # Re-fit the spectral model.
            spectral_model.fit()

            # Update the table view for this row.
            table_model = self.table_view.model()
            table_model.dataChanged.emit(
                table_model.createIndex(index, 0),
                table_model.createIndex(
                    index, table_model.columnCount(0)))

            # Update the view of the spectral model.
            self.update_spectrum_figure()

        xy[:, 0] = np.nan
        for patch in self._lines["interactive_mask"]:
            patch.set_xy(xy)

        self.figure.mpl_disconnect(signal_cid)
        self.figure.draw()
        del self._interactive_mask_region_signal
        return None

