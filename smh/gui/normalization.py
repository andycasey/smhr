#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The normalization tab in Spectroscopy Made Hard """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import sys
from PySide import QtCore, QtGui

from matplotlib import gridspec

import mpl

__all__ = ["NormalizationTab"]


c = 299792458e-3 # km/s

class NormalizationTab(QtGui.QWidget):

    def __init__(self, parent=None):
        super(NormalizationTab, self).__init__(parent)
        self.parent = parent

        # Establish the GUI for this tab.
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)

        # Create a top-level horizontal layout to contain a MPL figure and
        # a vertical layout of settings..
        tab_layout = QtGui.QHBoxLayout(self)
        tab_layout.setContentsMargins(10, 10, 10, 10)
        
        settings_widget = QtGui.QWidget()
        settings_layout = QtGui.QVBoxLayout(settings_widget)
        settings_widget.setFixedWidth(300)
        
        # Start the grid layout for the normalization tab.
        settings_grid_layout = QtGui.QGridLayout()

        # Normalization function.
        label = QtGui.QLabel(self)
        label.setText("Function")
        settings_grid_layout.addWidget(label, 0, 0, 1, 1)
        
        # Put the normalization function combo box in a horizontal layout with 
        # a spacer.
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_function = QtGui.QComboBox(self)
        self.norm_function.setObjectName("rv_norm_function")
        hbox.addWidget(self.norm_function)
        settings_grid_layout.addLayout(hbox, 0, 1, 1, 1)

        norm_functions = ("polynomial", "spline")
        for each in norm_functions:
            self.norm_function.addItem(each.title())

        # Normalization function order.
        label = QtGui.QLabel(self)
        label.setText("Order")
        settings_grid_layout.addWidget(label, 1, 0, 1, 1)
        
        # Put the normalization order combo box in a horizontal layout with a
        # spacer
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_order = QtGui.QComboBox(self)
        self.norm_order.setMaximumSize(QtCore.QSize(50, 16777215))
        self.norm_order.setObjectName("rv_norm_order")
        hbox.addWidget(self.norm_order)
        settings_grid_layout.addLayout(hbox, 1, 1, 1, 1)

        norm_orders = range(1, 10)
        for order in norm_orders:
            self.norm_order.addItem("{0:.0f}".format(order))

        # Maximum number of iterations.
        label = QtGui.QLabel(self)
        label.setText("Maximum iterations")
        settings_grid_layout.addWidget(label, 2, 0, 1, 1)

        # Put the maxium number of iterations in a horizontal layout with a 
        # spacer.
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_max_iter = QtGui.QComboBox(self)
        self.norm_max_iter.setMaximumSize(QtCore.QSize(50, 16777215))
        self.norm_max_iter.setObjectName("rv_norm_max_iter")
        hbox.addWidget(self.norm_max_iter)
        settings_grid_layout.addLayout(hbox, 2, 1, 1, 1)

        norm_max_iters = range(1, 10)
        for iteration in norm_max_iters:
            self.norm_max_iter.addItem("{0:.0f}".format(iteration))


        # Low sigma clipping.
        label = QtGui.QLabel(self)
        label.setText("Low sigma clip")
        settings_grid_layout.addWidget(label, 3, 0, 1, 1)

        # Put the low sigma line edit box in a horizontal layout with a spacer.
        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(-1, -1, 5, -1)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_low_sigma = QtGui.QLineEdit(self)
        self.norm_low_sigma.setMaximumSize(QtCore.QSize(40, 16777215))
        self.norm_low_sigma.setAlignment(QtCore.Qt.AlignCenter)
        self.norm_low_sigma.setObjectName("rv_norm_low_sigma")
        hbox.addWidget(self.norm_low_sigma)
        settings_grid_layout.addLayout(hbox, 3, 1, 1, 1)


        # High sigma clipping.
        label = QtGui.QLabel(self)
        label.setText("High sigma clip")
        settings_grid_layout.addWidget(label, 4, 0, 1, 1)

        # Put the high sigma line edit box in a horizontal layout with a spacer.
        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(-1, -1, 5, -1)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_high_sigma = QtGui.QLineEdit(self)
        self.norm_high_sigma.setMaximumSize(QtCore.QSize(40, 16777215))
        self.norm_high_sigma.setAlignment(QtCore.Qt.AlignCenter)
        self.norm_high_sigma.setObjectName("rv_norm_high_sigma")
        hbox.addWidget(self.norm_high_sigma)
        settings_grid_layout.addLayout(hbox, 4, 1, 1, 1)
        

        # Knot spacing.
        label = QtGui.QLabel(self)
        settings_grid_layout.addWidget(label, 5, 0, 1, 1)
        label.setText(u"Knot spacing (Å)")

        # Put the knot spacing lint edit box in a horizontal layout with a spacer
        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(-1, -1, 5, -1)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_knot_spacing = QtGui.QLineEdit(self)
        self.norm_knot_spacing.setMaximumSize(QtCore.QSize(40, 16777215))
        self.norm_knot_spacing.setAlignment(QtCore.Qt.AlignCenter)
        self.norm_knot_spacing.setObjectName("rv_norm_knot_spacing")
        hbox.addWidget(self.norm_knot_spacing)
        settings_grid_layout.addLayout(hbox, 5, 1, 1, 1)


        # End of the grid in the normalization tab.
        settings_layout.addLayout(settings_grid_layout)

        # Add a label.
        label = QtGui.QLabel(self)
        label.setText("Continuum mask:")
        settings_layout.addWidget(label)

        # Add options for continuum mask.
        hbox = QtGui.QHBoxLayout()
        self.continuum_mask = QtGui.QComboBox(self)
        self.continuum_mask.setObjectName("contiuum_mask")
        hbox.addWidget(self.continuum_mask)
        hbox.addItem(
            QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
                QtGui.QSizePolicy.Minimum))
        self.edit_continuum_masks = QtGui.QPushButton(self)
        self.edit_continuum_masks.setObjectName("edit_continuum_masks")
        self.edit_continuum_masks.setText("Edit masks..")
        hbox.addWidget(self.edit_continuum_masks)

        settings_layout.addLayout(hbox)


        # Add a spacer.
        settings_layout.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding))

        # Add a 'normalize and stitch button'
        norm_stitch_btn = QtGui.QPushButton(self)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(norm_stitch_btn.sizePolicy().hasHeightForWidth())
        norm_stitch_btn.setSizePolicy(sp)
        norm_stitch_btn.setMinimumSize(QtCore.QSize(300, 0))
        norm_stitch_btn.setMaximumSize(QtCore.QSize(300, 16777215))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        norm_stitch_btn.setFont(font)
        norm_stitch_btn.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        norm_stitch_btn.setDefault(True)
        norm_stitch_btn.setObjectName("norm_stitch_btn")
        norm_stitch_btn.setText("Normalize and stitch orders")
        if sys.platform == "darwin":
            norm_stitch_btn.setStyleSheet('QPushButton {color: white}')

        settings_layout.addWidget(norm_stitch_btn)


        tab_layout.addWidget(settings_widget)

        # Create a matplotlib widget.
        blank_widget = QtGui.QWidget(self)
        blank_widget.setStatusTip(
            "Use left/right arrows to move between orders, "
            "double-click for all keyboard shortcuts")
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(blank_widget.sizePolicy().hasHeightForWidth())
        blank_widget.setSizePolicy(sp)
        blank_widget.setObjectName("norm_plot")


        self.norm_plot = mpl.MPLWidget(blank_widget, tight_layout=True,
            autofocus=True)

        layout = QtGui.QVBoxLayout(blank_widget)
        layout.addWidget(self.norm_plot, 1)
        tab_layout.addWidget(blank_widget)

        # Set up the plot.
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        self.ax_order = self.norm_plot.figure.add_subplot(gs[0])
        # Line for the data.
        self.ax_order.plot([], [], c='k', zorder=3)#, drawstyle='steps-mid')
        # Line for the continuum.
        self.ax_order.plot([], [], c='r', zorder=4)

        # Line for the neighbouring order(s) (joined by a NaN).
        self.ax_order.plot([], [], c='#666666', zorder=1, drawstyle='steps-mid')
        # Line for the neighbouring order(s) continuum (joined by a NaN)
        self.ax_order.plot([], [], c='b', zorder=2)

        # Additional point markers.
        # TODO
        # Regions

        self.ax_order.set_xticklabels([])
        self.ax_order.set_yticklabels([])
        self.ax_order.set_ylabel("Flux")

        self.ax_order_norm = self.norm_plot.figure.add_subplot(gs[1])
        self.ax_order_norm.set_ylim(0, 1.2)
        self.ax_order_norm.set_yticks([0, 0.5, 1.0])
        self.ax_order_norm.set_xlabel(u"Wavelength (Å)")

        #self.norm_plot.figure.tight_layout(w_pad=0, h_pad=0, pad=0.4)
        self.norm_plot.draw()

        # Create signals.
        self.norm_plot.canvas.mpl_connect(
            "key_press_event", self.figure_key_press)

        return None


    def figure_key_press(self, event):
        """
        Key press event in the normalization figure.
        """

        print(event, event.key)

        if event.key in ("left", "right"):
            offset = 1 if event.key == "right" else -1

            self.update_order_index(self.current_order_index + offset)
            self.draw_order()
            self.update_continuum()
            self.draw_continuum(True)

            print("DID IT", self.current_order_index)
            return None

        # TODO: deal with discarded order indices, etc.


    def _populate_widgets(self):
        """
        Populate the widgets with default settings, and a default figure view.
        """

        if self.parent.session is None:
            # No point populating the widgets with the default values from the
            # SMH file because these will be updated when a session is loaded.
            return

        self._cache = {
            "input": {}
        }

        # A session exists. Load up everything.
        self.update_order_index(0)
        self.update_continuum()
        self.draw_order()
        self.draw_continuum(True)


    def update_order_index(self, index=None):
        """
        Update the currently selected order.
        """
        if index is None:
            index = self.current_order_index

        self.current_order_index = index
        self.current_order \
            = self.parent.session.input_spectra[self.current_order_index].copy()

        # Apply any RV correction.
        try:
            rv_applied = self.parent.session.metadata["rv"]["rv_applied"]
        except (AttributeError, KeyError):
            rv_applied = 0

        self.current_order._dispersion *= (1 - rv_applied/c)

        return None



    def draw_order(self, refresh=False):
        """
        Draw the current order.
        """

        x, y = self.current_order.dispersion, self.current_order.flux
        self.ax_order.lines[0].set_data([x, y])
        self.ax_order.set_xlim(x[0], x[-1])
        self.ax_order.set_ylim(np.nanmin(y), np.nanmax(y))

        self.ax_order.set_title("Order {0} of {1}".format(
            1 + self.current_order_index, 
            len(self.parent.session.input_spectra)))

        if refresh:
            self.norm_plot.draw()
        return None


    def update_continuum(self):
        """
        Update continuum for the current order.
        """

        kwds = self._cache["input"].copy()
        kwds["full_output"] = True

        normalized_spectrum, continuum, left, right \
            = self.current_order.fit_continuum(**kwds)

        index = self.current_order_index
        self.parent.session.metadata["normalization"]["continuum"][index] \
            = continuum

        return None


    def draw_continuum(self, refresh=False):
        """
        Draw the continuum for the current order.
        """

        meta = self.parent.session.metadata["normalization"]
        continuum = meta["continuum"][self.current_order_index]

        self.ax_order.lines[1].set_data([
            self.current_order.dispersion, continuum])

        if refresh:
            self.norm_plot.draw()
        return None

