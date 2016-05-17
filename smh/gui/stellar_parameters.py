#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The stellar parameters tab in Spectroscopy Made Hard """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["StellarParametersTab"]

import logging
import numpy as np
import sys
from PySide import QtCore, QtGui
from matplotlib import gridspec

import mpl
from smh.photospheres import available as available_photospheres


logger = logging.getLogger(__name__)


class StellarParametersTab(QtGui.QWidget):

    def __init__(self, parent=None):
        super(StellarParametersTab, self).__init__(parent)
        self.parent = parent

        # Establish the GUI for this tab.
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding,
            QtGui.QSizePolicy.MinimumExpanding
        )
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)

        # Create a top-level horizontal layout to contain a MPL figure and
        # a vertical layout of settings.
        tab_layout = QtGui.QHBoxLayout(self)
        tab_layout.setContentsMargins(10, 10, 10, 10)
        
        input_parameters = QtGui.QWidget()
        input_parameters_layout = QtGui.QVBoxLayout(input_parameters)
        input_parameters.setFixedWidth(300)
        
        # Top-level button to measure transitions.



        # Start the grid layout for the stellar parameters tab.
        input_parameters_grid = QtGui.QGridLayout()

        # Photospheres.
        self.photospheres_label = QtGui.QLabel(self)
        self.photospheres_label.setText("Photospheres")
        input_parameters_grid.addWidget(self.photospheres_label, 0, 0, 1, 1)
        
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.photospheres = QtGui.QComboBox(self)
        hbox.addWidget(self.photospheres)
        input_parameters_grid.addLayout(hbox, 0, 1, 1, 1)

        for description, kind, basename in available_photospheres:
            self.photospheres.addItem(description)

        # Effective temperature.
        self.effective_temperature_label = QtGui.QLabel(self)
        self.effective_temperature_label.setText("Effective temperature (K)")
        input_parameters_grid.addWidget(
            self.effective_temperature_label, 1, 0, 1, 1)
        
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))

        self.effective_temperature = QtGui.QLineEdit(self)
        self.effective_temperature.setMaximumSize(QtCore.QSize(40, 16777215))
        self.effective_temperature.setAlignment(QtCore.Qt.AlignCenter)
        self.effective_temperature.setValidator(
            QtGui.QDoubleValidator(0, 1000, 2, self.effective_temperature))
        hbox.addWidget(self.effective_temperature)
        input_parameters_grid.addLayout(hbox, 1, 1, 1, 1)



        input_parameters_layout.addLayout(input_parameters_grid)

        # Add a spacer.
        input_parameters_layout.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding))



        # Add a 'Measure abundances' button.
        self.measure_abundances = QtGui.QPushButton(self)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(
            self.measure_abundances.sizePolicy().hasHeightForWidth())
        self.measure_abundances.setSizePolicy(sp)
        self.measure_abundances.setMinimumSize(QtCore.QSize(300, 0))
        self.measure_abundances.setMaximumSize(QtCore.QSize(300, 16777215))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.measure_abundances.setFont(font)
        self.measure_abundances.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.measure_abundances.setDefault(True)
        self.measure_abundances.setObjectName("measure_abundances")
        self.measure_abundances.setText("Measure abundances")
        if sys.platform == "darwin":
            self.measure_abundances.setStyleSheet('QPushButton {color: white}')

        input_parameters_layout.addWidget(self.measure_abundances)


        tab_layout.addWidget(input_parameters)


        # Create a matplotlib widget.
        blank_widget = QtGui.QWidget(self)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(blank_widget.sizePolicy().hasHeightForWidth())
        blank_widget.setSizePolicy(sp)
        blank_widget.setObjectName("sp_plot")


        self.figure = mpl.MPLWidget(blank_widget, tight_layout=True,
            autofocus=True)

        layout = QtGui.QVBoxLayout(blank_widget)
        layout.addWidget(self.figure, 1)
        tab_layout.addWidget(blank_widget)

        # Set up the plot.
        gs = gridspec.GridSpec(2, 1)
        self.ax_first = self.figure.figure.add_subplot(gs[0])

        # Line for the data.
        self.ax_first.plot([0, 0.5], [0.5, 1], c='k', zorder=3)#, drawstyle='steps-mid')
        # Line for the continuum.
        self.ax_first.plot([], [], linestyle="--", linewidth=2, c='r', zorder=4)

        # Line for the neighbouring order(s) (joined by a NaN).
        self.ax_first.plot([], [], c='#666666', zorder=1, drawstyle='steps-mid')
        # Line for the neighbouring order(s) continuum (joined by a NaN)
        self.ax_first.plot([], [], c='b', zorder=2)

        # Additional point markers.
        self.ax_first.scatter([], [], facecolor="k", zorder=5, picker=5)

        # Regions

        self.ax_first.set_xticklabels([])
        self.ax_first.set_yticklabels([])
        self.ax_first.set_ylabel("Flux")

        self.ax_second = self.figure.figure.add_subplot(gs[1])
        self.ax_second.axhline(1, linestyle=":", c="#666666", zorder=1)
        self.ax_second.plot([], [], c='k', zorder=2)

        # TODO: Make (0, 1.2) a default view setting.
        self.ax_second.set_ylim(0, 1.2)
        self.ax_second.set_yticks([0, 0.5, 1.0])
        self.ax_second.set_xlabel(u"Wavelength (Å)")

        #self.figure.figure.tight_layout(w_pad=0, h_pad=0, pad=0.4)
        self.figure.draw()

        """
        # Create signals.
        self.measure_abundances.clicked.connect(self.normalize_and_stitch)

        # Note that key_press_event is linked to figure.canvas, while the
        # mouse events are linked to figure.
        # I don't know why, but that's how it works.
        self.figure.canvas.mpl_connect(
            "key_press_event", self.figure_key_press)
        self.figure.mpl_connect(
            "button_press_event", self.figure_mouse_press)
        self.figure.mpl_connect(
            "button_release_event", self.figure_mouse_release)
        
        self.function.currentIndexChanged.connect(
            self.update_normalization_function)
        self.order.currentIndexChanged.connect(
            self.update_normalization_order)
        self.norm_max_iter.currentIndexChanged.connect(
            self.update_normalization_max_iterations)
        self.low_sigma_clip.textChanged.connect(
            self.update_low_sigma_clip)
        self.high_sigma_clip.textChanged.connect(
            self.update_high_sigma_clip)
        self.knot_spacing.textChanged.connect(self.update_knot_spacing)

        self.low_sigma_clip.textChanged.connect(self.check_state)
        self.high_sigma_clip.textChanged.connect(self.check_state)
        self.knot_spacing.textChanged.connect(self.check_state)
        """



