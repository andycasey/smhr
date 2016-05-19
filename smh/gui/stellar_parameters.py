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
from linelist_manager import TransitionsDialog

logger = logging.getLogger(__name__)

if sys.platform == "darwin":
        
    # See http://successfulsoftware.net/2013/10/23/fixing-qt-4-for-mac-os-x-10-9-mavericks/
    substitutes = [
        (".Lucida Grande UI", "Lucida Grande"),
        (".Helvetica Neue DeskInterface", "Helvetica Neue")
    ]
    for substitute in substitutes:
        QtGui.QFont.insertSubstitution(*substitute)


class StellarParametersTab(QtGui.QWidget):

    def __init__(self, parent):
        """
        Create a tab for the determination of stellar parameters by excitation
        and ionization equalibrium.

        :param parent:
            The parent widget.
        """

        super(StellarParametersTab, self).__init__(parent)

        panel_size = 350
        self.parent = parent

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
        input_parameters.setFixedWidth(panel_size)
        
        # Top-level button to measure transitions.
        self.btn_measure_transitions = QtGui.QPushButton(self)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(
            self.btn_measure_transitions.sizePolicy().hasHeightForWidth())
        self.btn_measure_transitions.setSizePolicy(sp)
        self.btn_measure_transitions.setMinimumSize(
            QtCore.QSize(panel_size, 0))
        self.btn_measure_transitions.setMaximumSize(
            QtCore.QSize(panel_size, 16777215))

        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.btn_measure_transitions.setFont(font)
        self.btn_measure_transitions.setCursor(
            QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.btn_measure_transitions.setDefault(True)
        self.btn_measure_transitions.setObjectName("btn_measure_transitions")
        self.btn_measure_transitions.setText("Measure transitions..")
        if sys.platform == "darwin":
            self.btn_measure_transitions.setStyleSheet(
                'QPushButton {color: white}')

        input_parameters_layout.addWidget(self.btn_measure_transitions)


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
            QtGui.QDoubleValidator(3000, 8000, 0, self.effective_temperature))
        hbox.addWidget(self.effective_temperature)
        input_parameters_grid.addLayout(hbox, 1, 1, 1, 1)

        # Surface gravity.
        self.surface_gravity_label = QtGui.QLabel(self)
        self.surface_gravity_label.setText("Surface gravity")
        input_parameters_grid.addWidget(
            self.surface_gravity_label, 2, 0, 1, 1)
        
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))

        self.surface_gravity = QtGui.QLineEdit(self)
        self.surface_gravity.setMaximumSize(QtCore.QSize(40, 16777215))
        self.surface_gravity.setAlignment(QtCore.Qt.AlignCenter)
        self.surface_gravity.setValidator(
            QtGui.QDoubleValidator(-1, 6, 2, self.surface_gravity))
        hbox.addWidget(self.surface_gravity)
        input_parameters_grid.addLayout(hbox, 2, 1, 1, 1)

        # Metallicity.
        self.metallicity_label = QtGui.QLabel(self)
        self.metallicity_label.setText("Metallicity")
        input_parameters_grid.addWidget(
            self.metallicity_label, 3, 0, 1, 1)
        
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))

        self.metallicity = QtGui.QLineEdit(self)
        self.metallicity.setMaximumSize(QtCore.QSize(40, 16777215))
        self.metallicity.setAlignment(QtCore.Qt.AlignCenter)
        self.metallicity.setValidator(
            QtGui.QDoubleValidator(-7, 1, 2, self.metallicity))
        hbox.addWidget(self.metallicity)
        input_parameters_grid.addLayout(hbox, 3, 1, 1, 1)

        # Depending on the photospheres: alpha-enhancement.


        # Depending on the radiative transfer code used: microturbulence.
        self.microturbulence_label = QtGui.QLabel(self)
        self.microturbulence_label.setText("Microturbulence (km/s)")
        input_parameters_grid.addWidget(
            self.microturbulence_label, 4, 0, 1, 1)
        
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))

        self.microturbulence = QtGui.QLineEdit(self)
        self.microturbulence.setMaximumSize(QtCore.QSize(40, 16777215))
        self.microturbulence.setAlignment(QtCore.Qt.AlignCenter)
        self.microturbulence.setValidator(
            QtGui.QDoubleValidator(0, 5, 2, self.microturbulence))
        hbox.addWidget(self.microturbulence)
        input_parameters_grid.addLayout(hbox, 4, 1, 1, 1)


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
        self.measure_abundances.setMinimumSize(
            QtCore.QSize(panel_size, 0))
        self.measure_abundances.setMaximumSize(
            QtCore.QSize(panel_size, 16777215))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.measure_abundances.setFont(font)
        self.measure_abundances.setCursor(
            QtGui.QCursor(QtCore.Qt.PointingHandCursor))
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

        # Scatter transitions.
        self.ax_first.scatter([0.5, 0.6], [0, 1], facecolor="k")
    
        # Line of best fit? Error regions?
        self.ax_first.set_xlabel("Excitation potential (eV)")
        self.ax_first.set_ylabel("[X/M]")

        self.ax_second = self.figure.figure.add_subplot(gs[1])
        self.ax_second.scatter([0.4, 0.2], [0.1, 0.3], facecolor="k")

        self.ax_second.set_xlabel(r"$\log_{e}({\rm EW}/\lambda)$")
        self.ax_second.set_ylabel("[X/M]")
        self.figure.draw()


        # Connect buttons.
        self.btn_measure_transitions.clicked.connect(self.measure_transitions)

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


    def measure_transitions(self):
        """ Trigger for when the 'Message transitions..' button is clicked. """

        # Is there any line list?
        # TODO: Don't check just for lines, check for spectral models associated
        #       with the stellar parameter determination.
        if len(self.parent.session.metadata.get("line_list", [])) == 0:

            reply = QtGui.QMessageBox.information(self,
                "No spectral models found",
                "No spectral models are currently associated with the "
                "determination of stellar parameters.\n\n"
                "Click 'OK' to load the transitions manager.")

            if reply == QtGui.QMessageBox.Ok:
                # Load line list manager.
                dialog = TransitionsDialog(self.parent.session)
                dialog.exec_()

                # Do we even have transitions now?
                # TODO: as above.
                if len(self.parent.session.metadata.get("line_list", [])) == 0:
                    return None
                else:
                    self.measure_transitions()
            else:
                return None

        else:

            print("OK show measure transitions dialog")

            # TODO HACK
            dialog = TransitionsDialog(self.parent.session)
            dialog.exec_()

    def _populate_widgets(self):
        """
        Populate the widgets in this tab with the default parameters.
        """
        return None



