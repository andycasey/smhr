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


def relim(ax, percent=0.05):

    data = ax.collections[0].get_offsets()
    x, y = data[:,0], data[:, 1]
    xlim = [
        np.min(x) - np.ptp(x) * percent,
        np.max(x) + np.ptp(x) * percent,
    ]
    ylim = [
        np.min(y) - np.ptp(y) * percent,
        np.max(y) + np.ptp(y) * percent
    ]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    return (xlim, ylim)





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
        


        # Start the grid layout for the stellar parameters tab.
        input_parameters_grid = QtGui.QGridLayout()

        """
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
        """

        # Effective temperature.
        label = QtGui.QLabel(self)
        label.setText("Effective temperature (K)")
        input_parameters_grid.addWidget(label, 1, 0, 1, 1)
        
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))

        self.edit_effective_temperature = QtGui.QLineEdit(self)
        self.edit_effective_temperature.setMaximumSize(QtCore.QSize(40, 16777215))
        self.edit_effective_temperature.setAlignment(QtCore.Qt.AlignCenter)
        self.edit_effective_temperature.setValidator(
            QtGui.QDoubleValidator(3000, 8000, 0, self.edit_effective_temperature))
        hbox.addWidget(self.edit_effective_temperature)
        input_parameters_grid.addLayout(hbox, 1, 1, 1, 1)

        # Surface gravity.
        label = QtGui.QLabel(self)
        label.setText("Surface gravity")
        input_parameters_grid.addWidget(label, 2, 0, 1, 1)
        
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))

        self.edit_surface_gravity = QtGui.QLineEdit(self)
        self.edit_surface_gravity.setMaximumSize(QtCore.QSize(40, 16777215))
        self.edit_surface_gravity.setAlignment(QtCore.Qt.AlignCenter)
        self.edit_surface_gravity.setValidator(
            QtGui.QDoubleValidator(-1, 6, 2, self.edit_surface_gravity))
        hbox.addWidget(self.edit_surface_gravity)
        input_parameters_grid.addLayout(hbox, 2, 1, 1, 1)

        # Metallicity.
        label = QtGui.QLabel(self)
        label.setText("Metallicity")
        input_parameters_grid.addWidget(label, 3, 0, 1, 1)
        
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))

        self.edit_metallicity = QtGui.QLineEdit(self)
        self.edit_metallicity.setMaximumSize(QtCore.QSize(40, 16777215))
        self.edit_metallicity.setAlignment(QtCore.Qt.AlignCenter)
        self.edit_metallicity.setValidator(
            QtGui.QDoubleValidator(-7, 1, 2, self.edit_metallicity))
        hbox.addWidget(self.edit_metallicity)
        input_parameters_grid.addLayout(hbox, 3, 1, 1, 1)

        # Depending on the photospheres: alpha-enhancement.


        # Depending on the radiative transfer code used: microturbulence.
        label = QtGui.QLabel(self)
        label.setText("Microturbulence (km/s)")
        input_parameters_grid.addWidget(label, 4, 0, 1, 1)
        
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))

        self.edit_microturbulence = QtGui.QLineEdit(self)
        self.edit_microturbulence.setMaximumSize(QtCore.QSize(40, 16777215))
        self.edit_microturbulence.setAlignment(QtCore.Qt.AlignCenter)
        self.edit_microturbulence.setValidator(
            QtGui.QDoubleValidator(0, 5, 2, self.edit_microturbulence))
        hbox.addWidget(self.edit_microturbulence)
        input_parameters_grid.addLayout(hbox, 4, 1, 1, 1)


        input_parameters_layout.addLayout(input_parameters_grid)

        # Add a spacer.
        input_parameters_layout.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding))

        
        # Add a 'Measure abundances' button.
        self.btn_measure_abundances = QtGui.QPushButton(self)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(
            self.btn_measure_abundances.sizePolicy().hasHeightForWidth())
        self.btn_measure_abundances.setSizePolicy(sp)
        self.btn_measure_abundances.setMinimumSize(
            QtCore.QSize(panel_size, 0))
        self.btn_measure_abundances.setMaximumSize(
            QtCore.QSize(panel_size, 16777215))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.btn_measure_abundances.setFont(font)
        self.btn_measure_abundances.setCursor(
            QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.btn_measure_abundances.setDefault(True)
        self.btn_measure_abundances.setObjectName("measure_abundances")
        self.btn_measure_abundances.setText("Measure abundances")
        if sys.platform == "darwin":
            self.btn_measure_abundances.setStyleSheet('QPushButton {color: white}')

        input_parameters_layout.addWidget(self.btn_measure_abundances)


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
        N = 3
        gs = gridspec.GridSpec(3, 1)
        self.ax_excitation = self.figure.figure.add_subplot(gs[0])

        # Scatter transitions.
        self.ax_excitation.scatter([], [], facecolor="k")
    
        # Line of best fit? Error regions?
        self.ax_excitation.set_xlabel("Excitation potential (eV)")
        self.ax_excitation.set_ylabel("[X/M]")

        self.ax_line_strength = self.figure.figure.add_subplot(gs[1])
        self.ax_line_strength.scatter([], [], facecolor="k")

        self.ax_line_strength.set_xlabel(r"$\log_{e}({\rm EW}/\lambda)$")
        self.ax_line_strength.set_ylabel("[X/M]")

        if N == 3:
            self.ax_opacity = self.figure.figure.add_subplot(gs[2])
            self.ax_opacity.scatter([], [], facecolor="k")
            self.ax_opacity.set_xlabel(r"Wavelength")
            self.ax_opacity.set_ylabel(r"Abundance")

        else:
            self.ax_opacity = None

        self.figure.draw()



        # Connect buttons.
        self.btn_measure_abundances.clicked.connect(self.measure_abundances)

        self.populate_widgets()

        return None


    def populate_widgets(self):
        """ Update the stellar parameter edit boxes from the session. """

        if self.parent.session is None: return None

        metadata = self.parent.session.metadata["stellar_parameters"]
        self.edit_effective_temperature.setText("{0:.0f}".format(
            metadata["effective_temperature"]))
        self.edit_surface_gravity.setText("{0:.2f}".format(
            metadata["surface_gravity"]))
        self.edit_metallicity.setText("{0:+.2f}".format(
            metadata["metallicity"]))
        self.edit_microturbulence.setText("{0:.2f}".format(
            metadata["microturbulence"]))
        return None



    def measure_abundances(self):
        """ Trigger for when the 'Message abundances' button is clicked. """

        # Are there any spectral models to be used for the determination of
        # stellar parameters?

        """
        for sm in self.parent.session.metadata.get("spectral_models", []):
            if sm.use_for_stellar_parameter_inference: break

        else:
            reply = QtGui.QMessageBox.information(self,
                "No spectral models found",
                "No spectral models are currently associated with the "
                "determination of stellar parameters.\n\n"
                "Click 'OK' to load the transitions manager.")

            if reply == QtGui.QMessageBox.Ok:
                # Load line list manager.
                dialog = TransitionsDialog(self.parent.session)
                dialog.exec_()

                # Do we even have any spectral models now?
                for sm in self.parent.session.metadata.get("spectral_models", []):
                    if sm.use_for_stellar_parameter_inference: break
                else:
                    return None
            else:
                return None
        """

        if self.parent.session is None:
            return None

        # Collate the transitions from spectral models that are profiles.
        indices = []
        equivalent_widths = []
        for model in self.parent.session.metadata["spectral_models"]:
            if not model.use_for_stellar_parameter_inference \
            or not model.is_acceptable: continue

            indices.extend(model._transition_indices)
            equivalent_widths.append(
                1e3 * model.metadata["fitted_result"][2]["equivalent_width"][0])
        indices = np.array(indices)

        # For any spectral models to be used for SPs that are not profiles,
        # re-fit them.
        transitions = self.parent.session.metadata["line_list"][indices].copy()
        transitions["equivalent_width"] = np.array(equivalent_widths)


        self.parent.session.metadata["stellar_parameters"].update({
            "effective_temperature": float(self.edit_effective_temperature.text()),
            "surface_gravity": float(self.edit_surface_gravity.text()),
            "metallicity": float(self.edit_metallicity.text()),
            "microturbulence": float(self.edit_microturbulence.text())
        })

        abundances = self.parent.session.rt.abundance_cog(
            self.parent.session.stellar_photosphere, transitions)

        print("indices", indices)
        print(abundances)

        # Update figures.
        self.ax_excitation.collections[0].set_offsets(np.array([
            transitions["expot"], abundances]).T)
        relim(self.ax_excitation)
        

        rew = np.log10(1e-3 * transitions["equivalent_width"] \
            / transitions["wavelength"])
        self.ax_line_strength.collections[0].set_offsets(
            np.array([rew, abundances]).T)
        relim(self.ax_line_strength)


        if self.ax_opacity is not None:
            self.ax_opacity.collections[0].set_offsets(np.array([
                transitions["wavelength"], abundances]).T)
            relim(self.ax_opacity)

        self.figure.draw()

        return None



