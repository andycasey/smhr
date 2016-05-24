#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The stellar parameters tab in Spectroscopy Made Hard """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["StellarParametersTab"]

import logging
import matplotlib.gridspec
import numpy as np
import sys
from PySide import QtCore, QtGui
from matplotlib.ticker import MaxNLocator

import mpl, style_utils
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


        self.parent_layout = QtGui.QHBoxLayout(self)
        lhs_layout = QtGui.QVBoxLayout()
        grid_layout = QtGui.QGridLayout()

        # Effective temperature.
        label = QtGui.QLabel(self)
        label.setText("Effective temperature")
        grid_layout.addWidget(label, 0, 0, 1, 1)
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, 
            QtGui.QSizePolicy.Minimum))
        self.edit_teff = QtGui.QLineEdit(self)
        self.edit_teff.setMinimumSize(QtCore.QSize(40, 0))
        self.edit_teff.setMaximumSize(QtCore.QSize(50, 16777215))
        self.edit_teff.setValidator(
            QtGui.QDoubleValidator(3000, 8000, 0, self.edit_teff))

        hbox.addWidget(self.edit_teff)
        grid_layout.addLayout(hbox, 0, 1, 1, 1)
        
        # Surface gravity.
        label = QtGui.QLabel(self)
        label.setText("Surface gravity")
        grid_layout.addWidget(label, 1, 0, 1, 1)
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, 
            QtGui.QSizePolicy.Minimum))
        self.edit_logg = QtGui.QLineEdit(self)
        self.edit_logg.setMinimumSize(QtCore.QSize(40, 0))
        self.edit_logg.setMaximumSize(QtCore.QSize(50, 16777215))
        self.edit_logg.setValidator(
            QtGui.QDoubleValidator(-1, 6, 3, self.edit_logg))
        hbox.addWidget(self.edit_logg)
        grid_layout.addLayout(hbox, 1, 1, 1, 1)

        # Metallicity.
        label = QtGui.QLabel(self)
        label.setText("Metallicity")
        grid_layout.addWidget(label, 2, 0, 1, 1)
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, 
            QtGui.QSizePolicy.Minimum))
        self.edit_metallicity = QtGui.QLineEdit(self)
        self.edit_metallicity.setMinimumSize(QtCore.QSize(40, 0))
        self.edit_metallicity.setMaximumSize(QtCore.QSize(50, 16777215))
        self.edit_metallicity.setValidator(
            QtGui.QDoubleValidator(-5, 1, 3, self.edit_metallicity))
        hbox.addWidget(self.edit_metallicity)
        grid_layout.addLayout(hbox, 2, 1, 1, 1)


        # Microturbulence.
        label = QtGui.QLabel(self)
        label.setText("Microturbulence")
        grid_layout.addWidget(label, 3, 0, 1, 1)
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, 
            QtGui.QSizePolicy.Minimum))
        self.edit_xi = QtGui.QLineEdit(self)
        self.edit_xi.setMinimumSize(QtCore.QSize(40, 0))
        self.edit_xi.setMaximumSize(QtCore.QSize(50, 16777215))
        self.edit_xi.setValidator(QtGui.QDoubleValidator(0, 5, 3, self.edit_xi))
        hbox.addWidget(self.edit_xi)
        grid_layout.addLayout(hbox, 3, 1, 1, 1)

        # Optionally TODO: alpha-enhancement.

        lhs_layout.addLayout(grid_layout)

        # Buttons for solving/measuring.        
        hbox = QtGui.QHBoxLayout()
        self.btn_measure = QtGui.QPushButton(self)
        self.btn_measure.setAutoDefault(False)
        self.btn_measure.setDefault(True)
        self.btn_measure.setText("Measure")
        hbox.addWidget(self.btn_measure)

        self.btn_options = QtGui.QPushButton(self)
        self.btn_options.setText("Options..")
        hbox.addWidget(self.btn_options)

        hbox.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, 
            QtGui.QSizePolicy.Minimum))

        self.btn_solve = QtGui.QPushButton(self)
        self.btn_solve.setText("Solve")
        hbox.addWidget(self.btn_solve)
        lhs_layout.addLayout(hbox)

        line = QtGui.QFrame(self)
        line.setFrameShape(QtGui.QFrame.HLine)
        line.setFrameShadow(QtGui.QFrame.Sunken)
        lhs_layout.addWidget(line)

        self.spectral_models = []
        self.table_view = QtGui.QTableView(self)
        self.table_view.setModel(SpectralModelsTableModel(self))
        self.table_view.setSelectionBehavior(
            QtGui.QAbstractItemView.SelectRows)
        self.table_view.setSortingEnabled(True)
        self.table_view.resizeColumnsToContents()
        self.table_view.setColumnWidth(0, 30) # MAGIC
        self.table_view.horizontalHeader().setStretchLastSection(True)


        lhs_layout.addWidget(self.table_view)
        self.parent_layout.addLayout(lhs_layout)


        # Matplotlib figure.
        self.figure = mpl.MPLWidget(None, tight_layout=True, autofocus=True)
        self.figure.setMinimumSize(QtCore.QSize(800, 300))
        self.figure.figure.patch.set_facecolor([(_ - 10)/255. for _ in \
            self.palette().color(QtGui.QPalette.Window).getRgb()[:3]])

        gs_top = matplotlib.gridspec.GridSpec(4, 1)
        gs_bottom = matplotlib.gridspec.GridSpec(4, 1, 
            height_ratios=[2, 2, 1, 2])
        gs_bottom.update(hspace=0)

        self.ax_excitation = self.figure.figure.add_subplot(gs_top[0])
        self.ax_excitation.scatter([], [], facecolor="k")
        self.ax_excitation.xaxis.get_major_formatter().set_useOffset(False)
        self.ax_excitation.yaxis.set_major_locator(MaxNLocator(5))
        self.ax_excitation.yaxis.set_major_locator(MaxNLocator(3))
        self.ax_excitation.set_xlabel("Excitation potential (eV)")
        self.ax_excitation.set_ylabel("[X/M]")

        self.ax_line_strength = self.figure.figure.add_subplot(gs_top[1])
        self.ax_line_strength.scatter([], [], facecolor="k")
        self.ax_line_strength.xaxis.get_major_formatter().set_useOffset(False)
        self.ax_line_strength.yaxis.set_major_locator(MaxNLocator(5))
        self.ax_line_strength.yaxis.set_major_locator(MaxNLocator(3))
        self.ax_line_strength.set_xlabel(r"$\log({\rm EW}/\lambda)$")
        self.ax_line_strength.set_ylabel("[X/M]")

        self.ax_residuals = self.figure.figure.add_subplot(gs_bottom[2])
        self.ax_residuals.axhline(0, c="#666666")
        self.ax_residuals.set_xticklabels([])
        self.ax_residuals.yaxis.set_major_locator(MaxNLocator(2))
        self.ax_residuals.set_ylabel(r"$\Delta$")

        self.ax_spectrum = self.figure.figure.add_subplot(gs_bottom[3])
        self.ax_spectrum.xaxis.get_major_formatter().set_useOffset(False)
        self.ax_spectrum.xaxis.set_major_locator(MaxNLocator(5))
        self.ax_spectrum.yaxis.set_major_locator(MaxNLocator(3))
        self.ax_spectrum.set_xlabel(u"Wavelength (Å)")
        self.ax_spectrum.set_ylabel(r"Normalized flux")

        self.parent_layout.addWidget(self.figure)

        # Connect buttons.
        self.btn_measure.clicked.connect(self.measure_abundances)
        self.btn_options.clicked.connect(self.options)
        self.btn_solve.clicked.connect(self.solve_parameters)

        return None


    def populate_widgets(self):
        """ Update the stellar parameter edit boxes from the session. """

        if not hasattr(self.parent, "session") or self.parent.session is None:
            return None

        metadata = self.parent.session.metadata["stellar_parameters"]
        self.edit_teff.setText("{0:.0f}".format(
            metadata["effective_temperature"]))
        self.edit_logg.setText("{0:.2f}".format(
            metadata["surface_gravity"]))
        self.edit_metallicity.setText("{0:+.2f}".format(
            metadata["metallicity"]))
        self.edit_microturbulence.setText("{0:.2f}".format(
            metadata["microturbulence"]))
        return None


    def update_stellar_parameters(self):
        """ Update the stellar parameters with the values in the GUI. """

        self.parent.session.metadata["stellar_parameters"].update({
            "effective_temperature": float(self.edit_teff.text()),
            "surface_gravity": float(self.edit_logg.text()),
            "metallicity": float(self.edit_metallicity.text()),
            "microturbulence": float(self.edit_xi.text())
        })
        return True



    def _check_for_spectral_models(self):
        """
        Check the session for any valid spectral models that are associated with
        the determination of stellar parameters.
        """

        # Are there any spectral models to be used for the determination of
        # stellar parameters?
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
                    return False
            else:
                return False

        return True


    def measure_abundances(self):
        """ The measure abundances button has been clicked. """

        if self.parent.session is None or not self._check_for_spectral_models():
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

        # Update the session with the stellar parameters in the GUI.
        self.update_stellar_parameters()

        # Calculate abundances.
        abundances = self.parent.session.rt.abundance_cog(
            self.parent.session.stellar_photosphere, transitions)

        # Update figures.
        self.ax_excitation.collections[0].set_offsets(
            np.array([transitions["expot"], abundances]).T)

        rew = 1e-3 * transitions["equivalent_width"]/transitions["wavelength"]
        self.ax_line_strength.collections[0].set_offsets(
            np.array([np.log10(rew), abundances]).T)
        
        # Update limits on the excitation and line strength figures.
        style_utils.relim(self.ax_excitation)
        style_utils.relim(self.ax_line_strength)

        self.figure.draw()

        # Update the table?
        return None


    def options(self):
        """ Open a GUI for the radiative transfer and solver options. """
        raise NotImplementedError


    def solve_parameters(self):
        """ Solve the stellar parameters. """
        raise NotImplementedError




class SpectralModelsTableModel(QtCore.QAbstractTableModel):

    header = ["", u"λ\n(Å)", "Element\n", u"E. W.\n(mÅ)",
        "log ε\n(dex)"]
    attrs = ("is_acceptable", "_repr_wavelength", "_repr_element", 
        "equivalent_width", "abundance")

    def __init__(self, parent, *args):
        """
        An abstract table model for spectral models.

        :param parent:
            The parent widget.
        """

        super(SpectralModelsTableModel, self).__init__(parent, *args)
        self._parent = parent
        return None


    def rowCount(self, parent):
        """ Return the number of rows in the table. """
        return len(self._parent.spectral_models)


    def columnCount(self, parent):
        """ Return the number of columns in the table. """
        return len(self.header)


    def data(self, index, role):
        """
        Display the data.

        :param index:
            The table index.

        :param role:
            The display role.
        """

        if not index.isValid():
            return None

        column = index.column()
        if  column == 0 \
        and role in (QtCore.Qt.DisplayRole, QtCore.Qt.CheckStateRole):
            value = self._parent.spectral_models[index.row()].is_acceptable
            if role == QtCore.Qt.CheckStateRole:
                return QtCore.Qt.Checked if value else QtCore.Qt.Unchecked
            else:
                return None

        elif column == 1:
            value = self._parent.spectral_models[index.row()]._repr_wavelength

        elif column == 2:
            value = self._parent.spectral_models[index.row()]._repr_element

        elif column == 3:
            try:
                spectral_model = self._parent.spectral_models[index.row()]
                result = spectral_model.metadata["fitted_result"][2]
                equivalent_width = result["equivalent_width"][0]
            except:
                equivalent_width = np.nan

            value = "{0:.1f}".format(1000 * equivalent_width) \
                if np.isfinite(equivalent_width) else ""

        elif column == 4:
            try:
                spectral_model = self._parent.spectral_models[index.row()]
                abundances \
                    = spectral_model.metadata["fitted_result"][2]["abundances"]

            except (IndexError, KeyError):
                value = ""

            else:
                # How many elements were measured?
                value = "; ".join(["{0:.2f}".format(abundance) \
                    for abundance in abundances])

        return value if role == QtCore.Qt.DisplayRole else None


    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal \
        and role == QtCore.Qt.DisplayRole:
            return self.header[col]
        return None


    def setData(self, index, value, role=QtCore.Qt.DisplayRole):
        if index.column() != 0 or value:
            return False

        model = self._parent.spectral_models[index.row()]
        model.metadata["is_acceptable"] = value

        # Emit data change for this row.
        self.dataChanged.emit(
            self.createIndex(index.row(), 0),
            self.createIndex(index.row(), self.columnCount(0)))

        # Refresh the view.
        try:
            del model.metadata["fitted_result"]

        except KeyError:
            None

        self._parent.redraw_figure()

        return value

    """
    # TODO
    def sort(self, column, order):
        raise NotImplementedError
        self.emit(QtCore.SIGNAL("layoutAboutToBeChanged()"))
        sorter = {
            0: lambda model: model.is_acceptable,
            1: lambda model: model._repr_wavelength,
            2: lambda model: model._repr_element,
            3: lambda model: model.metadata.get("fitted_result", [])[-1]
        }

        self.spectral_models = sorted(self.spectral_models,
            key=lambda sm: getattr(sm, self.attrs[column]))
        
        if order == QtCore.Qt.DescendingOrder:
            self.spectral_models.reverse()

        self.dataChanged.emit(self.createIndex(0, 0),
            self.createIndex(self.rowCount(0), self.columnCount(0)))
        self.emit(QtCore.SIGNAL("layoutChanged()"))
    """


    def flags(self, index):
        if not index.isValid(): return
        return  QtCore.Qt.ItemIsSelectable|\
                QtCore.Qt.ItemIsEnabled|\
                QtCore.Qt.ItemIsUserCheckable


