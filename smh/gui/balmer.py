#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Balmer line widget. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)


import numpy as np
import sys
from matplotlib.ticker import MaxNLocator
from PySide import QtCore, QtGui

from smh.specutils import Spectrum1D

import mpl


if sys.platform == "darwin":
        
    # See http://successfulsoftware.net/2013/10/23/fixing-qt-4-for-mac-os-x-10-9-mavericks/
    substitutes = [
        (".Lucida Grande UI", "Lucida Grande"),
        (".Helvetica Neue DeskInterface", "Helvetica Neue")
    ]
    for substitute in substitutes:
        QtGui.QFont.insertSubstitution(*substitute)



class BalmerLineFittingDialog(QtGui.QDialog):

    # If the peak-to-peak wavelength range of the observed spectrum is greater
    # than `wavelength_range_deciver`, then we will only show +/- the 
    # `wavelength_window_default` around the Balmer line.

    __wavelength_range_decider = 1000
    __wavelength_window_default = 200

    __balmer_line_names = ("H-α", "H-β", "H-γ", "H-δ")
    __balmer_line_wavelengths = (6563, 4861, 4341, 4102)

    _default_option_metadata = {
        "redshift": True,
        "smoothing": False,
        "continuum": False,
        "bounds": {}
    }

    def __init__(self, observed_spectra, observed_spectra_labels=None, 
        session=None, callbacks=None, **kwargs):
        """
        Initialise a dialog to infer stellar parameters from Balmer line models.

        :param observed_spectra:
            A list-like of observed spectra that may contain any Balmer lines.

        :param observed_spectra_labels: [optional]
            If given, this should be a list-like object the same length as
            `observed_spectra` and give labels to describe the orders in those
            spectra. If `None` is given, the observed spectra will be named as
            "Order N".

        :param session: [optional]
            If the observed_spectra are associated with a parent session, then
            providing that session here will allow the results of any Balmer-
            line fitting to be saved to the session.

        :param callbacks: [optional]
            A list-like of functions to execute when the dialog is closed.
        """

        super(BalmerLineFittingDialog, self).__init__(**kwargs)
        
        if isinstance(observed_spectra, Spectrum1D):
            observed_spectra = [observed_spectra]

        elif not isinstance(observed_spectra, (list, tuple, np.ndarray)) \
        or not all([isinstance(s, Spectrum1D) for s in observed_spectra]):
            raise TypeError(
                "observed spectra must be a Spectrum1D "
                "or a list-like of Spectrum1D")


        if observed_spectra_labels is not None:
            if len(observed_spectra_labels) != len(observed_spectra):
                raise ValueError("number of observed spectra does not match the"
                    " number of observed spectra labels ({} != {})".format(
                        len(observed_spectra_labels), len(observed_spectra)))

            if len(set(observed_spectra_labels)) < len(observed_spectra_labels):
                raise ValueError("observed spectra labels must be unique")

        else:
            observed_spectra_labels = ["Order {}".format(i) \
                for i in range(1, 1 + len(observed_spectra))]

        # Save information to object.
        self.observed_spectra = observed_spectra
        self.observed_spectra_labels = observed_spectra_labels
        self.callbacks = callbacks or []
        self.session = session

        # Identify Balmer lines in the given data.
        self._identify_balmer_lines()


        # Start creating GUI.
        self.setGeometry(800, 600, 800, 600)
        self.move(QtGui.QApplication.desktop().screen().rect().center() \
            - self.rect().center())
        self.setWindowTitle("Balmer-line fitting")

        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, 
            QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)


        # Panel 1
        self.p1_vbox = QtGui.QVBoxLayout(self)

        hbox = QtGui.QHBoxLayout()
        self.combo_balmer_line_selected = QtGui.QComboBox(self)

        hbox.addWidget(self.combo_balmer_line_selected)
        hbox.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Minimum))

        self.combo_spectrum_selected = QtGui.QComboBox(self)
        hbox.addWidget(self.combo_spectrum_selected)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.p1_vbox.addLayout(hbox)

        line = QtGui.QFrame(self)
        line.setFrameShape(QtGui.QFrame.HLine)
        line.setFrameShadow(QtGui.QFrame.Sunken)
        self.p1_vbox.addWidget(line)

        # Model options table and group.
        hbox = QtGui.QHBoxLayout()
        gp = QtGui.QGroupBox(self)
        gp.setTitle("Model options")
        gp.setMinimumSize(QtCore.QSize(0, 0))
        vbox = QtGui.QVBoxLayout(gp)

        self.p1_model_options = QtGui.QTableView(gp)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Expanding)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(self.p1_model_options.sizePolicy().hasHeightForWidth())
        self.p1_model_options.setSizePolicy(sp)
        self.p1_model_options.setMinimumSize(QtCore.QSize(300, 16777215))
        self.p1_model_options.setMaximumSize(QtCore.QSize(300, 16777215))
        self.p1_model_options.setEditTriggers(
            QtGui.QAbstractItemView.CurrentChanged)
        self.p1_model_options.setModel(BalmerLineOptionsTableModel(self))
        self.p1_model_options.horizontalHeader().setResizeMode(
            QtGui.QHeaderView.Stretch)

        # First column should be fixed for the checkbox.
        self.p1_model_options.horizontalHeader().setResizeMode(
            0, QtGui.QHeaderView.Fixed)
        self.p1_model_options.horizontalHeader().resizeSection(0, 30) # MAGIC


        vbox.addWidget(self.p1_model_options)
        hbox.addWidget(gp)

        self.p1_figure = mpl.MPLWidget(None, tight_layout=True, matchbg=self)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(self.p1_figure.sizePolicy().hasHeightForWidth())
        self.p1_figure.setSizePolicy(sp)
        hbox.addWidget(self.p1_figure)


        self.p1_vbox.addLayout(hbox)


        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))

        self.p1_btn_next = QtGui.QPushButton(self)
        self.p1_btn_next.setText("Next")
        self.p1_btn_next.setFocusPolicy(QtCore.Qt.NoFocus)
        hbox.addWidget(self.p1_btn_next)
        
        self.p1_vbox.addLayout(hbox)


        # Initialize widgets that do not depend on the input spectra.
        for name, wavelength \
        in zip(self.__balmer_line_wavelengths, self.__balmer_line_names):
            self.combo_balmer_line_selected.addItem(
                u"{} ({} Å)".format(name, wavelength))

        self.combo_balmer_line_selected.setFocus()

        ax = self.p1_figure.figure.add_subplot(111)
        ax.plot([], [], c="k", drawstyle="steps-mid", zorder=10)
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.set_xlabel(u"Wavelength [Å]")
        ax.set_ylabel(u"Flux")

        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        # Enable drag to mask regions.
        self.p1_figure.enable_drag_to_mask(ax)

        self.metadata = {}
        self.metadata.update(self._default_option_metadata)

        # Signals.
        self.combo_balmer_line_selected.currentIndexChanged.connect(
            self.updated_balmer_line_selected)
        self.combo_spectrum_selected.currentIndexChanged.connect(
            self.updated_spectrum_selected)

        # Populate widgets.
        self.populate_widgets()
        
        return None




    def _identify_balmer_lines(self):
        """ Identify the Balmer lines in the spectra provided. """

        spectra_indices = []
        for wavelength in self.__balmer_line_wavelengths:

            spectrum_indices = []
            # Look for this wavelength in all of the spectra we have.
            for i, spec in enumerate(self.observed_spectra):    
                if spec.dispersion.max() >= wavelength >= spec.dispersion.min():
                    spectrum_indices.append(i)
            spectra_indices.append(spectrum_indices)

        self._balmer_line_indices = spectra_indices
        return spectra_indices



    def populate_widgets(self):
        """ Populate widgets based on information available in the session. """

        # Which Balmer lines are available, and in which spectral orders?
        first_available = None
        for i, indices in enumerate(self._balmer_line_indices):

            # If there are no orders containing this Balmer line, set the option
            # to disabled.
            item = self.combo_balmer_line_selected.model().item(i)
            is_available = len(indices) > 0
            item.setEnabled(is_available)

            if is_available:
                first_available = first_available or i

        # Select the first available Balmer line, which will trigger the
        # spectra available to be updated.
        if first_available is not None:
            self.combo_balmer_line_selected.setCurrentIndex(first_available)



        return None


    def updated_balmer_line_selected(self):
        """ The Balmer line selected has been updated. """

        # If we had already selected a spectrum (e.g., a rest-frame normalized
        # one) and the new Balmer line is in that spectrum too, we should keep
        # that selection.
        current_spectrum_selected = self.combo_spectrum_selected.currentText()

        self.combo_spectrum_selected.clear()

        selected_spectrum_index = None
        selected_balmer_index = self.combo_balmer_line_selected.currentIndex()
        for j, idx in enumerate(self._balmer_line_indices[selected_balmer_index]):

            spectrum_text = self.observed_spectra_labels[idx]
            self.combo_spectrum_selected.addItem(spectrum_text)

            if spectrum_text == current_spectrum_selected:
                selected_spectrum_index = j

        self.combo_spectrum_selected.setCurrentIndex(
            selected_spectrum_index or 0)

        return None


    def updated_spectrum_selected(self):
        """ The spectrum to use for model fitting has been updated. """

        # Get the spectrum.
        try:
            spectrum_index = self.observed_spectra_labels.index(
                self.combo_spectrum_selected.currentText())

        except ValueError:
            return None

        spectrum = self.observed_spectra[spectrum_index]

        # Either show the entire spectrum, or just around the Balmer line if the
        # order goes on for >1000 Angstroms.
        balmer_line_index = self.combo_balmer_line_selected.currentIndex()
        if np.ptp(spectrum.dispersion) > self.__wavelength_range_decider:

            balmer_wavelength = self.__balmer_line_wavelengths[balmer_line_index]
            window = self.__wavelength_window_default
            limits = (balmer_wavelength - window, balmer_wavelength + window)

        else:
            limits = (spectrum.dispersion.min(), spectrum.dispersion.max())

        view_mask \
            = (limits[1] >= spectrum.dispersion) * (spectrum.dispersion >= limits[0])

        # Update the spectrum figure.
        ax = self.p1_figure.figure.axes[0]
        ax.lines[0].set_data(np.array([
            spectrum.dispersion[view_mask],
            spectrum.flux[view_mask],
        ]))

        ax.set_xlim(limits)
        ax.set_ylim(0, 1.1 * np.nanmax(spectrum.flux[view_mask]))

        ax.xaxis.set_visible(True)
        ax.yaxis.set_visible(True)

        # TODO: Draw uncertainties.

        # TODO: Draw masks.
        self.p1_figure.draw()

        return None






class BalmerLineOptionsTableModel(QtCore.QAbstractTableModel):

    _max_continuum_order = 30
    _parameters = ["redshift", "smoothing", "continuum"]

    def __init__(self, parent, *args):
        super(BalmerLineOptionsTableModel, self).__init__(parent, *args)
        self.parent = parent


    def rowCount(self, parent):
        return len(self._parameters)

    def columnCount(self, parent):
        return 3


    def setData(self, index, value, role):

        if index.column() == 0:
            try:
                parameter = self._parameters[index.row()]

            except IndexError:
                return False

            value = bool(value)

            # Continuum is a special case.
            if index.row() == 2:
                if value:

                    N, is_ok = QtGui.QInputDialog.getItem(None, 
                        "Continuum order", "Specify continuum order:", 
                        ["{:.0f}".format(i) \
                            for i in range(1 + self._max_continuum_order)])

                    if not is_ok or int(N) > self._max_continuum_order:
                        return False

                    self._parameters.extend(
                        ["c{}".format(i) for i in range(1 + int(N))])

                else:
                    self._parameters = self._parameters[:3]

                self.reset()

            self.parent.metadata[parameter] = bool(value)
            return True

        else:
            parameter = self._parameters[index.row()]
            self.parent.metadata["bounds"].setdefault(parameter, [None, None])

            try:
                value = float(value)

            except:
                return False

            else:
                # Check bound is less than other bound.
                bounds = self.parent.metadata["bounds"][parameter]
                if (index.column() == 1 and bounds[1] is not None \
                    and value >= bounds[1]) \
                or (index.column() == 2 and bounds[0] is not None \
                    and value <= bounds[0]):
                    raise ValueError("bounds must be (lower, upper)")

                if not np.isfinite(value):
                    return False

                self.parent.metadata["bounds"][parameter][index.column() - 1] \
                    = value
                
            return True

        return False


    def headerData(self, index, orientation, role):

        if orientation == QtCore.Qt.Vertical \
        and role == QtCore.Qt.DisplayRole:

            translation = {
                "redshift": "Radial velocity [km/s]",
                "smoothing": "Macroturbulence [km/s]",
                "continuum": "Continuum",
            }.get(self._parameters[index], None)

            return translation or "        c{}".format(index - 3)

        elif orientation == QtCore.Qt.Horizontal \
        and role == QtCore.Qt.DisplayRole:
            return ["", "Lower\nbound", "Upper\nbound"][index]
        return None


    def flags(self, index):
        if not index.isValid():
            return None

        if index.column() == 0 and index.row() < 3:
            return  QtCore.Qt.ItemIsEnabled|\
                    QtCore.Qt.ItemIsUserCheckable

        elif index.column() == 0 and index.row() >= 3:
            return QtCore.Qt.NoItemFlags

        else:
            parameter = self._parameters[index.row()]
            if index.row() > 2 or self.parent.metadata[parameter]:
                return QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable
            else:
                return QtCore.Qt.NoItemFlags
        return None



    def data(self, index, role):
        if not index.isValid():
            return None

        if index.column() == 0 and index.row() < 3 \
        and role in (QtCore.Qt.DisplayRole, QtCore.Qt.CheckStateRole):

            value = self.parent.metadata[self._parameters[index.row()]]
            if role == QtCore.Qt.CheckStateRole:
                return QtCore.Qt.Checked if value else QtCore.Qt.Unchecked

            else:
                return None
        
        elif index.column() == 0 and index.row() >= 3 \
        and role == QtCore.Qt.DisplayRole:
            return ""


        if role != QtCore.Qt.DisplayRole:
            return None

        # Look for bound information.
        bounds = self.parent.metadata["bounds"].get(
            self._parameters[index.row()], (None, None))

        value = bounds[index.column() - 1]
        if value is None:
            return "None"

        elif value == np.inf:
            return u"+inf"

        elif value == -np.inf:
            return u"-inf"

        else:
            return "{}".format(value)

        return "None"










if __name__ == "__main__":

    import sys

    # This is just for development testing.
    try:
        app = QtGui.QApplication(sys.argv)

    except RuntimeError:
        None

    # Load some spectra.
    spectra = [] + \
        Spectrum1D.read("/Users/arc/Downloads/hd122563_1blue_multi_090205_oldbutgood.fits") + \
        Spectrum1D.read("/Users/arc/Downloads/hd122563_1red_multi_090205_oldbutgood.fits") + \
        [Spectrum1D.read("smh/balmer/hd122563.fits")]


    window = BalmerLineFittingDialog(spectra,
        observed_spectra_labels=["Order {}".format(i) for i in range(1, len(spectra))] + ["Normalized rest-frame spectrum"])

    window.exec_()

    



