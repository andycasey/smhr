#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A widget for fitting spectral models. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import operator
import numpy as np
from PySide import QtCore, QtGui
from time import time
from matplotlib.ticker import MaxNLocator

import mpl
import style_utils

from smh.spectral_models import (ProfileFittingModel, SpectralSynthesisModel)

DOUBLE_CLICK_INTERVAL = 0.1 # MAGIC HACK


class SpectralModelsTableModel(QtCore.QAbstractTableModel):

    header = [" ", u"Wavelength\n(Å)", "Element", u"E. W.\n(mÅ)",
        "Abundance\n(dex)"]
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
        if index.column() != 0:
            return False

        self._parent.spectral_models[index.row()].metadata["is_acceptable"] \
            = value
        self.dataChanged.emit(index, index)
        return value

    """
    def sort(self, column, order):

        self.emit(QtCore.SIGNAL("layoutAboutToBeChanged()"))
        self._data = sorted(self._data,
            key=lambda sm: getattr(sm, self.attrs[column]))
        
        if order == QtCore.Qt.DescendingOrder:
            self._data.reverse()

        self.dataChanged.emit(self.createIndex(0, 0),
            self.createIndex(self.rowCount(0), self.columnCount(0)))
        self.emit(QtCore.SIGNAL("layoutChanged()"))
    """

    def flags(self, index):
        if not index.isValid(): return
        return  QtCore.Qt.ItemIsSelectable|\
                QtCore.Qt.ItemIsEnabled|\
                QtCore.Qt.ItemIsUserCheckable





class SpectralModelsDialog(QtGui.QDialog):

    def __init__(self, session, stellar_parameter_models=True,
        stellar_abundance_models=True, *args):
        super(SpectralModelsDialog, self).__init__(*args)

        self.session = session

        # Create a list of the spectral models that we will display here.
        self.spectral_models = []
        for spectral_model in self.session.metadata.get("spectral_models", []):
            if (stellar_parameter_models and 
                spectral_model.use_for_stellar_parameter_inference) \
            or (stellar_abundance_models and
                spectral_model.use_for_stellar_composition_inference):
                self.spectral_models.append(spectral_model)

        self.setGeometry(1000, 500, 1000, 500)
        self.move(QtGui.QApplication.desktop().screen().rect().center() \
            - self.rect().center())
        self.setWindowTitle("Fit spectral models")

        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, 
            QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)

        hbox_parent = QtGui.QHBoxLayout(self)

        self.table_view = QtGui.QTableView(self)
        self.table_view.setMinimumSize(QtCore.QSize(300, 0))
        self.table_view.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.table_view.setModel(SpectralModelsTableModel(self))
        self.table_view.setSelectionBehavior(
            QtGui.QAbstractItemView.SelectRows)
        self.table_view.setSortingEnabled(True)
        self.table_view.resizeColumnsToContents()
        self.table_view.setColumnWidth(0, 30) # MAGIC
        self.table_view.setColumnWidth(1, 80) # MAGIC
        self.table_view.setColumnWidth(3, 60) # MAGIC

        self.table_view.horizontalHeader().setStretchLastSection(True)


        _ = self.table_view.selectionModel()
        _.selectionChanged.connect(self.selected_model_changed)

        hbox_parent.addWidget(self.table_view)


        vbox_rhs = QtGui.QVBoxLayout()

        self.mpl_figure = mpl.MPLWidget(None, tight_layout=True,
            autofocus=True)
        self.mpl_figure.figure.patch.set_facecolor([_/255. for _ in \
            self.palette().color(QtGui.QPalette.Window).getRgb()[:3]])
        self.mpl_figure.setMinimumSize(QtCore.QSize(0, 250))
        vbox_rhs.addWidget(self.mpl_figure)


        # Model options.
        group_box = QtGui.QGroupBox(self)
        group_box.setTitle("Model options")
        self.verticalLayout = QtGui.QVBoxLayout(group_box)
        self.verticalLayout.setContentsMargins(6, 12, 6, 6)
        self.tabWidget = QtGui.QTabWidget(group_box)

        # Common model options.
        self.tab_common = QtGui.QWidget()
        vbox_common = QtGui.QVBoxLayout(self.tab_common)
        grid_common = QtGui.QGridLayout()
        grid_common.addItem(
            QtGui.QSpacerItem(40, 20, 
                QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum),
            1, 2, 1, 1)

        label = QtGui.QLabel(self.tab_common)
        label.setText("Data fitting window")
        grid_common.addWidget(label, 0, 1, 1, 1)
        self.edit_window = QtGui.QLineEdit(self.tab_common)
        self.edit_window.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_window.setMaximumSize(QtCore.QSize(60, 16777215))
        self.edit_window.setValidator(
            QtGui.QDoubleValidator(0, 1000, 1, self.edit_window))
        grid_common.addWidget(self.edit_window, 0, 3, 1, 1)

        self.checkbox_continuum = QtGui.QCheckBox(self.tab_common)
        self.checkbox_continuum.setText("")
        grid_common.addWidget(self.checkbox_continuum, 1, 0, 1, 1)
        label = QtGui.QLabel(self.tab_common)
        label.setText("Model continuum with polynomial order")
        grid_common.addWidget(label, 1, 1, 1, 1)
        self.combo_continuum = QtGui.QComboBox(self.tab_common)
        self.combo_continuum.setMinimumSize(QtCore.QSize(60, 0))
        self.combo_continuum.setMaximumSize(QtCore.QSize(60, 16777215))
        grid_common.addWidget(self.combo_continuum, 1, 3, 1, 1)


        for i in range(10):
            self.combo_continuum.addItem("{:.0f}".format(i))


        self.checkbox_vrad_tolerance = QtGui.QCheckBox(self.tab_common)
        self.checkbox_vrad_tolerance.setText("")
        grid_common.addWidget(self.checkbox_vrad_tolerance, 2, 0, 1, 1)
        label = QtGui.QLabel(self.tab_common)
        label.setText("Set tolerance on residual radial velocity")
        grid_common.addWidget(label, 2, 1, 1, 1)
        self.edit_vrad_tolerance = QtGui.QLineEdit(self.tab_common)
        self.edit_vrad_tolerance.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_vrad_tolerance.setMaximumSize(QtCore.QSize(60, 16777215))
        self.edit_vrad_tolerance.setValidator(
            QtGui.QDoubleValidator(0, 100, 2, self.edit_vrad_tolerance))
        grid_common.addWidget(self.edit_vrad_tolerance, 2, 3, 1, 1)

        grid_common.addItem(QtGui.QSpacerItem(40, 20, 
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum), 2, 2, 1, 1)
        vbox_common.addLayout(grid_common)
        vbox_common.addItem(QtGui.QSpacerItem(20, 40, 
            QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding))
        self.tabWidget.addTab(self.tab_common, "Common")


        # Profile model options.
        self.tab_profile = QtGui.QWidget()
        vbox_profile = QtGui.QVBoxLayout(self.tab_profile)
        grid_profile = QtGui.QGridLayout()


        label = QtGui.QLabel(self.tab_profile)
        label.setText("Profile type")
        grid_profile.addWidget(label, 0, 1, 1, 1)
        grid_profile.addItem(
            QtGui.QSpacerItem(40, 20, 
                QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum),
            0, 2, 1, 1)
        self.combo_profile = QtGui.QComboBox(self.tab_profile)
        grid_profile.addWidget(self.combo_profile, 0, 3, 1, 1)

        for each in ("Gaussian", "Lorentzian", "Voigt"):
            self.combo_profile.addItem(each)

        label = QtGui.QLabel(self.tab_profile)
        label.setText("Detection sigma for nearby absorption lines")
        grid_profile.addWidget(label, 1, 1, 1, 1)
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(40, 20, 
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.edit_detection_sigma = QtGui.QLineEdit(self.tab_profile)
        self.edit_detection_sigma.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_detection_sigma.setMaximumSize(QtCore.QSize(60, 16777215))
        self.edit_detection_sigma.setValidator(
            QtGui.QDoubleValidator(0, 100, 1, self.edit_detection_sigma))
        hbox.addWidget(self.edit_detection_sigma)
        grid_profile.addLayout(hbox, 1, 3, 1, 1)


        label = QtGui.QLabel(self.tab_profile)
        label.setText("Neighbouring pixels required to detect nearby lines")
        grid_profile.addWidget(label, 2, 1, 1, 1)
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(40, 20, 
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.edit_detection_pixels = QtGui.QLineEdit(self.tab_profile)
        self.edit_detection_pixels.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_detection_pixels.setMaximumSize(QtCore.QSize(60, 16777215))
        self.edit_detection_pixels.setValidator(
            QtGui.QIntValidator(0, 100, self.edit_detection_pixels))
        hbox.addWidget(self.edit_detection_pixels)
        grid_profile.addLayout(hbox, 2, 3, 1, 1)


        label = QtGui.QLabel(self.tab_profile)
        label.setText("Use central pixel weighting")
        grid_profile.addWidget(label, 4, 1, 1, 1)
        self.checkbox_use_central_weighting = QtGui.QCheckBox(self.tab_profile)
        self.checkbox_use_central_weighting.setText("")
        grid_profile.addWidget(self.checkbox_use_central_weighting, 4, 0, 1, 1)

    
        label = QtGui.QLabel(self.tab_profile)
        label.setText("Tolerance in wavelength position")
        grid_profile.addWidget(label, 3, 1, 1, 1)
        self.checkbox_wavelength_tolerance = QtGui.QCheckBox(self.tab_profile)
        self.checkbox_wavelength_tolerance.setText("")
        grid_profile.addWidget(self.checkbox_wavelength_tolerance, 3, 0, 1, 1)
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum))
        self.edit_wavelength_tolerance = QtGui.QLineEdit(self.tab_profile)
        self.edit_wavelength_tolerance.setMinimumSize(QtCore.QSize(50, 0))
        self.edit_wavelength_tolerance.setMaximumSize(QtCore.QSize(60, 16777215))
        self.edit_wavelength_tolerance.setValidator(
            QtGui.QDoubleValidator(0, 10, 2, self.edit_wavelength_tolerance))
        hbox.addWidget(self.edit_wavelength_tolerance)
        grid_profile.addLayout(hbox, 3, 3, 1, 1)

        vbox_profile.addLayout(grid_profile)
        vbox_profile.addItem(QtGui.QSpacerItem(20, 10, 
            QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding))
        self.tabWidget.addTab(self.tab_profile, "Profile options")


        # Synthesis model options.
        self.tab_synthesis = QtGui.QWidget()
        vbox_synthesis = QtGui.QVBoxLayout(self.tab_synthesis)
        grid_synthesis = QtGui.QGridLayout()

        label = QtGui.QLabel(self.tab_synthesis)
        label.setText("Initial abundance boundary")
        grid_synthesis.addWidget(label, 0, 1, 1, 1)
        self.edit_initial_abundance_bound = QtGui.QLineEdit(self.tab_synthesis)
        self.edit_initial_abundance_bound.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_initial_abundance_bound.setMaximumSize(QtCore.QSize(60, 16777215))
        self.edit_initial_abundance_bound.setValidator(
            QtGui.QDoubleValidator(0, 2, 1, self.edit_initial_abundance_bound))
        grid_synthesis.addWidget(self.edit_initial_abundance_bound, 0, 3, 1, 1)

        self.checkbox_model_smoothing = QtGui.QCheckBox(self.tab_synthesis)
        self.checkbox_model_smoothing.setText("")
        grid_synthesis.addWidget(self.checkbox_model_smoothing, 1, 0, 1, 1)
        label = QtGui.QLabel(self.tab_synthesis)
        label.setText("Model observed resolution by smoothing")
        grid_synthesis.addWidget(label, 1, 1, 1, 1)

        label = QtGui.QLabel(self.tab_synthesis)
        label.setText("Constrain smoothing to less than:")
        grid_synthesis.addWidget(label, 2, 1, 1, 1)
        self.edit_smoothing_bound = QtGui.QLineEdit(self.tab_synthesis)
        self.edit_smoothing_bound.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_smoothing_bound.setMaximumSize(QtCore.QSize(60, 16777215))
        self.edit_smoothing_bound.setValidator(
            QtGui.QDoubleValidator(0, 10, 1, self.edit_smoothing_bound))
        grid_synthesis.addWidget(self.edit_smoothing_bound, 2, 3, 1, 1)

        self.btn_specify_abundances = QtGui.QPushButton(self.tab_synthesis)
        self.btn_specify_abundances.setText("Specify explicit abundance table")
        grid_synthesis.addWidget(self.btn_specify_abundances, 3, 1, 1, 1)

        grid_synthesis.addItem(QtGui.QSpacerItem(40, 20, 
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum), 0, 2, 1, 1)
        vbox_synthesis.addLayout(grid_synthesis)
        vbox_synthesis.addItem(QtGui.QSpacerItem(20, 40, 
            QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding))
        self.tabWidget.addTab(self.tab_synthesis, "Synthesis options")


        # Buttons.
        self.verticalLayout.addWidget(self.tabWidget)
        hbox_btns = QtGui.QHBoxLayout()
        self.btn_save_as_default = QtGui.QPushButton(group_box)
        self.btn_save_as_default.setAutoDefault(False)
        self.btn_save_as_default.setDefault(False)
        self.btn_save_as_default.setText("Save as default options")
        hbox_btns.addWidget(self.btn_save_as_default)

        hbox_btns.addItem(QtGui.QSpacerItem(40, 20, 
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))

        self.btn_apply_to_all = QtGui.QPushButton(group_box)
        self.btn_apply_to_all.setDefault(True)
        self.btn_apply_to_all.setText("Apply options to all models")
        hbox_btns.addWidget(self.btn_apply_to_all)

        # Final layout placement.
        self.verticalLayout.addLayout(hbox_btns)
        self.tabWidget.raise_()
        vbox_rhs.addWidget(group_box)
        hbox_parent.addLayout(vbox_rhs)


        self.mpl_axis = self.mpl_figure.figure.add_subplot(111)
        spectrum = self.spectral_models[0].session.normalized_spectrum
        self.mpl_axis.plot(spectrum.dispersion, spectrum.flux, c="k",
            drawstyle="steps-mid")
        sigma = 1.0/np.sqrt(spectrum.ivar)
        style_utils.fill_between_steps(self.mpl_axis, spectrum.dispersion,
            spectrum.flux - sigma, spectrum.flux + sigma, facecolor="#CCCCCC",
            edgecolor="None", alpha=0.5)

        self._lines = {
            "transitions_center": self.mpl_axis.axvline(
                np.nan, c="#666666", linestyle=":", zorder=-1),
            "model_fit": self.mpl_axis.plot([], [], c="r")[0],
            "nearby_lines": [],
            "model_masks": [],
        }
        self.mpl_axis.set_xlim(spectrum.dispersion[0], spectrum.dispersion[-1])
        self.mpl_axis.set_ylim(0, 1.2)

        self.mpl_axis.set_xlabel(u"Wavelength (Å)")
        self.mpl_axis.set_ylabel("Normalized flux")
        self.mpl_axis.xaxis.set_major_locator(MaxNLocator(5))
        self.mpl_axis.yaxis.set_major_locator(MaxNLocator(5))
        self.mpl_axis.xaxis.get_major_formatter().set_useOffset(False)
        self.mpl_figure.draw()

        # Select the first entry.
        self.spectral_models[0].fit()
        self.table_view.selectRow(0)

        # TODO: AND CUSTOM TEXT REPRS TO SHOW UNITS

        # Connect signals.
        self.btn_save_as_default.clicked.connect(self.clicked_save_as_default)
        self.btn_apply_to_all.clicked.connect(self.clicked_apply_to_all)

        # Common options.
        self.edit_window.textChanged.connect(self.update_edit_window)
        self.checkbox_continuum.stateChanged.connect(
            self.clicked_checkbox_continuum)
        self.combo_continuum.currentIndexChanged.connect(
            self.update_continuum_order)
        self.checkbox_vrad_tolerance.stateChanged.connect(
            self.clicked_checkbox_vrad_tolerance)
        self.edit_vrad_tolerance.textChanged.connect(
            self.update_vrad_tolerance)

        # Profile options.
        self.combo_profile.currentIndexChanged.connect(
            self.update_combo_profile)
        self.edit_detection_sigma.textChanged.connect(
            self.update_detection_sigma)
        self.edit_detection_pixels.textChanged.connect(
            self.update_detection_pixels)
        self.checkbox_use_central_weighting.stateChanged.connect(
            self.clicked_checkbox_use_central_weighting)
        self.checkbox_wavelength_tolerance.stateChanged.connect(
            self.clicked_checkbox_wavelength_tolerance)
        self.edit_wavelength_tolerance.textChanged.connect(
            self.update_wavelength_tolerance)

        # Synthesis options.
        self.edit_initial_abundance_bound.textChanged.connect(
            self.update_initial_abundance_bound)
        self.checkbox_model_smoothing.stateChanged.connect(
            self.clicked_checkbox_model_smoothing)
        self.edit_smoothing_bound.textChanged.connect(
            self.update_smoothing_bound)
        self.btn_specify_abundances.clicked.connect(
            self.clicked_btn_specify_abundances)

        # Connect matplotlib.
        self.mpl_figure.mpl_connect(
            "button_press_event", self.figure_mouse_press)
        self.mpl_figure.mpl_connect(
            "button_release_event", self.figure_mouse_release)

        return None


    def clicked_save_as_default(self):
        """ The 'Save as default' button has been clicked. """
        raise NotImplementedError


    def clicked_apply_to_all(self):
        """ The 'Apply to all models' button has been clicked. """
        raise NotImplementedError


    def update_edit_window(self):
        """ The window value has been updated. """

        model = self._get_selected_model()
        try:
            window = float(self.edit_window.text())

        except:
            return None

        else:
            model.metadata["window"] = window

            # Just update the axis limits.
            model = self._get_selected_model()
            transitions = model.transitions

            self.mpl_axis.set_xlim(
                transitions["wavelength"][0] - window,
                transitions["wavelength"][-1] + window,
            )
            self.mpl_figure.draw()

        return None


    def clicked_checkbox_continuum(self):
        """ The checkbox for modeling the continuum was clicked. """

        if self.checkbox_continuum.isChecked():
            self.combo_continuum.setEnabled(True)
            self.update_continuum_order()
        else:
            self._get_selected_model().metadata["continuum_order"] = -1
            self.combo_continuum.setEnabled(False)

        return None


    def update_continuum_order(self):
        """ The continuum order to use in model fitting was changed. """

        self._get_selected_model().metadata["continuum_order"] \
            = int(self.combo_continuum.currentText())
        return None


    def clicked_checkbox_vrad_tolerance(self):
        """ The checkbox for velocity tolerance was clicked. """

        if self.checkbox_vrad_tolerance.isChecked():
            self.edit_vrad_tolerance.setEnabled(True)
            self.update_vrad_tolerance()
        else:
            self.edit_vrad_tolerance.setEnabled(False)
            self._get_selected_model().metadata["velocity_tolerance"] = None
        return None


    def update_vrad_tolerance(self):
        """ The tolerance on radial velocity was updated. """
        try:
            value = float(self.edit_vrad_tolerance.text())
        except:
            value = None
        self._get_selected_model().metadata["velocity_tolerance"] = value
        return None


    def update_combo_profile(self):
        """ Update the profile that is used for fitting atomic transitions. """
        self._get_selected_model().metadata["profile"] \
            = self.update_combo_profile.currentText().lower()
        return None


    def update_detection_sigma(self):
        """ The detection sigma for nearby lines has been updated. """
        self._get_selected_model().metadata["detection_sigma"] \
            = float(self.edit_detection_sigma.text())
        return None


    def update_detection_pixels(self):
        """ The number of pixels to qualify a detection has been updated. """
        self._get_selected_model().metadata["detection_pixels"] \
            = int(self.edit_detection_pixels.text())
        return None


    def clicked_checkbox_use_central_weighting(self):
        """ The checkbox to use central weighting has been clicked. """

        self._get_selected_model().metadata["central_weighting"] \
            = self.checkbox_use_central_weighting.isChecked()
        return None


    def clicked_checkbox_wavelength_tolerance(self):
        """ The checkbox to set a wavelength tolerance has been clicked. """
        if self.checkbox_wavelength_tolerance.isChecked():
            self.edit_wavelength_tolerance.setEnabled(True)
            self.update_wavelength_tolerance()
        else:
            self.edit_wavelength_tolerance.setEnabled(False)
            self._get_selected_model().metadata["wavelength_tolerance"] = None
        return None


    def update_wavelength_tolerance(self):
        """ The wavelength tolerance for a profile centroid has been updated. """

        self._get_selected_model().metadata["wavelength_tolerance"] \
            = float(self.edit_wavelength_tolerance.text())
        return None


    def update_initial_abundance_bound(self):
        """ The initial abundance bound has been updated. """
        self._get_selected_model().metadata["initial_abundance_bounds"] \
            = float(self.edit_initial_abundance_bound.text())
        return None


    def clicked_checkbox_model_smoothing(self):
        """ The checkbox to smooth the model spectrum has been clicked. """

        if self.checkbox_model_smoothing.isChecked():
            self._get_selected_model().metadata["smoothing_kernel"] = True
            self.edit_smoothing_bound.setEnabled(True)
            self.update_smoothing_bound()

        else:
            self._get_selected_model().metadata["smoothing_kernel"] = False
            self.edit_smoothing_bound.setEnabled(False)
        return None


    def update_smoothing_bound(self):
        """ The limits on the smoothing kernel have been updated. """
        value = float(self.edit_smoothing_bound.text())
        self._get_selected_model().metadata["sigma_smooth"] = (-value, value)
        return None


    def clicked_btn_specify_abundances(self):
        """ Button to specify abundances for a synthesis setup was clicked. """

        raise NotImplementedError


    def _get_selected_model(self, full_output=False):
        index = self.table_view.selectionModel().selectedRows()[0].row()
        model = self.spectral_models[index]
        return (model, index) if full_output else model


    def figure_mouse_press(self, event):
        """
        Mouse button was pressed on the matplotlib widget.

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
                    self.redraw_figure()
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
                self.redraw_figure()
                return None

        else:
            # Single click.
            xmin, xmax, ymin, ymax = (event.xdata, np.nan, -1e8, +1e8)
            try:
                self._interactive_mask_region

            except AttributeError:
                self._interactive_mask_region = self.mpl_axis.axvspan(**{
                    "xmin": xmin,
                    "xmax": xmax,
                    "ymin": ymin,
                    "ymax": ymax,
                    "facecolor": "r",
                    "edgecolor": "none", # go home matplotlib, you're drunk
                    "alpha": 0.25,
                    "zorder": -1
                })

            else:
                self._interactive_mask_region.set_xy([
                    [xmin, ymin],
                    [xmin, ymax],
                    [xmax, ymax],
                    [xmax, ymin],
                    [xmin, ymin]
                ])

            # Set the signal and the time.
            self._interactive_mask_region_signal = (
                time(),
                self.mpl_figure.mpl_connect(
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

            data = self._interactive_mask_region.get_xy()

            # Update xmax.
            data[2:4, 0] = event.xdata
            self._interactive_mask_region.set_xy(data)
            self.mpl_figure.draw()

        return None


    def figure_mouse_release(self, event):
        """
        Mouse button was released from the matplotlib widget.

        :param event:
            The matplotlib event.
        """

        try:
            signal_time, signal_cid = self._interactive_mask_region_signal

        except AttributeError:
            return None

        xy = self._interactive_mask_region.get_xy()

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
            self.redraw_figure()

        xy[:, 0] = np.nan
        self._interactive_mask_region.set_xy(xy)

        self.mpl_figure.mpl_disconnect(signal_cid)
        self.mpl_figure.draw()
        del self._interactive_mask_region_signal
        return None


    def selected_model_changed(self):
        """
        Populate the widgets because the selected spectral model has changed.
        """

        try:
            selected_model = self._get_selected_model()
        except IndexError:
            return None

        # Common model.
        self.edit_window.setText("{}".format(selected_model.metadata["window"]))

        # Continuum order.
        continuum_order = selected_model.metadata["continuum_order"]
        if continuum_order < 0:
            self.checkbox_continuum.setChecked(False)
            self.combo_continuum.setEnabled(False)
        else:
            self.checkbox_continuum.setChecked(True)
            self.combo_continuum.setEnabled(True)
            self.combo_continuum.setCurrentIndex(continuum_order)

        # Radial velocity tolerance.
        vrad_tolerance = selected_model.metadata.get("velocity_tolerance", None)
        if vrad_tolerance is None:
            self.checkbox_vrad_tolerance.setChecked(False)
            self.edit_vrad_tolerance.setEnabled(False)
        else:
            self.checkbox_vrad_tolerance.setChecked(True)
            self.edit_vrad_tolerance.setEnabled(True)
            self.edit_vrad_tolerance.setText("{}".format(vrad_tolerance))

        # Profile options.
        if isinstance(selected_model, ProfileFittingModel):
            self.tab_profile.setEnabled(True)

            self.combo_profile.setCurrentIndex(
                ["gaussian", "lorentzian", "voight"].index(
                    selected_model.metadata["profile"]))

            self.edit_detection_sigma.setText("{}".format(
                selected_model.metadata["detection_sigma"]))
            self.edit_detection_pixels.setText("{}".format(
                selected_model.metadata["detection_pixels"]))

            self.checkbox_use_central_weighting.setEnabled(
                selected_model.metadata["central_weighting"])

            tolerance = selected_model.metadata.get("wavelength_tolerance", None)
            if tolerance is None:
                self.checkbox_wavelength_tolerance.setEnabled(False)
            else:
                self.checkbox_wavelength_tolerance.setEnabled(True)
                self.edit_wavelength_tolerance.setText(
                    "{}".format(tolerance))

        else:
            self.tab_profile.setEnabled(False)

        # Synthesis options.
        if isinstance(selected_model, SpectralSynthesisModel):
            self.tab_synthesis.setEnabled(True)

            # Update widgets.
            self.edit_initial_abundance_bound.setText(
                "{}".format(selected_model.metadata["initial_abundance_bounds"]))
            
            self.checkbox_model_smoothing.setEnabled(
                selected_model.metadata["smoothing_kernel"])

            # TODO sigma smooth tolerance needs implementing.
        else:
            self.tab_synthesis.setEnabled(False)


        # Update figure view.
        self.redraw_figure()


    def redraw_figure(self):
        """ Re-draw the matplotlib window based on a triggered changed. """

        # Set the axes.
        selected_model = self._get_selected_model()
        transitions = selected_model.transitions

        window = selected_model.metadata["window"]
        limits = [
            transitions["wavelength"][0] - window,
            transitions["wavelength"][-1] + window,
        ]

        # Zoom to region.
        self.mpl_axis.set_xlim(limits)

        if isinstance(selected_model, ProfileFittingModel):
            self._lines["transitions_center"].set_data(
                [transitions["wavelength"][0], transitions["wavelength"][0]],
                [0, 1.2])


        # Model masks specified by the user.
        # (These should be shown regardless of whether there is a fit or not.)
        for i, (start, end) in enumerate(selected_model.metadata["mask"]):
            try:
                patch = self._lines["model_masks"][i]
            except IndexError:
                self._lines["model_masks"].append(self.mpl_axis.axvspan(
                    np.nan, np.nan, facecolor="r", edgecolor="none",
                    alpha=0.25))
                patch = self._lines["model_masks"][-1]

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
        for unused_patch in self._lines["model_masks"][N:]:
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

        except KeyError:
            meta = {}
            self._lines["model_fit"].set_data([], [])
           
        else:
            self._lines["model_fit"].set_data(
                meta["model_x"], meta["model_y"])

            # Model yerr.
            if np.any(np.isfinite(meta["model_yerr"])):
                self._lines["model_yerr"] = self.mpl_axis.fill_between(
                    meta["model_x"],
                    meta["model_y"] + meta["model_yerr"],
                    meta["model_y"] - meta["model_yerr"],
                    facecolor="r", edgecolor="none", alpha=0.5)


            # Model masks due to nearby lines.
            if "nearby_lines" in meta:
                for i, (_, (start, end)) in enumerate(meta["nearby_lines"]):
                    try:
                        patch = self._lines["nearby_lines"][i]
                
                    except IndexError:
                        self._lines["nearby_lines"].append(
                            self.mpl_axis.axvspan(np.nan, np.nan,
                                facecolor="b", edgecolor="none", alpha=0.25))
                        patch = self._lines["nearby_lines"][-1]

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
        for each in self._lines["nearby_lines"][N:]:
            each.set_visible(False)

        # Hide any other model spectra?
        self.mpl_figure.draw()

        return None


if __name__ == "__main__":

    import sys

    if sys.platform == "darwin":
            
        # See http://successfulsoftware.net/2013/10/23/fixing-qt-4-for-mac-os-x-10-9-mavericks/
        substitutes = [
            (".Lucida Grande UI", "Lucida Grande"),
            (".Helvetica Neue DeskInterface", "Helvetica Neue")
        ]
        for substitute in substitutes:
            QtGui.QFont.insertSubstitution(*substitute)


    from smh import linelists, Session, specutils, spectral_models

    session = Session([
        "/Users/arc/codes/smh/hd44007red_multi.fits",
        "/Users/arc/codes/smh/hd44007blue_multi.fits",
    ])
    session.normalized_spectrum = specutils.Spectrum1D.read("../smh/hd44007-rest.fits")

    session.metadata["line_list"] = linelists.LineList.read("/Users/arc/research/ges/linelist/vilnius.ew.fe")

    session.metadata["spectral_models"] = []
    for i in range(3):
        session.metadata["spectral_models"].append(
            spectral_models.ProfileFittingModel(
                session,
                session.metadata["line_list"]["hash"][[i]]))

    session.metadata["spectral_models"].append(
        spectral_models.SpectralSynthesisModel(
            session,
            session.metadata["line_list"]["hash"][[i+1]],
            "Fe"))

    print("opening")
    try:
        app = QtGui.QApplication(sys.argv)
    except RuntimeError:
        None

    window = SpectralModelsDialog(session)
    window.exec_()


