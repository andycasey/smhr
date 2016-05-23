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

DOUBLE_CLICK_INTERVAL = 0.1 # MAGIC HACK



class SpectralModelsTableModel(QtCore.QAbstractTableModel):

    header = [" ", u"Wavelength\n(Å)", "Element", u"Equivalent Width\n(mÅ)",
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
                result = self._parent.spectral_models[index.row()]._result[2]
                equivalent_width = result["equivalent_width"]
            except:
                equivalent_width = (np.nan, np.nan, np.nan)

            if not np.any(np.isfinite(equivalent_width)):
                value = ""
            else:
                value = "{0:.1f} +/- {1:.1f}".format(
                    1000 * equivalent_width[0],
                    1000 * np.max(np.abs(equivalent_width[1:])))

        elif column == 4:
            value = ""

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

        self.setGeometry(800, 500, 800, 500)
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
        self.table_view.setMinimumSize(QtCore.QSize(0, 0))
        self.table_view.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.table_view.setModel(SpectralModelsTableModel(self))
        self.table_view.setSelectionBehavior(
            QtGui.QAbstractItemView.SelectRows)
        self.table_view.setSortingEnabled(True)
        self.table_view.resizeColumnsToContents()

        hbox_parent.addWidget(self.table_view)


        vbox_rhs = QtGui.QVBoxLayout()

        self.mpl_figure = mpl.MPLWidget(None, tight_layout=True)
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

        self.checkbox_vrad_tolerance = QtGui.QCheckBox(self.tab_common)
        self.checkbox_vrad_tolerance.setText("")
        grid_common.addWidget(self.checkbox_vrad_tolerance, 2, 0, 1, 1)
        label = QtGui.QLabel(self.tab_common)
        label.setText("Set tolerance on residual radial velocity")
        grid_common.addWidget(label, 2, 1, 1, 1)
        self.edit_vrad_tolerance = QtGui.QLineEdit(self.tab_common)
        self.edit_vrad_tolerance.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_vrad_tolerance.setMaximumSize(QtCore.QSize(60, 16777215))
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


        label = QtGui.QLabel(self.tab_profile)
        label.setText("Detection sigma for nearby absorption lines")
        grid_profile.addWidget(label, 1, 1, 1, 1)
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(40, 20, 
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.edit_detection_sigma = QtGui.QLineEdit(self.tab_profile)
        self.edit_detection_sigma.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_detection_sigma.setMaximumSize(QtCore.QSize(60, 16777215))
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
        grid_synthesis.addWidget(self.edit_initial_abundance_bound, 0, 3, 1, 1)

        self.check_model_smoothing = QtGui.QCheckBox(self.tab_synthesis)
        self.check_model_smoothing.setText("")
        grid_synthesis.addWidget(self.check_model_smoothing, 1, 0, 1, 1)
        label = QtGui.QLabel(self.tab_synthesis)
        label.setText("Model observed resolution by smoothing")
        grid_synthesis.addWidget(label, 1, 1, 1, 1)

        label = QtGui.QLabel(self.tab_synthesis)
        label.setText("Constrain smoothing to less than:")
        grid_synthesis.addWidget(label, 2, 1, 1, 1)
        self.edit_smoothing_bound = QtGui.QLineEdit(self.tab_synthesis)
        self.edit_smoothing_bound.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_smoothing_bound.setMaximumSize(QtCore.QSize(60, 16777215))
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
        self.mpl_axis.set_xlim(spectrum.dispersion[0], spectrum.dispersion[-1])
        self.mpl_axis.set_ylim(0, 1.2)

        self.mpl_axis.set_xlabel(u"Wavelength (Å)")
        self.mpl_axis.set_ylabel("Normalized flux")
        self.mpl_axis.xaxis.set_major_locator(MaxNLocator(5))
        self.mpl_axis.yaxis.set_major_locator(MaxNLocator(5))
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
    session.normalized_spectrum = specutils.Spectrum1D.read("test-spectrum.txt")
    #    "../smh/hd44007-rest.fits")

    session.metadata["line_list"] = linelists.LineList.read("/Users/arc/research/ges/linelist/vilnius.ew.fe")

    session.metadata["spectral_models"] = []
    for i in range(3):
        session.metadata["spectral_models"].append(
            spectral_models.ProfileFittingModel(
                session,
                session.metadata["line_list"]["hash"][[i]]))

    print("opening")
    try:
        app = QtGui.QApplication(sys.argv)
    except RuntimeError:
        None

    window = SpectralModelsDialog(session)
    window.exec_()


