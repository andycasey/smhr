#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The stellar parameters tab in Spectroscopy Made Hard """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["ChemicalAbundancesTab"]

import logging
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
from spectral_models_table import SpectralModelsTableView, SpectralModelsFilterProxyModel, SpectralModelsTableModel
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
        
        ################
        # LEFT HAND SIDE
        ################
        lhs_widget = QtGui.QWidget(self)
        lhs_layout = QtGui.QVBoxLayout()
        
        self.elem_combo_box = QtGui.QComboBox(self)
        lhs_layout.addWidget(self.elem_combo_box)
        # TODO add all species needed here whenever transitions are edited
        # TODO create different filters for each species
        # TODO summarize

        header = ["", u"λ\n(Å)", "Element\n", u"E. W.\n(mÅ)",
                  "log ε\n(dex)"]
        attrs = ("is_acceptable", "_repr_wavelength", "_repr_element", 
                 "equivalent_width", "abundance")
        self.table_view = SpectralModelsTableView(self)
        lhs_layout.addWidget(self.table_view)
        # Set up a proxymodel.
        self.proxy_spectral_models = SpectralModelsFilterProxyModel(self)
        self.proxy_spectral_models.add_filter_function(
            "use_for_stellar_composition_inference",
            lambda model: model.use_for_stellar_composition_inference)

        self.proxy_spectral_models.setDynamicSortFilter(True)
        self.proxy_spectral_models.setSourceModel(SpectralModelsTableModel(self, header, attrs))

        self.table_view.setModel(self.proxy_spectral_models)
        self.table_view.setSelectionBehavior(
            QtGui.QAbstractItemView.SelectRows)

        # TODO: Re-enable sorting.
        self.table_view.setSortingEnabled(False)
        self.table_view.resizeColumnsToContents()
        self.table_view.setColumnWidth(0, 30) # MAGIC
        self.table_view.setColumnWidth(1, 70) # MAGIC
        self.table_view.setColumnWidth(2, 70) # MAGIC
        self.table_view.setColumnWidth(3, 70) # MAGIC
        self.table_view.setMinimumSize(QtCore.QSize(240, 0))
        self.table_view.horizontalHeader().setStretchLastSection(True)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, 
            QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.table_view.sizePolicy().hasHeightForWidth())
        self.table_view.setSizePolicy(sp)
        lhs_layout.addWidget(self.table_view)
        # Abund tree
        #self.abundtree = AbundTreeView(self)
        #self.abundtree.setModel(AbundTreeModel(self))
        #self.abundtree.span_cols()
        #sp = QtGui.QSizePolicy(
        #    QtGui.QSizePolicy.MinimumExpanding, 
        #    QtGui.QSizePolicy.MinimumExpanding)
        #sp.setHeightForWidth(self.abundtree.sizePolicy().hasHeightForWidth())
        #self.abundtree.setSizePolicy(sp)
        #lhs_layout.addWidget(self.abundtree)

        # Buttons
        hbox = QtGui.QHBoxLayout()
        self.btn_fit_all = QtGui.QPushButton(self)
        self.btn_fit_all.setText("Fit and Measure all")
        self.btn_fit_one = QtGui.QPushButton(self)
        self.btn_fit_one.setText("Fit and Measure One")
        hbox.addWidget(self.btn_fit_all)
        hbox.addWidget(self.btn_fit_one)
        lhs_layout.addLayout(hbox)

        hbox = QtGui.QHBoxLayout()
        self.btn_refresh = QtGui.QPushButton(self)
        self.btn_refresh.setText("Refresh from session")
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
        
        # Model fitting options
        self._create_fitting_options_widget()
        lhs_layout.addWidget(self.fitting_options)

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
        #self.ax_residual.yaxis.set_major_locator(MaxNLocator(2))
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
        _ = self.table_view.selectionModel()
        _.selectionChanged.connect(self.selected_model_changed)

        # Connect buttons
        self.btn_fit_all.clicked.connect(self.fit_all)
        self.btn_fit_one.clicked.connect(self.fit_one)
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

    def _create_fitting_options_widget(self):
        group_box = QtGui.QGroupBox(self)
        group_box.setTitle("Fitting options")
        opt_layout = QtGui.QVBoxLayout(group_box)
        opt_layout.setContentsMargins(6, 12, 6, 6)
        self.opt_tabs = QtGui.QTabWidget(group_box)
        self.opt_tab_common = QtGui.QWidget()
        
        # Common options
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
        self.opt_tabs.addTab(self.tab_common, "Common")
        
        # Profile model options.
        self.tab_profile = QtGui.QWidget()
        grid_profile = QtGui.QGridLayout(self.tab_profile)


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

        self.opt_tabs.addTab(self.tab_profile, "Profile options")
        
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
        self.btn_specify_abundances.setText("Specify explicit abundance table TODO")
        grid_synthesis.addWidget(self.btn_specify_abundances, 3, 1, 1, 1)

        grid_synthesis.addItem(QtGui.QSpacerItem(40, 20, 
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum), 0, 2, 1, 1)
        vbox_synthesis.addLayout(grid_synthesis)
        vbox_synthesis.addItem(QtGui.QSpacerItem(20, 40, 
            QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding))
        self.opt_tabs.addTab(self.tab_synthesis, "Synthesis options")


        # Buttons.
        hbox_btns = QtGui.QHBoxLayout()
        self.auto_fit_checkbox = QtGui.QCheckBox(group_box)
        self.auto_fit_checkbox.setText("Autofit")
        hbox_btns.addWidget(self.auto_fit_checkbox)
        #hbox_btns = QtGui.QHBoxLayout()
        #self.btn_save_as_default = QtGui.QPushButton(group_box)
        #self.btn_save_as_default.setText("Save as default options")

        #hbox_btns.addItem(QtGui.QSpacerItem(40, 20, 
        #    QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))

        #hbox_btns.addWidget(self.btn_save_as_default)

        #self.btn_apply_to_all = QtGui.QPushButton(group_box)
        #self.btn_apply_to_all.setText("Apply options to all models")
        #hbox_btns.addWidget(self.btn_apply_to_all)

        # Final layout placement.
        opt_layout.addWidget(self.opt_tabs)
        opt_layout.addLayout(hbox_btns)
        self.opt_tabs.raise_()
        self.fitting_options = group_box

        # Connect Signals
        # Common options.
        self.edit_window.textChanged.connect(
            self.update_edit_window)
        self.edit_window.textChanged.connect(
            self.autofit)
        self.checkbox_continuum.stateChanged.connect(
            self.clicked_checkbox_continuum)
        self.checkbox_continuum.stateChanged.connect(
            self.autofit)
        self.combo_continuum.currentIndexChanged.connect(
            self.update_continuum_order)
        self.combo_continuum.currentIndexChanged.connect(
            self.autofit)
        self.checkbox_vrad_tolerance.stateChanged.connect(
            self.clicked_checkbox_vrad_tolerance)
        self.checkbox_vrad_tolerance.stateChanged.connect(
            self.autofit)
        self.edit_vrad_tolerance.textChanged.connect(
            self.update_vrad_tolerance)
        self.edit_vrad_tolerance.textChanged.connect(
            self.autofit)

        # Profile options.
        self.combo_profile.currentIndexChanged.connect(
            self.update_combo_profile)
        self.combo_profile.currentIndexChanged.connect(
            self.autofit)
        self.edit_detection_sigma.textChanged.connect(
            self.update_detection_sigma)
        self.edit_detection_sigma.textChanged.connect(
            self.autofit)
        self.edit_detection_pixels.textChanged.connect(
            self.update_detection_pixels)
        self.edit_detection_pixels.textChanged.connect(
            self.autofit)
        self.checkbox_use_central_weighting.stateChanged.connect(
            self.clicked_checkbox_use_central_weighting)
        self.checkbox_use_central_weighting.stateChanged.connect(
            self.autofit)
        self.checkbox_wavelength_tolerance.stateChanged.connect(
            self.clicked_checkbox_wavelength_tolerance)
        self.checkbox_wavelength_tolerance.stateChanged.connect(
            self.autofit)
        self.edit_wavelength_tolerance.textChanged.connect(
            self.update_wavelength_tolerance)
        self.edit_wavelength_tolerance.textChanged.connect(
            self.autofit)

        # Synthesis options.
        self.edit_initial_abundance_bound.textChanged.connect(
            self.update_initial_abundance_bound)
        self.edit_initial_abundance_bound.textChanged.connect(
            self.autofit)
        self.checkbox_model_smoothing.stateChanged.connect(
            self.clicked_checkbox_model_smoothing)
        self.checkbox_model_smoothing.stateChanged.connect(
            self.autofit)
        self.edit_smoothing_bound.textChanged.connect(
            self.update_smoothing_bound)
        self.edit_smoothing_bound.textChanged.connect(
            self.autofit)
        self.btn_specify_abundances.clicked.connect(
            self.clicked_btn_specify_abundances)
        self.btn_specify_abundances.clicked.connect(
            self.autofit)

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
        logger.debug("Resetting tree model from session")
        #self.abundtree.model().reset()
        self.proxy_spectral_models.reset()
        return None

    def refresh_plots(self):
        self.update_spectrum_figure(False)
        self.update_line_strength_figure(True)
        return None

    def fit_all(self):
        self._check_for_spectral_models()
        self.updated_spectral_models()
        # TODO order by complexity
        # TODO figure out is_acceptable
        logger.debug("Looping through spectral models...")
        # FIT ALL
        profile_measurements = []
        synthesis_measurements = []
        for i,m in enumerate(self.spectral_models):
            if isinstance(m, SpectralSynthesisModel):
                synthesis_measurements.append((i,m))
                try:
                    res = m.fit()
                    logger.debug(res)
                except (ValueError, RuntimeError) as e:
                    logger.debug("Fitting error",m)
                    logger.debug(e)
            if isinstance(m, ProfileFittingModel):
                profile_measurements.append((i,m))
                try:
                    res = m.fit()
                    logger.debug(res)
                except (ValueError, RuntimeError) as e:
                    logger.debug("Fitting error",m)
                    logger.debug(e)
        # CALL CURVE OF GROWTH FOR EW
        # TODO hook this into the session rather than the tab
        equivalent_widths = []
        transition_indices = []
        spectral_model_indices = []
        for i,m in profile_measurements:
            spectral_model_indices.append(i)
            transition_indices.extend(m._transition_indices)
            if m.is_acceptable:
                equivalent_widths.append(1000.* \
                    m.metadata["fitted_result"][-1]["equivalent_width"][0])
            else:
                equivalent_widths.append(np.nan)
        if len(equivalent_widths) == 0 \
        or np.isfinite(equivalent_widths).sum() == 0:
            raise ValueError("no measured transitions to calculate abundances")
        
        transition_indices = np.array(transition_indices)
        spectral_model_indices = np.array(spectral_model_indices)
        transitions = self.parent.session.metadata["line_list"][transition_indices].copy()
        transitions["equivalent_width"] = equivalent_widths
        finite = np.isfinite(transitions["equivalent_width"])
            
        # TODO right now it is not hooking into the session, be careful!
        abundances = self.parent.session.rt.abundance_cog(
            self.parent.session.stellar_photosphere, transitions[finite])
        for index, abundance in zip(spectral_model_indices[finite], abundances):
            self.spectral_models[index]\
                .metadata["fitted_result"][-1]["abundances"] = [abundance]

        #self.abundtree.model().reset()
        self.proxy_spectral_models.reset()
        self.selected_model_changed()

    def fit_one(self):
        spectral_model, proxy_index, index = self._get_selected_model(True)
        if spectral_model is None: return None
        try:
            res = spectral_model.fit()
            logger.debug(res)
        except (ValueError, RuntimeError) as e:
            logger.debug("Fitting error",spectral_model)
            logger.debug(e)
            return None
        try:
            ab = spectral_model.abundances
            logger.debug(ab)
        except rt.RTError as e:
            logger.debug("Abundance error",spectral_model)
            logger.debug(e)
            return None
        self.update_tree_data(index)
        self.selected_model_changed()
        return None

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
        #model = self.abundtree.model()
        model = self.proxy_spectral_models
        """
        ii = model.all_species == self.currently_plotted_species
        assert np.sum(ii) == 1, "{} {}".format(self.currently_plotted_species, model.all_species)
        summary_index = np.where(ii)[0]
        summary = model.summaries[summary_index]
        item = summary.subnodes[event.ind[0]]
        index = model.createIndex(event.ind[0],0,item)
        self.abundtree.setCurrentIndex(index)
        self.selected_model_changed()
        """
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
            spectral_model, proxy_index, index = self._get_selected_model(True)
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

                    # Update the view for this row.
                    self.update_tree_data(index)

                    # Update the view of the current model.
                    self.update_spectrum_figure(True)
                    break

            else:
                # No match with a masked region. 
                # TODO: Add a point that will be used for the continuum?
                # For the moment just refit the model.
                spectral_model.fit()

                # Update the view for this row.
                self.update_tree_data(index)

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
            spectral_model, proxy_index, index = self._get_selected_model(True)
            if spectral_model is None: 
                raise RuntimeError("""Must have a spectral model selected while making mask!
                                   Must have mouseover bug?""")

            # Add mask metadata.
            spectral_model.metadata["mask"].append([xy[0,0], xy[2, 0]])

            # Re-fit the spectral model.
            spectral_model.fit()

            # Update the table view for this row.
            self.update_tree_data(index)

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
        """
        index = self.abundtree.selectionModel().currentIndex()
        item = index.internalPointer()
        if isinstance(item,AbundTreeMeasurementItem):
            model = self.spectral_models[item.sm_ix]
            return (model, index) if full_output else model
        else:
            return (None, None) if full_output else None
        """
        proxy_index = self.table_view.selectionModel().selectedIndexes()[0]
        index = self.proxy_spectral_models.mapToSource(proxy_index).row()
        model = self.parent.session.metadata["spectral_models"][index]
        return (model, proxy_index, index) if full_output else model

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
        
        self.update_fitting_options()
        self.update_spectrum_figure(False)
        self.update_line_strength_figure(True)

        return None

    def update_spectrum_figure(self, redraw=False):
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
            three_sigma = 3*np.median(sigma[np.isfinite(sigma)])
            self.ax_residual.set_ylim(-three_sigma, three_sigma)

            if redraw: self.figure.draw()
        
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

        if redraw: self.figure.draw()

        return None

    def update_line_strength_figure(self, redraw=False):
        selected_model, proxy_index, index = self._get_selected_model(True)
        if selected_model is None:
            # TODO clear plot?
            return None
        """
        item = index.internalPointer() # AbundTreeMeasurementItem
        summary = item.parent

        # If the species is already plotted, don't replot
        assert isinstance(selected_model, ProfileFittingModel)
        ## TODO this doesn't update the plot when selecting/deselecting
        #if selected_model.transitions["species"][0] == self.currently_plotted_species:
        #    if redraw: self.figure.draw()
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
        
        if redraw: self.figure.draw()
        """
        return None

    def update_fitting_options(self):
        try:
            selected_model = self._get_selected_model()
        except IndexError:
            return None
        if selected_model is None: return None
    
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

        return None

    def update_tree_data(self, index):
        self.proxy_spectral_models.reset()
        """
        item = index.internalPointer()
        if not isinstance(item, AbundTreeMeasurementItem):
            raise RuntimeError(item)
        tree_model = self.abundtree.model()
        tree_model.dataChanged.emit(
            tree_model.createIndex(0, 0, item),
            tree_model.createIndex(0, item.columnCount(), item))
        """
        return None

    ###############################
    # FITTING OPTION UPDATE METHODS
    ###############################

    def update_edit_window(self):
        """ The wavelength window was updated """
        model = self._get_selected_model()
        try:
            window = float(self.edit_window.text())
        except:
            return None
        else:
            model.metadata["window"] = window
            # Just update the axis limits.
            transitions = model.transitions
            # TODO synth wavelength
            xlim = (transitions["wavelength"][0] - window,
                    transitions["wavelength"][-1] + window)
            self.ax_spectrum.set_xlim(xlim)
            self.ax_line_strength.set_xlim(xlim)
            self.figure.draw()
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
            = self.combo_profile.currentText().lower()
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
        if self.auto_fit_checkbox.isChecked(): self._get_selected_model().fit()
        return None
    def clicked_btn_specify_abundances(self):
        raise NotImplementedError

    def autofit(self):
        if self.auto_fit_checkbox.isChecked():
            m, pix, ix = self._get_selected_model(True)
            m.fit()
            #self.update_tree_data(ix)
            self.update_spectrum_figure(True)
            #self.update_line_strength_figure(True)
