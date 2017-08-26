#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The chemical abundances tab in Spectroscopy Made Hard """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["ChemicalAbundancesTab"]

import logging
import matplotlib.gridspec
import numpy as np
import sys
from PySide import QtCore, QtGui
import time

import smh
from smh.gui.base import SMHSpecDisplay
from smh.gui.base import MeasurementTableModelBase, MeasurementTableModelProxy, MeasurementTableView
from smh.gui import base
from smh import utils
import mpl, style_utils
from matplotlib.ticker import MaxNLocator
from smh.photospheres.abundances import asplund_2009 as solar_composition
from smh.spectral_models import (ProfileFittingModel, SpectralSynthesisModel)

from spectral_models_table import SpectralModelsTableViewBase, SpectralModelsFilterProxyModel, SpectralModelsTableModelBase
from linelist_manager import TransitionsDialog

logger = logging.getLogger(__name__)
logger.addHandler(smh.handler)

if sys.platform == "darwin":
        
    # See http://successfulsoftware.net/2013/10/23/fixing-qt-4-for-mac-os-x-10-9-mavericks/
    substitutes = [
        (".Lucida Grande UI", "Lucida Grande"),
        (".Helvetica Neue DeskInterface", "Helvetica Neue")
    ]
    for substitute in substitutes:
        QtGui.QFont.insertSubstitution(*substitute)

_QFONT = QtGui.QFont("Helvetica Neue", 10)
_ROWHEIGHT = 20
DOUBLE_CLICK_INTERVAL = 0.1 # MAGIC HACK
PICKER_TOLERANCE = 10 # MAGIC HACK


class ChemicalAbundancesTab(QtGui.QWidget):
    
    def __init__(self, parent):
        super(ChemicalAbundancesTab, self).__init__(parent)

        self.parent = parent
        self.FeH = np.nan

        #############################################################
        #### Create Layout of main widgets
        #############################################################

        self.parent_layout = QtGui.QVBoxLayout(self)
        ################
        # TOP: Spectrum
        ################
        self.figure = SMHSpecDisplay(None, self.parent.session, enable_masks=True,
                                     get_selected_model=self._get_selected_model)
        self.ax_spectrum = self.figure.ax_spectrum
        self.ax_residual = self.figure.ax_residual
        self.parent_layout.addWidget(self.figure)
        self.figure.add_callback_after_fit(self.refresh_current_model)
        self.figure.add_callback_after_fit(self.summarize_current_table)
        
        ################
        # BOTTOM
        ################
        bot_layout = QtGui.QHBoxLayout()
        # Measurement List
        bot_lhs_layout = self._create_measurement_list()
        bot_layout.addLayout(bot_lhs_layout)
        # Model fitting options
        self._create_fitting_options_widget()
        bot_layout.addWidget(self.opt_tabs)
        
        self.parent_layout.addLayout(bot_layout)

        #############################################################
        #### Connect signals
        #############################################################
        # Not: signals for spectrum plot already connected

        # Connect filter combo box
        self.filter_combo_box.currentIndexChanged.connect(self.filter_combo_box_changed)
        
        # Connect selection model
        #_ = self.table_view.selectionModel()
        #_.selectionChanged.connect(self.selected_model_changed)
        _ = self.synth_abund_table.selectionModel()
        _.selectionChanged.connect(self.update_spectrum_figure)

        # Connect buttons
        self.btn_fit_all.clicked.connect(self.fit_all_profiles)
        self.btn_measure_all.clicked.connect(self.measure_all)

        # TODO 

        # Set up things as if a fresh session
        self._currently_plotted_element = "All"
        self.new_session_loaded()

    def init_tab(self):
        """
        Call this to update the state of the tab from the session.
        """
        

    def new_session_loaded(self):
        """
        Call this whenever you have a new session or a new normalized spectrum
        """
        session = self.parent.session
        if session is None: return None
        #logger.debug("LOADING NEW SESSION")
        self.figure.new_session(session)
        self.refresh_table()
        self.summarize_current_table()
        self.refresh_plots()
        self.update_fitting_options()
        return None

    def _create_measurement_list(self):
        bot_lhs_layout = QtGui.QVBoxLayout()
        hbox = QtGui.QHBoxLayout()
        self.filter_combo_box = QtGui.QComboBox(self)
        self.filter_combo_box.setSizeAdjustPolicy(QtGui.QComboBox.AdjustToContents)
        self.filter_combo_box.addItem("All")
        self.element_summary_text = QtGui.QLabel(self)
        self.element_summary_text.setText("Please load spectral models (or fit all)")
        sp = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, 
                               QtGui.QSizePolicy.Minimum)
        self.element_summary_text.setSizePolicy(sp)
        hbox.addWidget(self.filter_combo_box)
        hbox.addWidget(self.element_summary_text)
        bot_lhs_layout.addLayout(hbox)

        self.full_measurement_model = MeasurementTableModelBase(self, self.parent.session, 
                                                    ["is_acceptable",
                                                     "species","wavelength",
                                                     "abundances",
                                                     "equivalent_width","reduced_equivalent_width",
                                                     "equivalent_width_uncertainty",
                                                     "expot","loggf",
                                                     "is_upper_limit","user_flag"])
        self.measurement_model = MeasurementTableModelProxy(self)
        self.measurement_model.setSourceModel(self.full_measurement_model)
        vbox, measurement_view, btn_filter, btn_refresh = base.create_measurement_table_with_buttons(
            self, self.measurement_model, self.parent.session, 
            callbacks_after_menu=[self.refresh_current_model,
                                  self.summarize_current_table,
                                  self.refresh_plots],
            display_fitting_options=True)
        self.measurement_view = measurement_view
        self.measurement_model.add_view_to_update(self.measurement_view)
        _ = self.measurement_view.selectionModel()
        _.selectionChanged.connect(self.selected_model_changed)
        self.measurement_view.setSizePolicy(QtGui.QSizePolicy(
            QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.MinimumExpanding))
        self.btn_filter_acceptable = btn_filter
        self.btn_filter_acceptable.clicked.connect(self.refresh_plots)
        self.btn_refresh = btn_refresh
        self.btn_refresh.clicked.connect(self.new_session_loaded)
        bot_lhs_layout.addLayout(vbox)

        #self.table_view = SpectralModelsTableView(self)
        #sp = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, 
        #                       QtGui.QSizePolicy.MinimumExpanding)
        #self.table_view.setSizePolicy(sp)
        ## Set up a proxymodel.
        #self.proxy_spectral_models = SpectralModelsFilterProxyModel(self)
        #self.proxy_spectral_models.add_filter_function(
        #    "use_for_stellar_composition_inference",
        #    lambda model: model.use_for_stellar_composition_inference)
        #self.proxy_spectral_models.setDynamicSortFilter(True)
        #header = ["", u"λ", "log ε", u"E. W.",
        #          "REW", "σ(X)", "σ(E.W.)", "loggf","Element"]
        #self.all_spectral_models = SpectralModelsTableModel(self, header, None)
        #self.proxy_spectral_models.setSourceModel(self.all_spectral_models)
        #self.table_view.setModel(self.proxy_spectral_models)
        #self.table_view.setSelectionBehavior(
        #    QtGui.QAbstractItemView.SelectRows)

        # TODO: Re-enable sorting.
        #self.table_view.setSortingEnabled(False)
        #self.table_view.resizeColumnsToContents()
        #self.table_view.resizeRowsToContents()
        #self.table_view.verticalHeader().setResizeMode(QtGui.QHeaderView.Fixed)
        #self.table_view.verticalHeader().setDefaultSectionSize(_ROWHEIGHT)
        #self.table_view.setColumnWidth(0, 25) # MAGIC
        #self.table_view.setColumnWidth(1, 50) # MAGIC
        #self.table_view.setColumnWidth(2, 50) # MAGIC
        #self.table_view.setColumnWidth(3, 50) # MAGIC
        #self.table_view.setColumnWidth(4, 50) # MAGIC
        #self.table_view.setColumnWidth(5, 50) # MAGIC
        #self.table_view.setColumnWidth(6, 50) # MAGIC
        #self.table_view.setColumnWidth(7, 50) # MAGIC
        #self.table_view.setMinimumSize(QtCore.QSize(240, 0))
        #self.table_view.horizontalHeader().setStretchLastSection(True)
        #bot_lhs_layout.addWidget(self.table_view)

        # Buttons
        hbox = QtGui.QHBoxLayout()
        sp = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, 
                               QtGui.QSizePolicy.Minimum)
        self.btn_fit_all = QtGui.QPushButton(self)
        self.btn_fit_all.setText("Fit all EW")
        self.btn_fit_all.setSizePolicy(sp)
        self.btn_measure_all = QtGui.QPushButton(self)
        self.btn_measure_all.setText("Measure all acceptable EW")
        self.btn_measure_all.setSizePolicy(sp)
        hbox.addWidget(self.btn_fit_all)
        hbox.addWidget(self.btn_measure_all)
        bot_lhs_layout.addLayout(hbox)
        
        return bot_lhs_layout

    def _create_fitting_options_widget(self):
        self.opt_tabs = QtGui.QTabWidget(self)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Minimum,
            QtGui.QSizePolicy.MinimumExpanding)
        self.opt_tabs.setSizePolicy(sp)

        def _create_line_in_hbox(parent, text, bot, top, dec, validate_int=False):
            hbox = QtGui.QHBoxLayout()
            hbox.setSpacing(0)
            hbox.setContentsMargins(0,0,0,0)
            label = QtGui.QLabel(parent)
            label.setText(text)
            label.setFont(_QFONT)
            line = QtGui.QLineEdit(parent)
            line.setMinimumSize(QtCore.QSize(60, 0))
            line.setMaximumSize(QtCore.QSize(60, _ROWHEIGHT))
            line.setFont(_QFONT)
            if validate_int:
                line.setValidator(QtGui.QIntValidator(bot, top, line))
            else:
                line.setValidator(QtGui.QDoubleValidator(bot, top, dec, line))
            hbox.addWidget(label)
            hbox.addItem(QtGui.QSpacerItem(40, _ROWHEIGHT, QtGui.QSizePolicy.Expanding,
                                           QtGui.QSizePolicy.Minimum))
            hbox.addWidget(line)
            return hbox, label, line

        def _create_checkline_in_hbox(parent, text, bot, top, dec):
            hbox = QtGui.QHBoxLayout()
            hbox.setSpacing(0)
            hbox.setContentsMargins(0,0,0,0)
            checkbox = QtGui.QCheckBox(parent)
            label = QtGui.QLabel(parent)
            label.setText(text)
            label.setFont(_QFONT)
            line = QtGui.QLineEdit(parent)
            line.setMinimumSize(QtCore.QSize(60, 0))
            line.setMaximumSize(QtCore.QSize(60, _ROWHEIGHT))
            line.setFont(_QFONT)
            line.setValidator(QtGui.QDoubleValidator(bot, top, dec, line))
            hbox.addWidget(checkbox)
            hbox.addItem(QtGui.QSpacerItem(10, 10, QtGui.QSizePolicy.Fixed,
                                           QtGui.QSizePolicy.Minimum))
            hbox.addWidget(label)
            hbox.addItem(QtGui.QSpacerItem(20, _ROWHEIGHT, QtGui.QSizePolicy.Expanding,
                                           QtGui.QSizePolicy.Minimum))
            hbox.addWidget(line)
            return hbox, checkbox, label, line

        def _create_combo_in_hbox(parent, text):
            hbox = QtGui.QHBoxLayout()
            hbox.setSpacing(0)
            hbox.setContentsMargins(0,0,0,0)
            label = QtGui.QLabel(parent)
            label.setText(text)
            label.setFont(_QFONT)
            combo = QtGui.QComboBox(parent)
            combo.setFont(_QFONT)
            combo.setSizeAdjustPolicy(QtGui.QComboBox.AdjustToContents)
            combo.setMinimumSize(QtCore.QSize(60, 0))
            combo.setMaximumSize(QtCore.QSize(1000, _ROWHEIGHT))
            hbox.addWidget(label)
            hbox.addItem(QtGui.QSpacerItem(20, _ROWHEIGHT, QtGui.QSizePolicy.Expanding,
                                           QtGui.QSizePolicy.Minimum))
            hbox.addWidget(combo)
            return hbox, label, combo

        def _create_checkcombo_in_hbox(parent, text):
            hbox = QtGui.QHBoxLayout()
            hbox.setSpacing(0)
            hbox.setContentsMargins(0,0,0,0)
            checkbox = QtGui.QCheckBox(parent)
            label = QtGui.QLabel(parent)
            label.setText(text)
            label.setFont(_QFONT)
            combo = QtGui.QComboBox(parent)
            combo.setFont(_QFONT)
            combo.setMinimumSize(QtCore.QSize(60, 0))
            combo.setMaximumSize(QtCore.QSize(60, _ROWHEIGHT))
            hbox.addWidget(checkbox)
            hbox.addItem(QtGui.QSpacerItem(10, 10, QtGui.QSizePolicy.Fixed,
                                           QtGui.QSizePolicy.Minimum))
            hbox.addWidget(label)
            hbox.addItem(QtGui.QSpacerItem(20, _ROWHEIGHT, QtGui.QSizePolicy.Expanding,
                                           QtGui.QSizePolicy.Minimum))
            hbox.addWidget(combo)
            return hbox, checkbox, label, combo

        ###################
        ### Profile options
        self.tab_profile = QtGui.QWidget(self.opt_tabs)
        tab_hbox = QtGui.QHBoxLayout(self.tab_profile)
        tab_hbox.setSpacing(0)
        tab_hbox.setContentsMargins(0,0,0,0)

        ### LHS
        vbox_lhs = QtGui.QVBoxLayout()
        vbox_lhs.setSpacing(0)
        vbox_lhs.setContentsMargins(0,0,0,0)
        hbox, label, line = _create_line_in_hbox(self.tab_profile, "View window",
                                                 0, 1000, 1)
        self.edit_view_window = line
        vbox_lhs.addLayout(hbox)

        hbox, label, line = _create_line_in_hbox(self.tab_profile, "Fit window",
                                                 0, 1000, 1)
        self.edit_fit_window = line
        vbox_lhs.addLayout(hbox)

        hbox, checkbox, label, combo = _create_checkcombo_in_hbox(self.tab_profile, 
                                                                  "Use poly for cont.")
        self.checkbox_continuum = checkbox
        self.combo_continuum = combo
        for i in range(10):
            self.combo_continuum.addItem("{:.0f}".format(i))
        vbox_lhs.addLayout(hbox)

        hbox, checkbox, label, line = _create_checkline_in_hbox(self.tab_profile, "RV tol",
                                                                0, 100, 2)
        self.checkbox_vrad_tolerance = checkbox
        self.edit_vrad_tolerance = line
        vbox_lhs.addLayout(hbox)

        hbox, checkbox, label, line = _create_checkline_in_hbox(self.tab_profile, "WL tol",
                                                                0, 10, 2)
        self.checkbox_wavelength_tolerance = checkbox
        self.edit_wavelength_tolerance = line
        vbox_lhs.addLayout(hbox)

        self.checkbox_use_central_weighting = QtGui.QCheckBox(self.tab_profile)
        self.checkbox_use_central_weighting.setText("Central pixel weighting")
        self.checkbox_use_central_weighting.setFont(_QFONT)
        vbox_lhs.addWidget(self.checkbox_use_central_weighting)

        self.checkbox_use_antimasks = QtGui.QCheckBox(self.tab_profile)
        self.checkbox_use_antimasks.setText("Use Antimasks")
        self.checkbox_use_antimasks.setEnabled(False) # Editable by shift clicking only
        self.checkbox_use_antimasks.setFont(_QFONT)
        vbox_lhs.addWidget(self.checkbox_use_antimasks)

        ### RHS
        vbox_rhs = QtGui.QVBoxLayout()
        vbox_rhs.setSpacing(0)
        vbox_rhs.setContentsMargins(0,0,0,0)
        hbox, label, combo = _create_combo_in_hbox(self.tab_profile, "Type")
        self.combo_profile = combo
        for each in ("Gaussian", "Lorentzian", "Voigt"):
            self.combo_profile.addItem(each)
        vbox_rhs.addLayout(hbox)
        
        hbox, label, line = _create_line_in_hbox(self.tab_profile, "Automask sigma",
                                                 0, 100, 1)
        self.edit_detection_sigma = line
        vbox_rhs.addLayout(hbox)
        
        hbox, label, line = _create_line_in_hbox(self.tab_profile, "Automask pixels",
                                                 0, 100, np.nan, validate_int=True)
        self.edit_detection_pixels = line
        vbox_rhs.addLayout(hbox)
        
        vbox_rhs.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Minimum,
                                           QtGui.QSizePolicy.Expanding))
        
        self.checkbox_upper_limit = QtGui.QCheckBox(self.tab_profile)
        self.checkbox_upper_limit.setText("Upper Limit")
        self.checkbox_upper_limit.setFont(_QFONT)
        vbox_rhs.addWidget(self.checkbox_upper_limit)

        self.btn_fit_one = QtGui.QPushButton(self.tab_profile)
        self.btn_fit_one.setText("Fit One")
        vbox_rhs.addWidget(self.btn_fit_one)
        
        self.btn_clear_masks = QtGui.QPushButton(self.tab_profile)
        self.btn_clear_masks.setText("Clear Masks")
        vbox_rhs.addWidget(self.btn_clear_masks)
        
        ### Finish Profile Tab
        tab_hbox.addLayout(vbox_lhs)
        tab_hbox.addLayout(vbox_rhs)
        self.opt_tabs.addTab(self.tab_profile, "Profile")
        
        # Connect Signals for Profile: after synthesis

        ###################
        ### Synthesis options
        self.tab_synthesis = QtGui.QWidget(self.opt_tabs)
        tab_hbox = QtGui.QHBoxLayout(self.tab_synthesis)
        tab_hbox.setSpacing(0)
        tab_hbox.setContentsMargins(0,0,0,0)

        ### LHS
        vbox_lhs = QtGui.QVBoxLayout()
        vbox_lhs.setSpacing(0)
        vbox_lhs.setContentsMargins(0,0,0,0)
        hbox, label, line = _create_line_in_hbox(self.tab_synthesis, "View window",
                                                 0, 1000, 1)
        self.edit_view_window_2 = line
        vbox_lhs.addLayout(hbox)

        hbox, label, line = _create_line_in_hbox(self.tab_synthesis, "Fit window",
                                                 0, 1000, 1)
        self.edit_fit_window_2 = line
        vbox_lhs.addLayout(hbox)

        hbox, checkbox, label, combo = _create_checkcombo_in_hbox(self.tab_synthesis, 
                                                                  "Use poly for cont.")
        self.checkbox_continuum_2 = checkbox
        self.combo_continuum_2 = combo
        for i in range(10):
            self.combo_continuum_2.addItem("{:.0f}".format(i))
        vbox_lhs.addLayout(hbox)

        hbox, label, line = _create_line_in_hbox(self.tab_synthesis, "Manual cont.", 
                                                 -10, 10, 4)
        self.edit_manual_continuum = line
        vbox_lhs.addLayout(hbox)

        hbox, checkbox, label, line = _create_checkline_in_hbox(self.tab_synthesis, "RV tol",
                                                                0, 100, 2)
        self.checkbox_vrad_tolerance_2 = checkbox
        self.edit_vrad_tolerance_2 = line
        vbox_lhs.addLayout(hbox)
        
        hbox, label, line = _create_line_in_hbox(self.tab_synthesis, "Manual RV", 
                                                 -1000, 1000, 4)
        self.edit_manual_rv = line
        vbox_lhs.addLayout(hbox)

        hbox, checkbox, label, line = _create_checkline_in_hbox(self.tab_synthesis, "Smoothing",
                                                                0, 10, 3)
        self.checkbox_model_smoothing = checkbox
        self.edit_manual_smoothing = line
        vbox_lhs.addLayout(hbox)
        
        hbox, label, line = _create_line_in_hbox(self.tab_synthesis, "Initial abund bound",
                                                 0.01, 2, 2)
        self.edit_initial_abundance_bound = line
        vbox_lhs.addLayout(hbox)

        ### RHS
        vbox_rhs = QtGui.QVBoxLayout()
        vbox_rhs.setSpacing(0)
        vbox_rhs.setContentsMargins(0,0,0,0)

        # Element abundance table
        self.synth_abund_table = SynthesisAbundanceTableView(self.tab_synthesis, )
        self.synth_abund_table_model = SynthesisAbundanceTableModel(self)
        self.synth_abund_table.setModel(self.synth_abund_table_model)
        self.synth_abund_table.resizeColumnsToContents()
        self.synth_abund_table.verticalHeader().setResizeMode(QtGui.QHeaderView.Fixed)
        self.synth_abund_table.verticalHeader().setDefaultSectionSize(_ROWHEIGHT)
        self.synth_abund_table.setColumnWidth(0, 40) # MAGIC
        self.synth_abund_table.setColumnWidth(1, 55) # MAGIC
        self.synth_abund_table.horizontalHeader().setStretchLastSection(True)
        sp = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, 
                               QtGui.QSizePolicy.MinimumExpanding)
        self.synth_abund_table.setSizePolicy(sp)
        vbox_rhs.addWidget(self.synth_abund_table)
        
        self.btn_synthesize = QtGui.QPushButton(self.tab_synthesis)
        self.btn_synthesize.setText("Synthesize")
        vbox_rhs.addWidget(self.btn_synthesize)

        self.checkbox_upper_limit_2 = QtGui.QCheckBox(self.tab_synthesis)
        self.checkbox_upper_limit_2.setText("Upper Limit")
        self.checkbox_upper_limit_2.setFont(_QFONT)
        vbox_rhs.addWidget(self.checkbox_upper_limit_2)

        #vbox_rhs.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Minimum,
        #                                   QtGui.QSizePolicy.Minimum))

        self.btn_fit_synth = QtGui.QPushButton(self.tab_synthesis)
        self.btn_fit_synth.setText("Fit Model")
        vbox_rhs.addWidget(self.btn_fit_synth)
        
        self.btn_update_abund_table = QtGui.QPushButton(self.tab_synthesis)
        self.btn_update_abund_table.setText("Update Abundance Table")
        vbox_rhs.addWidget(self.btn_update_abund_table)
        
        
        hbox = QtGui.QHBoxLayout()
        self.btn_clear_masks_2 = QtGui.QPushButton(self.tab_synthesis)
        self.btn_clear_masks_2.setText("Clear Masks")
        hbox.addWidget(self.btn_clear_masks_2)
        self.btn_export_synth = QtGui.QPushButton(self.tab_synthesis)
        self.btn_export_synth.setText("Export")
        hbox.addWidget(self.btn_export_synth)
        vbox_rhs.addLayout(hbox)

        ### Finish Synthesis Tab
        tab_hbox.addLayout(vbox_lhs)
        tab_hbox.addLayout(vbox_rhs)
        self.opt_tabs.addTab(self.tab_synthesis, "Synthesis")

        ## opt_tabs settings
        ## TODO: for some reason, the window automatically expands ruthlessly unless I put this in...
        #self.opt_tabs.setMaximumSize(400,250)
        #self.opt_tabs.tabBar().setFont(_QFONT)
        # There's actually no need to show the tabs!
        self.opt_tabs.tabBar().setMaximumSize(0,0)

        # Connect signals for Profile and Synthesis
        self._connect_profile_signals()

    def _disconnect_profile_signals(self):
        for signal_obj, method in self._profile_signals:
            signal_obj.disconnect(method)
        for signal_obj, method in self._synth_signals:
            signal_obj.disconnect(method)
    def _connect_profile_signals(self):
        self.edit_view_window.textChanged.connect(
            self.update_edit_view_window)
        self.edit_fit_window.textChanged.connect(
            self.update_edit_fit_window)
        self.edit_fit_window.returnPressed.connect(
            self.fit_one)
        self.checkbox_continuum.stateChanged.connect(
            self.clicked_checkbox_continuum)
        self.checkbox_continuum.stateChanged.connect(
            self.fit_one)
        self.combo_continuum.currentIndexChanged.connect(
            self.update_continuum_order)
        self.combo_continuum.currentIndexChanged.connect(
            self.fit_one)
        self.checkbox_vrad_tolerance.stateChanged.connect(
            self.clicked_checkbox_vrad_tolerance)
        self.checkbox_vrad_tolerance.stateChanged.connect(
            self.fit_one)
        self.edit_vrad_tolerance.textChanged.connect(
            self.update_vrad_tolerance)
        self.edit_vrad_tolerance.returnPressed.connect(
            self.fit_one)
        self.combo_profile.currentIndexChanged.connect(
            self.update_combo_profile)
        self.combo_profile.currentIndexChanged.connect(
            self.fit_one)
        self.edit_detection_sigma.textChanged.connect(
            self.update_detection_sigma)
        self.edit_detection_sigma.returnPressed.connect(
            self.fit_one)
        self.edit_detection_pixels.textChanged.connect(
            self.update_detection_pixels)
        self.edit_detection_pixels.returnPressed.connect(
            self.fit_one)
        self.checkbox_use_central_weighting.stateChanged.connect(
            self.clicked_checkbox_use_central_weighting)
        self.checkbox_use_central_weighting.stateChanged.connect(
            self.fit_one)
        self.checkbox_use_antimasks.stateChanged.connect(
            self.clicked_checkbox_use_antimasks)
        self.checkbox_wavelength_tolerance.stateChanged.connect(
            self.clicked_checkbox_wavelength_tolerance)
        self.checkbox_wavelength_tolerance.stateChanged.connect(
            self.fit_one)
        self.edit_wavelength_tolerance.textChanged.connect(
            self.update_wavelength_tolerance)
        self.edit_wavelength_tolerance.returnPressed.connect(
            self.fit_one)
        self.checkbox_upper_limit.stateChanged.connect(
            self.clicked_checkbox_upper_limit)
        self.btn_fit_one.clicked.connect(
            self.fit_one)
        self.btn_clear_masks.clicked.connect(
            self.clicked_btn_clear_masks)
        self._profile_signals = [
            (self.edit_view_window.textChanged,self.update_edit_view_window),
            (self.edit_fit_window.textChanged,self.update_edit_fit_window),
            (self.edit_fit_window.returnPressed,self.fit_one),
            (self.checkbox_continuum.stateChanged,self.clicked_checkbox_continuum),
            (self.checkbox_continuum.stateChanged,self.fit_one),
            (self.combo_continuum.currentIndexChanged,self.update_continuum_order),
            (self.combo_continuum.currentIndexChanged,self.fit_one),
            (self.checkbox_vrad_tolerance.stateChanged,self.clicked_checkbox_vrad_tolerance),
            (self.checkbox_vrad_tolerance.stateChanged,self.fit_one),
            (self.edit_vrad_tolerance.textChanged,self.update_vrad_tolerance),
            (self.edit_vrad_tolerance.returnPressed,self.fit_one),
            (self.combo_profile.currentIndexChanged,self.update_combo_profile),
            (self.combo_profile.currentIndexChanged,self.fit_one),
            (self.edit_detection_sigma.textChanged,self.update_detection_sigma),
            (self.edit_detection_sigma.returnPressed,self.fit_one),
            (self.edit_detection_pixels.textChanged,self.update_detection_pixels),
            (self.edit_detection_pixels.returnPressed,self.fit_one),
            (self.checkbox_use_central_weighting.stateChanged,self.clicked_checkbox_use_central_weighting),
            (self.checkbox_use_central_weighting.stateChanged,self.fit_one),
            (self.checkbox_use_antimasks.stateChanged,self.clicked_checkbox_use_antimasks),
            (self.checkbox_wavelength_tolerance.stateChanged,self.clicked_checkbox_wavelength_tolerance),
            (self.checkbox_wavelength_tolerance.stateChanged,self.fit_one),
            (self.edit_wavelength_tolerance.textChanged,self.update_wavelength_tolerance),
            (self.edit_wavelength_tolerance.returnPressed,self.fit_one),
            (self.checkbox_upper_limit.stateChanged,self.clicked_checkbox_upper_limit),
            (self.btn_fit_one.clicked,self.fit_one),
            (self.btn_clear_masks.clicked,self.clicked_btn_clear_masks)
            ]

        #self.synthDefaultAction = self.fit_one

        self.edit_view_window_2.textChanged.connect(
            self.update_edit_view_window_2)
        self.edit_fit_window_2.textChanged.connect(
            self.update_edit_fit_window_2)
        self.checkbox_continuum_2.stateChanged.connect(
            self.clicked_checkbox_continuum_2)
        #self.checkbox_continuum_2.stateChanged.connect(
        #    self.synthDefaultAction)
        self.combo_continuum_2.currentIndexChanged.connect(
            self.update_continuum_order_2)
        #self.combo_continuum_2.currentIndexChanged.connect(
        #    self.synthDefaultAction)
        self.edit_manual_continuum.textChanged.connect(
            self.update_edit_manual_continuum)
        #self.edit_manual_continuum.returnPressed.connect(
        #    self.synthDefaultAction)
        self.checkbox_vrad_tolerance_2.stateChanged.connect(
            self.clicked_checkbox_vrad_tolerance_2)
        #self.checkbox_vrad_tolerance_2.stateChanged.connect(
        #    self.synthDefaultAction)
        self.edit_vrad_tolerance_2.textChanged.connect(
            self.update_vrad_tolerance_2)
        #self.edit_vrad_tolerance_2.returnPressed.connect(
        #    self.synthDefaultAction)
        self.edit_manual_rv.textChanged.connect(
            self.update_manual_rv)
        #self.edit_manual_rv.returnPressed.connect(
        #    self.synthDefaultAction)
        self.checkbox_model_smoothing.stateChanged.connect(
            self.clicked_checkbox_model_smoothing)
        #self.checkbox_model_smoothing.stateChanged.connect(
        #    self.synthDefaultAction)
        self.edit_manual_smoothing.textChanged.connect(
            self.update_manual_smoothing)
        #self.edit_manual_smoothing.returnPressed.connect(
        #    self.synthDefaultAction)
        self.edit_initial_abundance_bound.textChanged.connect(
            self.update_initial_abundance_bound)
        #self.edit_initial_abundance_bound.returnPressed.connect(
        #    self.synthDefaultAction)
        self.btn_synthesize.clicked.connect(
            self.synthesize_current_model)
        self.checkbox_upper_limit_2.stateChanged.connect(
            self.clicked_checkbox_upper_limit_2)
        self.btn_fit_synth.clicked.connect(
            self.fit_one)
        self.btn_update_abund_table.clicked.connect(
            self.clicked_btn_update_abund_table)
        self.btn_clear_masks_2.clicked.connect(
            self.clicked_btn_clear_masks)
        self.btn_export_synth.clicked.connect(
            self.clicked_export_synthesis)

        self._synth_signals = [
            (self.edit_view_window_2.textChanged,self.update_edit_view_window_2),
            (self.edit_fit_window_2.textChanged,self.update_edit_fit_window_2),
            (self.checkbox_continuum_2.stateChanged,self.clicked_checkbox_continuum_2),
            #(self.checkbox_continuum_2.stateChanged,self.synthDefaultAction),
            (self.combo_continuum_2.currentIndexChanged,self.update_continuum_order_2),
            #(self.combo_continuum_2.currentIndexChanged,self.synthDefaultAction),
            (self.edit_manual_continuum.textChanged,self.update_edit_manual_continuum),
            #(self.edit_manual_continuum.returnPressed,self.synthDefaultAction),
            (self.checkbox_vrad_tolerance_2.stateChanged,self.clicked_checkbox_vrad_tolerance_2),
            #(self.checkbox_vrad_tolerance_2.stateChanged,self.synthDefaultAction),
            (self.edit_vrad_tolerance_2.textChanged,self.update_vrad_tolerance_2),
            #(self.edit_vrad_tolerance_2.returnPressed,self.synthDefaultAction),
            (self.edit_manual_rv.textChanged,self.update_manual_rv),
            #(self.edit_manual_rv.returnPressed,self.synthDefaultAction),
            (self.checkbox_model_smoothing.stateChanged,self.clicked_checkbox_model_smoothing),
            #(self.checkbox_model_smoothing.stateChanged,self.synthDefaultAction),
            (self.edit_manual_smoothing.textChanged,self.update_manual_smoothing),
            #(self.edit_manual_smoothing.returnPressed,self.synthDefaultAction),
            (self.edit_initial_abundance_bound.textChanged,self.update_initial_abundance_bound),
            #(self.edit_initial_abundance_bound.returnPressed,self.synthDefaultAction),
            (self.btn_synthesize.clicked,self.synthesize_current_model),
            (self.checkbox_upper_limit_2.stateChanged,self.clicked_checkbox_upper_limit_2),
            (self.btn_fit_synth.clicked,self.fit_one),
            (self.btn_update_abund_table.clicked,self.clicked_btn_update_abund_table),
            (self.btn_clear_masks_2.clicked,self.clicked_btn_clear_masks),
            (self.btn_export_synth.clicked,self.clicked_export_synthesis)
            ]

    def populate_filter_combo_box(self):
        if self.parent.session is None: return None
        box = self.filter_combo_box
        box.clear()
        box.addItem("All")

        all_species = set([])
        for spectral_model in self.full_measurement_model.spectral_models:
            if isinstance(spectral_model, ProfileFittingModel):
                all_species.update(set(spectral_model.species))
            elif isinstance(spectral_model, SpectralSynthesisModel):
                for specie in spectral_model.species:
                    all_species.update(set(specie))
        if len(all_species)==0: return None
        all_species = np.sort(list(all_species))
        for species in all_species:
            elem = utils.species_to_element(species)
            assert species == utils.element_to_species(elem)
            box.addItem(elem)

    def filter_combo_box_changed(self):
        elem = self.filter_combo_box.currentText()
        # Update the filter
        if self._currently_plotted_element not in ["All", "", "None"]:
            try:
                self.measurement_model.delete_filter_function(self._currently_plotted_element)
            except KeyError as e:
                logger.debug(self._currently_plotted_element)
                logger.debug(e)
                logger.debug(self.measurement_model.filter_functions)
                raise            
        if elem in [None, "", "All"]:
            self.element_summary_text.setText("")
        else:
            species = utils.element_to_species(elem)
            def filter_function(model):
                if isinstance(model, ProfileFittingModel):
                    return species in model.species
                elif isinstance(model, SpectralSynthesisModel):
                    return np.any([species in specie for specie in model.species])
            self.measurement_model.add_filter_function(elem, filter_function)
        self._currently_plotted_element = elem
        self.measurement_model.reset()
        self.summarize_current_table()
        self.refresh_plots()
        self.measurement_view.selectRow(0)
        return None

    def calculate_FeH(self):
        summary_dict = self.parent.session.summarize_spectral_models()
        # TODO using Fe I right now
        try:
            self.FeH = summary_dict[26.0][4]
        except KeyError: # No Fe measured yet
            return np.nan
        return self.FeH

    def summarize_current_table(self):
        elem = self.filter_combo_box.currentText()
        if elem is None or elem == "" or elem == "All":
            N = self.measurement_model.rowCount()
            self.element_summary_text.setText("N={} lines".format(N))
            return None
        summary_dict = self.parent.session.summarize_spectral_models(organize_by_element=False)
        species = utils.element_to_species(elem)
        if species not in summary_dict:
            logger.debug("Cannot find species {} in summary_dict (keys={})".format(
                    species, summary_dict.keys()))
            self.element_summary_text.setText("ERROR {}".format(species))
        else:
            N, abund, stdev, stderr, XH, XFe = summary_dict[species]
            text = "N={1} A({0})={2:5.2f} σ({0})={5:4.2f} [{0}/H]={3:5.2f} [{0}/Fe I]={4:5.2f}"
            text = text.format(elem,N,abund,XH,XFe,stdev)
            self.element_summary_text.setText(text)
        
        return None

    def refresh_table(self):
        session = self.parent.session
        if session is None: return None
        #self._check_for_spectral_models()
        self.full_measurement_model.new_session(session)
        self.measurement_model.reset()
        self.measurement_view.update_session(session)
        self.populate_filter_combo_box()
        self.calculate_FeH()
        return None

    def refresh_plots(self):
        self.update_spectrum_figure(True)
        return None

    def fit_all_profiles(self):
        self._check_for_spectral_models()
        current_element_index = self.filter_combo_box.currentIndex()

        # Fit all acceptable
        num_unacceptable = 0
        for i,spectral_model in enumerate(self.full_measurement_model.spectral_models):
            if not spectral_model.is_acceptable:
                num_unacceptable += 1
                continue
            if isinstance(spectral_model, SpectralSynthesisModel):
                num_unacceptable += 1
                continue
                #try:
                #    res = spectral_model.fit()
                #except (ValueError, RuntimeError, TypeError) as e:
                #    logger.debug("Fitting error",spectral_model)
                #    logger.debug(e)
            elif isinstance(spectral_model, ProfileFittingModel):
                try:
                    res = spectral_model.fit()
                except (ValueError, RuntimeError, TypeError) as e:
                    logger.warn("Fitting error",spectral_model)
                    logger.warn(e)
        # If none are acceptable, then fit all
        if num_unacceptable == self.full_measurement_model.rowCount(None):
            logger.info("Found no acceptable spectral models, fitting all!")
            for i,spectral_model in enumerate(self.full_measurement_model.spectral_models):
                if isinstance(spectral_model, SpectralSynthesisModel):
                    continue
                    #try:
                    #    res = spectral_model.fit()
                    #except (ValueError, RuntimeError, TypeError) as e:
                    #    logger.debug("Fitting error",spectral_model)
                    #    logger.debug(e)
                elif isinstance(spectral_model, ProfileFittingModel):
                    try:
                        res = spectral_model.fit()
                    except (ValueError, RuntimeError, TypeError) as e:
                        logger.debug("Fitting error",spectral_model)
                        logger.debug(e)

        self.measurement_model.reset()
        self.populate_filter_combo_box()
        self.summarize_current_table()
        self.refresh_plots()
        self.filter_combo_box.setCurrentIndex(current_element_index)
        return None

    def measure_all(self):
        self._check_for_spectral_models()
        # Save this just to go back 
        current_element_index = self.filter_combo_box.currentIndex()
        try:
            #current_table_index = self.table_view.selectedIndexes()[-1]
            current_table_index = self.measurement_view.selectedIndexes()[-1]
        except:
            current_table_index = None

        # Gets abundances and uncertainties into session
        self.parent.session.measure_abundances()

        self.measurement_model.reset()
        self.populate_filter_combo_box()
        self.summarize_current_table()
        self.refresh_plots()

        self.filter_combo_box.setCurrentIndex(current_element_index)
        if current_table_index is not None:
            try:
                #self.table_view.selectRow(current_table_index)
                self.measurement_view.selectRow(current_table_index)
            except:
                logger.debug("Could not set index")
                pass
        return None

    def fit_one(self):
        spectral_model, proxy_index, index = self._get_selected_model(True)
        if spectral_model is None: return None
        try:
            res = spectral_model.fit()
        except (ValueError, RuntimeError) as e:
            logger.info("Fitting error",spectral_model)
            logger.info(e)
            return None
        #self.table_view.update_row(proxy_index.row())
        self.measurement_view.update_row(proxy_index.row())
        self.summarize_current_table()
        self.update_fitting_options()
        self.refresh_plots()
        return None

    def measure_one(self):
        spectral_model, proxy_index, index = self._get_selected_model(True)
        if spectral_model is None: return None

        self.parent.session.measure_abundances([spectral_model])

        #self.table_view.update_row(proxy_index.row())
        self.measurement_view.update_row(proxy_index.row())
        self.summarize_current_table()
        self.update_fitting_options()
        self.refresh_plots()
        return None

    def synthesize_current_model(self):
        spectral_model, proxy_index, index = self._get_selected_model(True)
        if spectral_model is None: return None
        spectral_model.update_fit_after_parameter_change()
        #self.table_view.update_row(proxy_index.row())
        self.measurement_view.update_row(proxy_index.row())
        self.summarize_current_table()
        self.update_fitting_options()
        self.refresh_plots()
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
                # Make sure callbacks are included in ui_mainwindow too!
                dialog = TransitionsDialog(self.parent.session,
                    callbacks=[self.parent.transition_dialog_callback])
                dialog.exec_()

                # Do we even have any spectral models now?
                for sm in self.parent.session.metadata.get("spectral_models", []):
                    if sm.use_for_stellar_composition_inference: break
                else:
                    return False
            else:
                return False
        return True


    def _get_selected_model(self, full_output=False):
        try:
            #proxy_index = self.table_view.selectionModel().selectedRows()[-1]
            proxy_index = self.measurement_view.selectionModel().selectedRows()[-1]
        except IndexError:
            return (None, None, None) if full_output else None
        index = self.measurement_model.mapToSource(proxy_index).row()
        model = self.parent.session.metadata["spectral_models"][index]
        return (model, proxy_index, index) if full_output else model

    def selected_model_changed(self):
        self.update_fitting_options()
        self.update_spectrum_figure(True)
        return None

    def update_spectrum_figure(self, redraw=False, reset_limits=True):
        ## If synthesis, label selected lines
        try:
            selected_model = self._get_selected_model()
        except IndexError:
            selected_transitions = None
            label_rv = None
        else:
            if isinstance(selected_model, SpectralSynthesisModel):
                selected_elem = self.synth_abund_table.get_selected_element()
                if selected_elem is not None:
                    transitions = selected_model.transitions
                    ii = np.logical_or(transitions["elem1"] == selected_elem,
                                       transitions["elem2"] == selected_elem)
                    if np.sum(ii) != 0:
                        selected_transitions = transitions[ii]
                        redraw=True # force redraw
                        reset_limits=False
                        label_rv = selected_model.metadata["manual_rv"]
                    else:
                        selected_transitions = None
                        label_rv = None
                else:
                    selected_transitions = None
                    label_rv = None
            else:
                selected_transitions = None
                label_rv = None
        # Update figure
        self.figure.update_spectrum_figure(redraw=redraw,
                                           reset_limits=reset_limits,
                                           label_transitions=selected_transitions,
                                           label_rv=label_rv)


    def update_fitting_options(self):
        try:
            selected_model = self._get_selected_model()
        except IndexError:
            return None
        if selected_model is None: return None
    
        if isinstance(selected_model, ProfileFittingModel):
            # Disconnect signals so it doesn't automatically fit
            self._disconnect_profile_signals()
            
            self.opt_tabs.setTabEnabled(0, True)
            self.opt_tabs.setCurrentIndex(0)

            # Increase view window by 2 percent for ease of masking edges
            self.edit_view_window.setText("{}".format(selected_model.metadata["window"]*1.02))
            self.edit_fit_window.setText("{}".format(selected_model.metadata["window"]))

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

            # Wavelength tolerance
            tolerance = selected_model.metadata.get("wavelength_tolerance", None)
            if tolerance is None:
                self.checkbox_wavelength_tolerance.setChecked(False)
                self.edit_wavelength_tolerance.setText("")
            else:
                self.checkbox_wavelength_tolerance.setChecked(True)
                self.edit_wavelength_tolerance.setText(
                    "{}".format(tolerance))

            self.checkbox_use_central_weighting.setChecked(
                selected_model.metadata["central_weighting"])

            self.combo_profile.setCurrentIndex(
                ["gaussian", "lorentzian", "voigt"].index(
                    selected_model.metadata["profile"]))

            self.edit_detection_sigma.setText("{}".format(
                selected_model.metadata["detection_sigma"]))
            self.edit_detection_pixels.setText("{}".format(
                selected_model.metadata["detection_pixels"]))

            try:
                # HACK
                self.checkbox_use_antimasks.setChecked(
                    selected_model.metadata["antimask_flag"])
            except KeyError: # HACK Old SMH sessions will not load with antimask_flag
                self.checkbox_use_antimasks.setChecked(False)

            try:
                # HACK
                self.checkbox_upper_limit.setChecked(
                    selected_model.metadata["is_upper_limit"])
            except KeyError: # HACK Old SMH sessions may not have upper limit flag
                self.checkbox_upper_limit.setChecked(False)

            # Reconnect signals
            self._connect_profile_signals()
        else:
            self.opt_tabs.setTabEnabled(0, False)

        # Synthesis options.
        if isinstance(selected_model, SpectralSynthesisModel):
            # Disconnect signals so it doesn't automatically fit
            self._disconnect_profile_signals()

            self.opt_tabs.setTabEnabled(1, True)
            self.opt_tabs.setCurrentIndex(1)

            self.edit_view_window.setText("{}".format(selected_model.metadata["window"]))
            self.edit_fit_window.setText("{}".format(selected_model.metadata["window"]))

            # Continuum order.
            continuum_order = selected_model.metadata["continuum_order"]
            if continuum_order < 0:
                self.checkbox_continuum_2.setChecked(False)
                self.combo_continuum_2.setEnabled(False)
                self.edit_manual_continuum.setEnabled(True)
            else:
                self.checkbox_continuum_2.setChecked(True)
                self.combo_continuum_2.setEnabled(True)
                self.combo_continuum_2.setCurrentIndex(continuum_order)
                self.edit_manual_continuum.setEnabled(False)
            self.edit_manual_continuum.setText(
                "{}".format(selected_model.metadata["manual_continuum"]))

            # Radial velocity tolerance.
            vrad_tolerance = selected_model.metadata.get("velocity_tolerance", None)
            if vrad_tolerance is None:
                self.checkbox_vrad_tolerance_2.setChecked(False)
                self.edit_vrad_tolerance_2.setEnabled(False)
                self.edit_manual_rv.setEnabled(True)
            else:
                self.checkbox_vrad_tolerance_2.setChecked(True)
                self.edit_vrad_tolerance_2.setEnabled(True)
                self.edit_vrad_tolerance_2.setText("{}".format(vrad_tolerance))
                self.edit_manual_rv.setEnabled(False)
            self.edit_manual_rv.setText("{:.4f}".format(selected_model.metadata["manual_rv"]))

            self.edit_initial_abundance_bound.setText(
                "{}".format(selected_model.metadata["initial_abundance_bounds"]))
            
            if selected_model.metadata["smoothing_kernel"]:
                self.checkbox_model_smoothing.setChecked(True)
                self.edit_manual_smoothing.setEnabled(False)
            else:
                self.checkbox_model_smoothing.setChecked(False)
                self.edit_manual_smoothing.setEnabled(True)
            self.edit_manual_smoothing.setText(
                "{:.4f}".format(selected_model.metadata["manual_sigma_smooth"]))

            # sigma smooth tolerance needs implementing.
            self.synth_abund_table_model.load_new_model(selected_model)

            try:
                # HACK
                self.checkbox_upper_limit_2.setChecked(
                    selected_model.metadata["is_upper_limit"])
            except KeyError: # HACK Old SMH sessions may not have upper limit flag
                self.checkbox_upper_limit_2.setChecked(False)

            # Reconnect signals
            self._connect_profile_signals()
        else:
            self.opt_tabs.setTabEnabled(1, False)

        return None

    ###############################
    # FITTING OPTION UPDATE METHODS
    ###############################

    ### Profile update
    def update_edit_view_window(self):
        """ The wavelength window was updated """
        try:
            window = float(self.edit_view_window.text())
        except:
            return None
        else:
            transitions = self._get_selected_model().transitions
            xlim = (transitions["wavelength"][0] - window,
                    transitions["wavelength"][-1] + window)
            self.ax_spectrum.set_xlim(xlim)
            self.ax_residual.set_xlim(xlim)
            self.figure.reset_zoom_limits()
            self.figure.draw()
        return None
    def update_edit_fit_window(self):
        """ The wavelength window was updated """
        model = self._get_selected_model()
        try:
            window = float(self.edit_fit_window.text())
        except:
            return None
        else:
            model.metadata["window"] = window
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
        try:
            value = float(self.edit_detection_sigma.text())
        except:
            value = None
        self._get_selected_model().metadata["detection_sigma"] = value
        return None
    def update_detection_pixels(self):
        """ The number of pixels to qualify a detection has been updated. """
        try:
            value = int(self.edit_detection_pixels.text())
        except:
            value = None
        self._get_selected_model().metadata["detection_pixels"] = value
        return None
    def clicked_checkbox_use_central_weighting(self):
        """ The checkbox to use central weighting has been clicked. """
        self._get_selected_model().metadata["central_weighting"] \
            = self.checkbox_use_central_weighting.isChecked()
        return None
    def clicked_checkbox_use_antimasks(self):
        """ The checkbox to use antimasks has been clicked. """
        # TODO
        self._get_selected_model().metadata["antimask_flag"] \
            = self.checkbox_use_antimasks.isChecked()
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
        try:
            value = float(self.edit_wavelength_tolerance.text())
        except:
            value = None
        self._get_selected_model().metadata["wavelength_tolerance"] = value
        return None
    def clicked_checkbox_upper_limit(self):
        """ The checkbox to set as upper limit has been clicked. """
        spectral_model, proxy_index, index = self._get_selected_model(True)
        spectral_model.metadata["is_upper_limit"] \
            = self.checkbox_upper_limit.isChecked()
        #self.table_view.update_row(proxy_index.row())
        self.measurement_view.update_row(proxy_index.row())
        self.summarize_current_table()
        self.refresh_plots()
        return None
    def clicked_btn_clear_masks(self):
        spectral_model = self._get_selected_model()
        spectral_model.metadata["mask"] = []
        spectral_model.metadata.pop("fitted_result", None)
        self.fit_one()
        return None

    ### Synthesis update
    def update_edit_view_window_2(self):
        """ The wavelength window was updated """
        try:
            window = float(self.edit_view_window_2.text())
        except:
            return None
        else:
            transitions = self._get_selected_model().transitions
            xlim = (transitions["wavelength"][0] - window,
                    transitions["wavelength"][-1] + window)
            self.ax_spectrum.set_xlim(xlim)
            self.ax_residual.set_xlim(xlim)
            self.figure.reset_zoom_limits()
            self.figure.draw()
        return None
    def update_edit_fit_window_2(self):
        """ The wavelength window was updated """
        model = self._get_selected_model()
        try:
            window = float(self.edit_fit_window_2.text())
        except:
            return None
        else:
            model.metadata["window"] = window
        return None
    def clicked_checkbox_continuum_2(self):
        """ The checkbox for modeling the continuum was clicked. """
        if self.checkbox_continuum_2.isChecked():
            self.combo_continuum_2.setEnabled(True)
            self.update_continuum_order()
            self.edit_manual_continuum.setEnabled(False)
        else:
            self._get_selected_model().metadata["continuum_order"] = -1
            self.combo_continuum_2.setEnabled(False)
            self.edit_manual_continuum.setEnabled(True)
        return None
    def update_continuum_order_2(self):
        """ The continuum order to use in model fitting was changed. """
        self._get_selected_model().metadata["continuum_order"] \
            = int(self.combo_continuum_2.currentText())
        return None
    def update_edit_manual_continuum(self):
        try:
            value = float(self.edit_manual_continuum.text())
        except:
            value = None
        else:
            selected_model = self._get_selected_model()
            selected_model.metadata["manual_continuum"] = value
            selected_model.update_fit_after_parameter_change()
            self.update_spectrum_figure(redraw=True)
        return None
    def clicked_checkbox_vrad_tolerance_2(self):
        """ The checkbox for velocity tolerance was clicked. """
        if self.checkbox_vrad_tolerance_2.isChecked():
            self.edit_vrad_tolerance_2.setEnabled(True)
            self.update_vrad_tolerance_2()
            self.edit_manual_rv.setEnabled(False)
        else:
            self.edit_vrad_tolerance_2.setEnabled(False)
            self._get_selected_model().metadata["velocity_tolerance"] = None
            self.edit_manual_rv.setEnabled(True)
        return None
    def update_vrad_tolerance_2(self):
        """ The tolerance on radial velocity was updated. """
        try:
            value = float(self.edit_vrad_tolerance_2.text())
        except:
            value = None
        self._get_selected_model().metadata["velocity_tolerance"] = value
        return None
    def update_manual_rv(self):
        try:
            value = float(self.edit_manual_rv.text())
        except:
            value = None
        else:
            selected_model = self._get_selected_model()
            selected_model.metadata["manual_rv"] = value
            selected_model.update_fit_after_parameter_change()
            self.update_spectrum_figure(redraw=True)
        return None
    def update_initial_abundance_bound(self):
        """ The initial abundance bound has been updated. """
        try:
            value = float(self.edit_initial_abundance_bound.text())
        except:
            value = None
        else:
            self._get_selected_model().metadata["initial_abundance_bounds"] \
            = value
        return None
    def clicked_checkbox_model_smoothing(self):
        """ The checkbox to smooth the model spectrum has been clicked. """
        if self.checkbox_model_smoothing.isChecked():
            self._get_selected_model().metadata["smoothing_kernel"] = True
            self.edit_manual_smoothing.setEnabled(False)
        else:
            self._get_selected_model().metadata["smoothing_kernel"] = False
            self.edit_manual_smoothing.setEnabled(True)
        return None
    def update_manual_smoothing(self):
        try:
            value = float(self.edit_manual_smoothing.text())
        except:
            value = None
        else:
            selected_model = self._get_selected_model()
            selected_model.metadata["manual_sigma_smooth"] = value
            selected_model.update_fit_after_parameter_change()
            self.update_spectrum_figure(redraw=True)
        return None
    def clicked_checkbox_upper_limit_2(self):
        """ The checkbox to set as upper limit has been clicked. """
        spectral_model, proxy_index, index = self._get_selected_model(True)
        spectral_model.metadata["is_upper_limit"] \
            = self.checkbox_upper_limit_2.isChecked()
        #self.table_view.update_row(proxy_index.row())
        self.measurement_view.update_row(proxy_index.row())
        self.summarize_current_table()
        self.refresh_plots()
        return None
    def clicked_btn_update_abund_table(self):
        selected_model = self._get_selected_model()
        if selected_model is None: return None
        assert isinstance(selected_model, SpectralSynthesisModel), selected_model
        summary_dict = self.parent.session.summarize_spectral_models(organize_by_element=True)
        
        self.synth_abund_table.model().beginResetModel()

        # Fill in fixed abundances
        for elem in selected_model.metadata["rt_abundances"]:
            try:
                selected_model.metadata["rt_abundances"][elem] = summary_dict[elem][1]
            except KeyError:
                selected_model.metadata["rt_abundances"][elem] = np.nan

        # Fill in fixed abundances
        try:
            fitted_result = selected_model.metadata["fitted_result"]
        except KeyError:
            logger.info("Run at least one fit before setting abundances of "
                  "fitted element {}!".format(elem))
        else:
            for i,elem in enumerate(selected_model.elements):
                try:
                    abund = summary_dict[elem][1]
                except KeyError:
                    logger.warn("No abundance found for {}, using nan".format(elem))
                    abund = np.nan
                key = "log_eps({})".format(elem)
                fitted_result[0][key] = abund
                fitted_result[2]["abundances"][i] = abund

        self.synth_abund_table.model().endResetModel()

        logger.debug(summary_dict)
        return None

    def clicked_export_synthesis(self):
        ## Get current spectral model, make sure it is a synthesis
        spectral_model = self._get_selected_model()
        if not isinstance(spectral_model, SpectralSynthesisModel): 
            logger.info("Must select a synthesis spectral model to export")
            return
        ## ask for synthesis output filename
        synth_path, _ = QtGui.QFileDialog.getSaveFileName(self,
                caption="Enter synthesis output filename", dir="") #, filter="*.txt")
        if not synth_path: return
        ## Ask for data output filename
        data_path, _ = QtGui.QFileDialog.getSaveFileName(self,
                caption="Enter data output filename", dir="") #, filter="*.txt")
        if not data_path: return
        ## Ask for parameter output filename
        param_path, _ = QtGui.QFileDialog.getSaveFileName(self,
                caption="Enter parameter output filename", dir="") #, filter="*.txt")
        if not param_path: return
        ## Export
        spectral_model.export_fit(synth_path, data_path, param_path)
        logger.info("Exported to {}, {}, and {}".format(synth_path, data_path, param_path))
        return

    def refresh_current_model(self):
        spectral_model, proxy_index, index = self._get_selected_model(True)
        #self.table_view.update_row(proxy_index.row())
        self.measurement_view.update_row(proxy_index.row())
        self.update_fitting_options()

#class SpectralModelsTableView(SpectralModelsTableViewBase):
#    def sizeHint(self):
#        #return QtCore.QSize(240,100)
#        return QtCore.QSize(125,100)
#
#    def minimumSizeHint(self):
#        return QtCore.QSize(125,0)
#
#    def refresh_gui(self):
#        self.parent.summarize_current_table()
#        self.parent.update_fitting_options()
#        self.parent.refresh_plots()
#        return None
#
#    def fit_selected_models(self, proxy_indices):
#        """ Fit the selected spectral models. """
#        # Fit the models one by one
#        for proxy_index in proxy_indices:
#            index = self.model().mapToSource(proxy_index).row()
#            self.parent.parent.session.metadata["spectral_models"][index].fit()
#
#            # Update the data model, view
#            self.update_row(proxy_index.row())
#
#        self.refresh_gui()
#        return None
#    
#    def measure_selected_models(self, proxy_indices):
#        """ Fit the selected spectral models. """
#
#        # Get list of selected spectral models
#        spectral_models = []
#        for proxy_index in proxy_indices:
#            index = self.model().mapToSource(proxy_index).row()
#            spectral_models.append(self.parent.parent.session.metadata["spectral_models"][index])
#
#        # Fit abundances
#        self.parent.parent.session.measure_abundances(spectral_models)
#
#        # Update the data model
#        start = time.time()
#        for proxy_index in proxy_indices:
#            self.update_row(proxy_index.row())
#        logger.debug("Time to update data model: {:.1f}".format(time.time()-start))
#
#        self.refresh_gui()
#        return None
#
#    def mark_selected_models_as_acceptable(self, proxy_indices):
#        proxy_model = self.parent.measurement_model
#        full_model = proxy_model.sourceModel()
#        for proxy_index in proxy_indices:
#            full_index = proxy_model.mapToSource(proxy_index)
#            full_model.setData(full_index, 2, refresh_view=False)
#        self.refresh_gui()
#        return None
#    def mark_selected_models_as_unacceptable(self, proxy_indices):
#        proxy_model = self.parent.measurement_model
#        full_model = proxy_model.sourceModel()
#        for proxy_index in proxy_indices:
#            full_index = proxy_model.mapToSource(proxy_index)
#            full_model.setData(full_index, 0, refresh_view=False)
#        self.refresh_gui()
#        return None
#
#    def set_fitting_option_value(self, proxy_indices, key, value,
#                                 valid_for_profile=False,
#                                 valid_for_synth=False):
#        num_fit = 0
#        num_unacceptable = 0
#        num_profile_models = 0
#        num_synthesis_models = 0
#        for proxy_index in proxy_indices:
#            idx = self.model().mapToSource(proxy_index).row()
#            spectral_model \
#                = self.parent.parent.session.metadata["spectral_models"][idx]
#            run_fit = False
#            if not spectral_model.is_acceptable: 
#                num_unacceptable += 1
#                continue
#            if valid_for_profile and isinstance(spectral_model,ProfileFittingModel):
#                num_profile_models += 1
#                spectral_model.metadata[key] = value
#                run_fit = True
#            if valid_for_synth and isinstance(spectral_model,SpectralSynthesisModel):
#                num_synthesis_models += 1
#                spectral_model.metadata[key] = value
#                run_fit = True
#                
#            if run_fit and "fitted_result" in spectral_model.metadata:
#                num_fit += 1
#                spectral_model.fit()
#                self.update_row(proxy_index.row())
#        logger.debug("Changed {0}={1}, fit {2} out of {3} models ({4} profile, {5} synth, skipped {6} unacceptable)".format(\
#                key, value, num_fit, len(proxy_indices), num_profile_models, num_synthesis_models, num_unacceptable))
#        self.refresh_gui()
#        return None
#
#
#class SpectralModelsTableModel(SpectralModelsTableModelBase):
#    def data(self, index, role):
#        """
#        Display the data.
#
#        :param index:
#            The table index.
#
#        :param role:
#            The display role.
#        """
#
#        if not index.isValid():
#            return None
#
#        if role==QtCore.Qt.FontRole:
#            return _QFONT
#
#        column = index.column()
#        spectral_model = self.spectral_models[index.row()]
#
#        if  column == 0 \
#        and role in (QtCore.Qt.DisplayRole, QtCore.Qt.CheckStateRole):
#            value = spectral_model.is_acceptable
#            if role == QtCore.Qt.CheckStateRole:
#                return QtCore.Qt.Checked if value else QtCore.Qt.Unchecked
#            else:
#                return None
#        elif column == 1:
#            value = spectral_model._repr_wavelength
#
#        elif column == 2: #abundance
#            try:
#                abundances \
#                    = spectral_model.metadata["fitted_result"][2]["abundances"]
#
#            except (IndexError, KeyError):
#                value = ""
#
#            else:
#                if len(abundances) == 1:
#                    value = "{0:.2f}".format(abundances[0])
#                else:
#                    assert isinstance(spectral_model, SpectralSynthesisModel), spectral_model
#                    current_element = self.parent.filter_combo_box.currentText()
#                    if current_element=="":
#                        value = ""
#                    elif current_element=="All":
#                        try:
#                            value = "; ".join(["{}".format(abund) \
#                                               for abund in spectral_model.abundances])
#                        except TypeError:
#                            value = ""
#                    else:
#                        # HACK for molecules
#                        if "-" in current_element: _elem = utils.species_to_element(utils.element_to_atomic_number(current_element)).split()[0]
#                        else: _elem = current_element.split()[0]
#                        for i,elem in enumerate(spectral_model.elements):
#                            if _elem == elem: break
#                        else:
#                            raise ValueError("{} ({}) not in {}".format(current_element, _elem, spectral_model.elements))
#                        value = "{0:.2f}".format(abundances[i])
#        elif column in [3, 4]: #EW, REW
#            if isinstance(spectral_model, ProfileFittingModel):
#                try:
#                    result = spectral_model.metadata["fitted_result"][2]
#                    equivalent_width = result["equivalent_width"][0]
#                except:
#                    equivalent_width = np.nan
#    
#                if column == 3:
#                    value = "{0:.1f}".format(1000 * equivalent_width) \
#                        if np.isfinite(equivalent_width) else ""
#                if column == 4:
#                    value = "{:.2f}".format(np.log10(equivalent_width/float(spectral_model._repr_wavelength))) \
#                        if np.isfinite(equivalent_width) else ""
#            elif isinstance(spectral_model, SpectralSynthesisModel):
#                if column == 3:
#                    value = ""
#                if column == 4:
#                    # HACK
#                    value = "-4.0"
#                
#        elif column == 5: #abundance err
#            if isinstance(spectral_model, ProfileFittingModel):
#                try:
#                    result = spectral_model.metadata["fitted_result"][2]
#                    err = result["abundance_uncertainties"][0]
#                    value = "{:.2f}".format(err)
#                except:
#                    value = ""
#            elif isinstance(spectral_model, SpectralSynthesisModel):
#                current_element = self.parent.filter_combo_box.currentText()
#                if current_element=="":
#                    value = ""
#                elif current_element=="All":
#                    value = ""
#                else:
#                    # HACK for molecules
#                    if "-" in current_element: _elem = utils.species_to_element(utils.element_to_atomic_number(current_element)).split()[0]
#                    else: _elem = current_element.split()[0]
#                    for i,elem in enumerate(spectral_model.elements):
#                        if _elem == elem: break
#                    else:
#                        raise ValueError("{} not in {}".format(current_element, spectral_model.elements))
#                    try:
#                        covar = spectral_model.metadata["fitted_result"][1]
#                    except:
#                        value = ""
#                    else:
#                        err = np.sqrt(covar[i,i])
#                        value = "{:.2f}".format(err)
#        elif column == 6: #EW err
#            try:
#                result = spectral_model.metadata["fitted_result"][2]
#                err = 1000.*np.nanmax(np.abs(result["equivalent_width"][1:3]))
#                value = "{:.2f}".format(err)
#            except:
#                value = ""
#        elif column == 7:
#            if isinstance(spectral_model, SpectralSynthesisModel):
#                value = ""
#            else:
#                try:
#                    loggf = spectral_model.transitions[0]['loggf']
#                    value = "{:6.3f}".format(loggf)
#                except:
#                    value = ""
#        elif column == 8:
#            value = "; ".join(["{}".format(element) \
#                      for element in spectral_model.elements])
#
#        return value if role == QtCore.Qt.DisplayRole else None
#
#    def setData(self, index, value, role=QtCore.Qt.DisplayRole, refresh_view=True):
#        value = super(SpectralModelsTableModel, self).setData(index, value, role)
#        if index.column() != 0: return False
#        
#        proxy_index = self.parent.table_view.model().mapFromSource(index)
#        proxy_row = proxy_index.row()
#        self.parent.table_view.rowMoved(proxy_row, proxy_row, proxy_row)
#
#        if refresh_view:
#            self.parent.summarize_current_table()
#            self.parent.refresh_plots()
#
#        return value

class SynthesisAbundanceTableView(QtGui.QTableView):
    """
    Make a small table view
    """
    def sizeHint(self):
        return QtCore.QSize(100,100)
    def minimumSizeHint(self):
        return QtCore.QSize(100,0)
    def get_selected_element(self):
        try:
            index = self.selectionModel().selectedIndexes()[-1]
            elem = self.model().elem_order[index.row()]
            return elem
        except:
            logger.exception("Could not get selected element")
            return None
class SynthesisAbundanceTableModel(QtCore.QAbstractTableModel):
    """
    Editable table of abundances to synthesize
    """
    def __init__(self, parent, *args):
        super(SynthesisAbundanceTableModel, self).__init__(parent, *args)
        self.parent = parent
        self.spectral_model = None
        self.elem_order = None
        self.num_fit_elems = 0

    def load_new_model(self, spectral_model):
        """
        Call this to reset the table with a new spectral model
        """
        self.parent.calculate_FeH()
        self.beginResetModel()
        self.spectral_model = spectral_model
        # Sort table by Z
        if spectral_model is not None:
            # First rows in table are fit elems
            # Other rows are rt_abundances
            self.num_fit_elems = len(spectral_model.elements)

            elems = spectral_model.metadata["rt_abundances"].keys()
            Zs = [utils.element_to_atomic_number(elem) for elem in elems]
            sorted_indices = np.argsort(Zs)
            # Put in rt_abundances indices
            self.elem_order = dict(zip(self.num_fit_elems + np.arange(len(elems)), \
                                       np.array(elems)[sorted_indices]))
            # Put in parameters indices
            for i,elem in enumerate(spectral_model.elements):
                self.elem_order[i] = elem
        else:
            self.num_fit_elems = 0
            self.elem_order = None
        self.endResetModel()
        return None
    def rowCount(self, parent):
        try:
            return self.num_fit_elems + \
                   len(self.spectral_model.metadata["rt_abundances"])
        except Exception as e:
            logger.info(e)
            return 0
    def columnCount(self, parent):
        return 3
    def data(self, index, role):
        if role==QtCore.Qt.FontRole:
            return _QFONT
        if not index.isValid() or role != QtCore.Qt.DisplayRole:
            return None
        if self.spectral_model is None: return None
        elem = self.elem_order[index.row()]
        if index.column()==0: 
            return elem

        if elem in self.spectral_model.metadata["rt_abundances"]:
            logeps = self.spectral_model.metadata["rt_abundances"][elem]
        elif "fitted_result" not in self.spectral_model.metadata:
            logeps = np.nan
        else:
            fitted_result = self.spectral_model.metadata["fitted_result"]
            key = "log_eps({})".format(elem)
            assert key in fitted_result[0], "{} {}".format(key,fitted_result[0])
            logeps = fitted_result[0][key]
        
        if index.column()==1:
            return "{:.3f}".format(logeps)
        elif index.column()==2:
            return "{:.3f}".format(logeps-solar_composition(elem)-self.parent.FeH)
        return None
    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal \
        and role == QtCore.Qt.DisplayRole:
            if col==0: return "El."
            if col==1: return "A(X)"
            if col==2: return "[X/Fe]"
        if role==QtCore.Qt.FontRole:
            return _QFONT
        return None
    def setData(self, index, value, role):
        # [X/Fe] and logeps appear to automatically update each other!
        if index.column()==0: return False
        if self.spectral_model is None: return False
        # Modify the spectral model abundance
        elem = self.elem_order[index.row()]
        if elem in self.spectral_model.metadata["rt_abundances"]:
            try:
                value = float(value)
            except ValueError:
                return False
            else:
                if index.column() == 1:
                    self.spectral_model.metadata["rt_abundances"][elem] = value
                elif index.column() == 2:
                    self.spectral_model.metadata["rt_abundances"][elem] = \
                        value + solar_composition(elem) + self.parent.FeH
                return True
        elif elem in self.spectral_model.elements:
            try:
                value = float(value)
            except ValueError:
                return False

            # HACK
            # Replace abundance in both fitted parameters and abundances
            # SpectralSynthesisModel.__call__ uses fitted parameters to synth
            try:
                fitted_result = self.spectral_model.metadata["fitted_result"]
            except KeyError:
                logger.info("Run at least one fit before setting abundances!")
                return False
            else:
                key = "log_eps({})".format(elem)
                if index.column() == 1:
                    pass #value = value
                elif index.column() == 2:
                    value = value + solar_composition(elem) + self.parent.FeH

                fitted_result[0][key] = value
                for i,_elem in enumerate(self.spectral_model.elements):
                    if _elem == elem: break
                else: raise ValueError(elem+" "+str(self.spectral_model.elements))
                fitted_result[2]["abundances"][i] = value
                return True
        else:
            raise ValueError(elem+" "+str(self.spectral_model.elements))
                
    def flags(self, index):
        if not index.isValid():
            return None
        flags = QtCore.Qt.ItemIsEnabled|\
                QtCore.Qt.ItemIsSelectable
        if index.column()==1 or index.column(): #abundance
            flags |= QtCore.Qt.ItemIsEditable
        return flags
