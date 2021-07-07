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
from PySide2 import (QtCore, QtGui as QtGui2, QtWidgets as QtGui)
import time
from copy import deepcopy

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
        QtGui2.QFont.insertSubstitution(*substitute)

_QFONT = QtGui2.QFont("Helvetica Neue", 10)
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
        self.figure.mpl_connect("key_press_event", self.key_press_selectcheck)
        ## Stuff for extra synthesis
        self.extra_spec_1 = self.ax_spectrum.plot([np.nan],[np.nan], ls='-', color='#cea2fd', lw=1.5, zorder=9999)[0]
        self.extra_spec_2 = self.ax_spectrum.plot([np.nan],[np.nan], ls='-', color='#ffb07c', lw=1.5, zorder=9999)[0]
        
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
                                                     "fwhm",
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
        self.measurement_model.add_callback_after_setData(self.summarize_current_table)
        _ = self.measurement_view.selectionModel()
        _.selectionChanged.connect(self.selected_model_changed)
        self.measurement_view.setSizePolicy(QtGui.QSizePolicy(
            QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.MinimumExpanding))
        self.btn_filter_acceptable = btn_filter
        self.btn_filter_acceptable.clicked.connect(self.refresh_plots)
        self.btn_refresh = btn_refresh
        self.btn_refresh.clicked.connect(self.new_session_loaded)
        bot_lhs_layout.addLayout(vbox)

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
                line.setValidator(QtGui2.QIntValidator(bot, top, line))
            else:
                line.setValidator(QtGui2.QDoubleValidator(bot, top, dec, line))
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
            line.setValidator(QtGui2.QDoubleValidator(bot, top, dec, line))
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
                                                 0, 100, 2)
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
                                                 -1000, 1000, 1)
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
        self.synth_abund_table.verticalHeader().setSectionResizeMode(QtGui.QHeaderView.Fixed)
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

        hbox = QtGui.QHBoxLayout()
        hbox.setSpacing(0)
        hbox.setContentsMargins(0,0,0,0)
        hbox2 = QtGui.QHBoxLayout()
        hbox2.setSpacing(0)
        hbox2.setContentsMargins(0,0,0,0)
        self.checkbox_upper_limit_2 = QtGui.QCheckBox(self.tab_synthesis)
        self.checkbox_upper_limit_2.setText("Upper Limit")
        self.checkbox_upper_limit_2.setFont(_QFONT)
        hbox.addWidget(self.checkbox_upper_limit_2)
        label = QtGui.QLabel(self.tab_synthesis)
        label.setText("Find")
        label.setFont(_QFONT)
        hbox2.addWidget(label)
        line = QtGui.QLineEdit(self.tab_synthesis)
        line.setMinimumSize(QtCore.QSize(20, 0))
        line.setMaximumSize(QtCore.QSize(25, _ROWHEIGHT))
        line.setFont(_QFONT)
        line.setValidator(QtGui2.QDoubleValidator(0, 10, 1, line))
        line.setText("3.0")
        self.edit_ul_sigma = line
        hbox2.addWidget(self.edit_ul_sigma)
        label = QtGui.QLabel(self.tab_synthesis)
        label.setText(u"\u03C3")
        label.setFont(_QFONT)
        hbox2.addWidget(label)
        hbox.addLayout(hbox2)
        self.btn_find_upper_limit = QtGui.QPushButton(self.tab_synthesis)
        self.btn_find_upper_limit.setText("Upper Limit")
        self.btn_find_upper_limit.setFont(_QFONT)
        hbox.addWidget(self.btn_find_upper_limit)
        vbox_rhs.addLayout(hbox)

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
        ## TODO NEXT: in PySide2, there is a bug with disconnecting signals.
        ## Recommended to switch to PyQt5 for this, which introduces other issues...
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
        self.btn_find_upper_limit.clicked.connect(
            self.clicked_btn_find_upper_limit)
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
            (self.btn_find_upper_limit.clicked,self.clicked_btn_find_upper_limit),
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
        self.measurement_model.beginResetModel()
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
        self.measurement_model.endResetModel()
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
        self.measurement_model.beginResetModel()
        self.full_measurement_model.new_session(session)
        self.measurement_model.endResetModel()
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
        self.measurement_model.beginResetModel()
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

        self.measurement_model.endResetModel()
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
            current_table_index = self.measurement_view.selectedIndexes()[-1]
        except:
            current_table_index = None

        # Gets abundances and uncertainties into session
        self.measurement_model.beginResetModel()
        self.parent.session.measure_abundances()

        self.measurement_model.endResetModel()
        self.populate_filter_combo_box()
        self.summarize_current_table()
        self.refresh_plots()

        self.filter_combo_box.setCurrentIndex(current_element_index)
        if current_table_index is not None:
            try:
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
        self.measurement_view.update_row(proxy_index.row())
        self.summarize_current_table()
        self.update_fitting_options()
        self.refresh_plots()
        if self.parent.session.setting("bring_to_top_after_fit", False):
            self.parent.raise_()
            self.parent.activateWindow()
            self.parent.showNormal()
        return None

    def measure_one(self):
        spectral_model, proxy_index, index = self._get_selected_model(True)
        if spectral_model is None: return None

        self.parent.session.measure_abundances([spectral_model])

        self.measurement_view.update_row(proxy_index.row())
        self.summarize_current_table()
        self.update_fitting_options()
        self.refresh_plots()
        return None

    def synthesize_current_model(self):
        # When we get selected model, it erases the extra_abundances.
        # So let's cache it then put it back in...
        extra_abundances = self.synth_abund_table_model.get_extra_abundances()
        
        spectral_model, proxy_index, index = self._get_selected_model(True)
        if spectral_model is None: return None
        spectral_model.update_fit_after_parameter_change()
        self.measurement_view.update_row(proxy_index.row())
        self.summarize_current_table()
        self.update_fitting_options()
        #self.refresh_plots()
        self.update_spectrum_figure(redraw=True,reset_limits=False)
        ### TODO
        try:
            fitted_result = spectral_model.metadata["fitted_result"]
        except:
            pass
        else:
            if extra_abundances is not None:
                logger.debug(extra_abundances)
                abundances1 = deepcopy(spectral_model.metadata["rt_abundances"])
                for i, elem in enumerate(spectral_model.elements):
                    abundances1[elem] = fitted_result[-1]["abundances"][i]
                abundances2 = deepcopy(abundances1)
                for elem, abunddiff in extra_abundances[0].items():
                    abundances1[elem] += abunddiff
                    self.synth_abund_table_model.extra_abundances[elem][0] = abunddiff
                for elem, abunddiff in extra_abundances[1].items():
                    abundances2[elem] += abunddiff
                    self.synth_abund_table_model.extra_abundances[elem][1] = abunddiff
                x, y = spectral_model.get_synth(abundances1)
                self.extra_spec_1.set_data([x,y])
                x, y = spectral_model.get_synth(abundances2)
                self.extra_spec_2.set_data([x,y])
            self.figure.draw()
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
        self.extra_spec_1.set_data([[np.nan], [np.nan]])
        self.extra_spec_2.set_data([[np.nan], [np.nan]])
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

            #self.edit_view_window_2.setText("{}".format(selected_model.metadata["window"]))
            self.edit_fit_window_2.setText("{}".format(selected_model.metadata["window"]))

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
            selected_model.update_fit_after_parameter_change(synthesize=False)
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
            selected_model.update_fit_after_parameter_change(synthesize=False)
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
            selected_model.update_fit_after_parameter_change(synthesize=False)
            self.update_spectrum_figure(redraw=True)
        return None
    def clicked_checkbox_upper_limit_2(self):
        """ The checkbox to set as upper limit has been clicked. """
        spectral_model, proxy_index, index = self._get_selected_model(True)
        spectral_model.metadata["is_upper_limit"] \
            = self.checkbox_upper_limit_2.isChecked()
        self.measurement_view.update_row(proxy_index.row())
        self.summarize_current_table()
        self.refresh_plots()
        return None
    def clicked_btn_find_upper_limit(self):
        """ The button to find an upper limit has been clicked. """
        spectral_model, proxy_index, index = self._get_selected_model(True)
        # Find the upper limit 
        try:
            sigma = round(float(self.edit_ul_sigma.text()),1)
        except:
            logger.debug("Invalid sigma for finding limit")
            return None
        upper_limit = spectral_model.find_upper_limit(sigma=sigma, start_at_current=True)
        # Refresh GUI
        self.measurement_view.update_row(proxy_index.row())
        self.summarize_current_table()
        self.update_fitting_options()
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
        ## Ask for linelist output filename
        linelist_path, _ = QtGui.QFileDialog.getSaveFileName(self,
                caption="Enter linelist output filename", dir="") #, filter="*.txt")
        if not param_path: return
        ## Export
        spectral_model.export_fit(synth_path, data_path, param_path)
        spectral_model.export_line_list(linelist_path)
        logger.info("Exported to {}, {}, {}, {}".format(synth_path, data_path, param_path, linelist_path))
        return

    def refresh_current_model(self):
        spectral_model, proxy_index, index = self._get_selected_model(True)
        if spectral_model is None: return None
        self.measurement_view.update_row(proxy_index.row())
        self.update_fitting_options()

    def key_press_selectcheck(self, event):
        if "'"+event.key+"'" == "' '":
            model, proxy_index, index = self._get_selected_model(True)
            model.is_acceptable = np.logical_not(model.is_acceptable)
            self.measurement_view.update_row(proxy_index.row())
            self.update_spectrum_figure(redraw=True)
            return None
        if event.key in ["f", "F"]:
            model, proxy_index, index = self._get_selected_model(True)
            model.user_flag = np.logical_not(model.user_flag)
            self.measurement_view.update_row(proxy_index.row())
            return None
        if event.key not in ["up","down","j", "J", "k", "K"]: return None
        try:
            proxy_index_row = self.measurement_view.selectionModel().selectedRows()[-1].row()
        except IndexError:
            return None
        if event.key in ["up", "j", "J"]: proxy_index_row -= 1
        if event.key in ["down", "k", "K"]: proxy_index_row += 1
        self.measurement_view.selectRow(proxy_index_row)
        return None

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
            index = self.selectionModel().selectedIndexes()
            if len(index)==0: return None
            index = index[-1]
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
        # Extra synthesis table info
        self.extra_abundances = {}

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

            #elems = [k.decode() for k in spectral_model.metadata["rt_abundances"].keys()]
            elems = list(spectral_model.metadata["rt_abundances"].keys())
            Zs = [utils.element_to_atomic_number(elem) for elem in elems]
            sorted_indices = np.argsort(Zs)
            # Put in rt_abundances indices
            self.elem_order = dict(zip(self.num_fit_elems + np.arange(len(elems)), \
                                       np.array(elems)[sorted_indices]))
            # Put in parameters indices and extra_abundances
            for elem in self.extra_abundances:
                self.extra_abundances[elem] = [None, None]
            for i,elem in enumerate(spectral_model.elements):
                self.elem_order[i] = elem
                self.extra_abundances[elem] = [None, None]
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
        return 5
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
        elif (index.column()==3) or (index.column()==4):
            if (elem not in self.extra_abundances) or (self.extra_abundances[elem][index.column()-3] is None): return ""
            return "{:.3f}".format(self.extra_abundances[elem][index.column()-3])
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
        if (index.column() == 3) or (index.column() == 4):
            if elem not in self.extra_abundances:
                self.extra_abundances[elem] = [None, None]
            if value == "":
                self.extra_abundances[elem][index.column()-3] = None
                return True
            try:
                value = float(value)
            except:
                return False
            else:
                self.extra_abundances[elem][index.column()-3] = value
                return True
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

    def get_extra_abundances(self):
        """
        Returns None if no extra abundances, otherwise two dictionaries of the extra abundances
        """
        if self.elem_order is None: return None
        extra_abundances_1 = {}
        extra_abundances_2 = {}
        for elem in self.elem_order.values():
            if elem not in self.extra_abundances:
                self.extra_abundances[elem] = [None, None]
                continue
            if self.extra_abundances[elem][0] is not None:
                extra_abundances_1[elem] = self.extra_abundances[elem][0]
            if self.extra_abundances[elem][1] is not None:
                extra_abundances_2[elem] = self.extra_abundances[elem][1]
        if (extra_abundances_1 == {}) and (extra_abundances_2 == {}):
            return None
        return extra_abundances_1, extra_abundances_2
