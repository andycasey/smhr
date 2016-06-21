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

from smh import utils
import mpl, style_utils
from matplotlib.ticker import MaxNLocator
from smh.photospheres.abundances import asplund_2009 as solar_composition
from smh.spectral_models import (ProfileFittingModel, SpectralSynthesisModel)
from spectral_models_table import SpectralModelsTableViewBase, SpectralModelsFilterProxyModel, SpectralModelsTableModelBase
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
        self.FeH = np.nan

        self.parent_layout = QtGui.QHBoxLayout(self)
        
        ################
        # LEFT HAND SIDE
        ################
        lhs_layout = QtGui.QVBoxLayout()
        
        hbox = QtGui.QHBoxLayout()
        self.filter_combo_box = QtGui.QComboBox(self)
        self.filter_combo_box.setSizeAdjustPolicy(QtGui.QComboBox.AdjustToContents)
        self.filter_combo_box.addItem("All")
        self.element_summary_text = QtGui.QLabel(self)
        self.element_summary_text.setText("Please load spectral models (or fit all)")
        sp = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, 
                               QtGui.QSizePolicy.Minimum)
        self.element_summary_text.setSizePolicy(sp)
        hbox.addWidget(self.filter_combo_box)
        hbox.addWidget(self.element_summary_text)
        lhs_layout.addLayout(hbox)

        self.table_view = SpectralModelsTableView(self)
        # Set up a proxymodel.
        self.proxy_spectral_models = SpectralModelsFilterProxyModel(self)
        self.proxy_spectral_models.add_filter_function(
            "use_for_stellar_composition_inference",
            lambda model: model.use_for_stellar_composition_inference)

        self.proxy_spectral_models.setDynamicSortFilter(True)
        header = ["", u"λ\n(Å)", "log ε\n(dex)", u"E. W.\n(mÅ)",
                  "REW", "σ(X)\n(dex)", "σ(E.W.)\n(mÅ)", "loggf","Element\n"]
        #attrs = ("is_acceptable", "_repr_wavelength", "abundance", "equivalent_width", 
        #         "reduced_equivalent_width", "_repr_element")
        self.all_spectral_models = SpectralModelsTableModel(self, header, None)
        self.proxy_spectral_models.setSourceModel(self.all_spectral_models)

        self.table_view.setModel(self.proxy_spectral_models)
        self.table_view.setSelectionBehavior(
            QtGui.QAbstractItemView.SelectRows)

        # TODO: Re-enable sorting.
        self.table_view.setSortingEnabled(False)
        self.table_view.resizeColumnsToContents()
        self.table_view.setColumnWidth(0, 30) # MAGIC
        self.table_view.setColumnWidth(1, 60) # MAGIC
        self.table_view.setColumnWidth(2, 50) # MAGIC
        self.table_view.setColumnWidth(3, 50) # MAGIC
        self.table_view.setColumnWidth(4, 50) # MAGIC
        self.table_view.setColumnWidth(5, 50) # MAGIC
        self.table_view.setColumnWidth(6, 50) # MAGIC
        self.table_view.setColumnWidth(7, 50) # MAGIC
        self.table_view.setMinimumSize(QtCore.QSize(240, 0))
        self.table_view.horizontalHeader().setStretchLastSection(True)
        lhs_layout.addWidget(self.table_view)

        # Buttons
        hbox = QtGui.QHBoxLayout()
        self.btn_fit_all = QtGui.QPushButton(self)
        self.btn_fit_all.setText("Fit all acceptable")
        self.btn_measure_all = QtGui.QPushButton(self)
        self.btn_measure_all.setText("Measure all acceptable")
        hbox.addWidget(self.btn_fit_all)
        hbox.addWidget(self.btn_measure_all)
        lhs_layout.addLayout(hbox)

        # Model fitting options
        self._create_fitting_options_widget()
        lhs_layout.addWidget(self.opt_tabs)
        
        self.parent_layout.addLayout(lhs_layout)

        #############################
        # RIGHT HAND SIDE: MPL WIDGET
        #############################
        rhs_layout = QtGui.QVBoxLayout()
        self.figure = mpl.MPLWidget(None, tight_layout=True, autofocus=True)
        self.figure.setMinimumSize(QtCore.QSize(300, 300))
        
        gs_top = matplotlib.gridspec.GridSpec(3,1,height_ratios=[1,2,1])
        gs_top.update(top=.95,bottom=.05,hspace=0)
        gs_bot = matplotlib.gridspec.GridSpec(3,1,height_ratios=[1,2,1])
        gs_bot.update(top=.95,bottom=.05,hspace=.3)
        
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
        
        self._points = [self.ax_line_strength.scatter([], [], s=30, \
             facecolor="k", edgecolor="k", picker=PICKER_TOLERANCE, \
             alpha=0.5)]
        self._trend_lines = None
        
        # Some empty figure objects that we will use later.
        self._lines = {
            "selected_point": [
                self.ax_line_strength.scatter([], [],
                    edgecolor="b", facecolor="none", s=150, linewidth=3, zorder=2)
            ],
            "spectrum": None,
            "spectrum_fill": None,
            "residual_fill": None,
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
            ],
            "dotted_line_at_one": self.ax_spectrum.plot([2000,10000],[1,1], 'k:')
        }
        
        rhs_layout.addWidget(self.figure)
        self.parent_layout.addLayout(rhs_layout)

        # Connect filter combo box
        self.filter_combo_box.currentIndexChanged.connect(self.filter_combo_box_changed)
        
        # Connect selection model
        _ = self.table_view.selectionModel()
        _.selectionChanged.connect(self.selected_model_changed)

        # Connect buttons
        self.btn_fit_all.clicked.connect(self.fit_all)
        self.btn_measure_all.clicked.connect(self.measure_all)

        # Connect matplotlib.
        self.figure.mpl_connect("button_press_event", self.figure_mouse_press)
        self.figure.mpl_connect("button_release_event", self.figure_mouse_release)
        self.figure.figure.canvas.callbacks.connect(
            "pick_event", self.figure_mouse_pick)
        # Zoom box and keyboard shortcuts
        self.figure.mpl_connect("button_press_event", self.figure.axis_right_mouse_press)
        self.figure.mpl_connect("button_release_event", self.figure.axis_right_mouse_release)
        self.figure.mpl_connect("key_press_event", self.figure.unzoom_on_z_press)
        self.figure.mpl_connect("key_press_event", self.key_press_zoom)
        # Check and uncheck
        self.figure.mpl_connect("key_press_event", self.key_press_check_uncheck)
        # Antimasks
        self.figure.mpl_connect("key_press_event", self.figure.key_press_flags)
        self.figure.mpl_connect("key_release_event", self.figure.key_release_flags)
        # Allow focusing figure for keyboard shortcuts
        self.figure.setFocusPolicy(QtCore.Qt.ClickFocus)
        
        self._currently_plotted_element = None
        self._rew_cache = []
        self._abund_cache = []
        self._err_cache = []
        self.refresh_table()

    def _create_fitting_options_widget(self):
        self.opt_tabs = QtGui.QTabWidget(self)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Expanding, 
            QtGui.QSizePolicy.MinimumExpanding)
        self.opt_tabs.setSizePolicy(sp)

        def _create_line_in_hbox(parent, text, bot, top, dec, validate_int=False):
            hbox = QtGui.QHBoxLayout()
            label = QtGui.QLabel(parent)
            label.setText(text)
            line = QtGui.QLineEdit(parent)
            line.setMinimumSize(QtCore.QSize(60, 0))
            line.setMaximumSize(QtCore.QSize(60, 16777215))
            if validate_int:
                line.setValidator(QtGui.QIntValidator(bot, top, line))
            else:
                line.setValidator(QtGui.QDoubleValidator(bot, top, dec, line))
            hbox.addWidget(label)
            hbox.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
                                           QtGui.QSizePolicy.Minimum))
            hbox.addWidget(line)
            return hbox, label, line

        def _create_checkline_in_hbox(parent, text, bot, top, dec):
            hbox = QtGui.QHBoxLayout()
            checkbox = QtGui.QCheckBox(parent)
            checkbox.setText("")
            label = QtGui.QLabel(parent)
            label.setText(text)
            line = QtGui.QLineEdit(parent)
            line.setMinimumSize(QtCore.QSize(60, 0))
            line.setMaximumSize(QtCore.QSize(60, 16777215))
            line.setValidator(QtGui.QDoubleValidator(bot, top, dec, line))
            hbox.addWidget(checkbox)
            hbox.addWidget(label)
            hbox.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
                                           QtGui.QSizePolicy.Minimum))
            hbox.addWidget(line)
            return hbox, checkbox, label, line

        def _create_combo_in_hbox(parent, text):
            hbox = QtGui.QHBoxLayout()
            label = QtGui.QLabel(parent)
            label.setText(text)
            combo = QtGui.QComboBox(parent)
            #combo.setMinimumSize(QtCore.QSize(60, 0))
            #combo.setMaximumSize(QtCore.QSize(60, 16777215))
            combo.setSizeAdjustPolicy(QtGui.QComboBox.AdjustToContents)
            hbox.addWidget(label)
            hbox.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
                                           QtGui.QSizePolicy.Minimum))
            hbox.addWidget(combo)
            return hbox, label, combo

        def _create_checkcombo_in_hbox(parent, text):
            hbox = QtGui.QHBoxLayout()
            checkbox = QtGui.QCheckBox(parent)
            checkbox.setText("")
            label = QtGui.QLabel(parent)
            label.setText(text)
            combo = QtGui.QComboBox(parent)
            combo.setMinimumSize(QtCore.QSize(60, 0))
            combo.setMaximumSize(QtCore.QSize(60, 16777215))
            hbox.addWidget(checkbox)
            hbox.addWidget(label)
            hbox.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
                                           QtGui.QSizePolicy.Minimum))
            hbox.addWidget(combo)
            return hbox, checkbox, label, combo

        ###################
        ### Profile options
        self.tab_profile = QtGui.QWidget()
        tab_hbox = QtGui.QHBoxLayout(self.tab_profile)

        ### LHS
        vbox_lhs = QtGui.QVBoxLayout()
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
        vbox_lhs.addWidget(self.checkbox_use_central_weighting)

        self.checkbox_use_antimasks = QtGui.QCheckBox(self.tab_profile)
        self.checkbox_use_antimasks.setText("Use Antimasks")
        self.checkbox_use_antimasks.setEnabled(False) # Editable by shift clicking only
        vbox_lhs.addWidget(self.checkbox_use_antimasks)

        ### RHS
        vbox_rhs = QtGui.QVBoxLayout()
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
        
        # Connect Signals for Profile
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
        self.btn_fit_one.clicked.connect(
            self.fit_one)
        self.btn_clear_masks.clicked.connect(
            self.clicked_btn_clear_masks)

        ###################
        ### Synthesis options
        self.tab_synthesis = QtGui.QWidget()
        tab_hbox = QtGui.QHBoxLayout(self.tab_synthesis)

        ### LHS
        vbox_lhs = QtGui.QVBoxLayout()
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

        # Element abundance table
        self.synth_abund_table = SynthesisAbundanceTableView(self.tab_synthesis)
        self.synth_abund_table_model = SynthesisAbundanceTableModel(self)
        self.synth_abund_table.setModel(self.synth_abund_table_model)
        self.synth_abund_table.resizeColumnsToContents()
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

        vbox_rhs.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Minimum,
                                           QtGui.QSizePolicy.Minimum))

        self.btn_fit_synth = QtGui.QPushButton(self.tab_synthesis)
        self.btn_fit_synth.setText("Fit Model")
        vbox_rhs.addWidget(self.btn_fit_synth)
        
        self.btn_update_abund_table = QtGui.QPushButton(self.tab_synthesis)
        self.btn_update_abund_table.setText("Update Abundance Table")
        vbox_rhs.addWidget(self.btn_update_abund_table)
        
        self.btn_clear_masks_2 = QtGui.QPushButton(self.tab_synthesis)
        self.btn_clear_masks_2.setText("Clear Masks")
        vbox_rhs.addWidget(self.btn_clear_masks_2)

        ### Finish Synthesis Tab
        tab_hbox.addLayout(vbox_lhs)
        tab_hbox.addLayout(vbox_rhs)
        self.opt_tabs.addTab(self.tab_synthesis, "Synthesis")

        # Connect signals for synthesis
        self.synthDefaultAction = self.fit_one

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
        self.btn_fit_synth.clicked.connect(
            self.fit_one)
        self.btn_update_abund_table.clicked.connect(
            self.clicked_btn_update_abund_table)
        self.btn_clear_masks_2.clicked.connect(
            self.clicked_btn_clear_masks)

    def new_session_loaded(self):
        """
        Call this whenever a new session is loaded
        TODO not tested
        """
        self.refresh_table()
        self.refresh_cache()
        self.summarize_current_table()
        self.refresh_plots()
        return None

    def populate_filter_combo_box(self):
        if self.parent.session is None: return None
        box = self.filter_combo_box
        box.clear()
        box.addItem("All")

        all_species = set([])
        for spectral_model in self.all_spectral_models.spectral_models:
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
        table_model = self.proxy_spectral_models
        table_model.delete_all_filter_functions()
        table_model.reset()
        if elem is None or elem == "" or elem == "All":
            self.element_summary_text.setText("")
        else:
            species = utils.element_to_species(elem)
            def filter_function(model):
                if isinstance(model, ProfileFittingModel):
                    return species in model.species
                elif isinstance(model, SpectralSynthesisModel):
                    return np.any([species in specie for specie in model.species])
            table_model.beginResetModel()
            table_model.add_filter_function(elem, filter_function)
            table_model.endResetModel()
        self.refresh_cache()
        self.summarize_current_table()
        self.refresh_plots()
        self.table_view.selectRow(0)
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
            N = self.proxy_spectral_models.rowCount()
            self.element_summary_text.setText("N={} lines".format(N))
            return None
        # Use cache to get abundance and avoid looping through proxy table
        ii = np.isfinite(self._abund_cache)
        N = np.sum(ii)
        _abund = self._abund_cache[ii]
        _errs = self._err_cache[ii]
        weights = 1/_errs**2
        total_weights = np.sum(weights)
        # TODO weights needed
        abund = np.mean(_abund)#np.sum(_abund*weights)/total_weights
        stdev = np.std(_abund)#np.sum(_errs*weights**2)/(total_weights**2)
        XH = abund - solar_composition(elem.split()[0])
        #if elem == "Fe I": self.FeH = XH
        self.calculate_FeH()
        XFe = XH - self.FeH
        text = "N={1} A({0})={2:5.2f} σ({0})={5:4.2f} [{0}/H]={3:5.2f} [{0}/Fe I]={4:5.2f}"
        text = text.format(elem,N,abund,XH,XFe,stdev)
        self.element_summary_text.setText(text)

        # TODO debug, checking against session.summarize_spectral_models()
        summary_dict = self.parent.session.summarize_spectral_models(organize_by_element=False)
        species = utils.element_to_species(elem)
        if species in summary_dict:
            print("From summary instead of cache",summary_dict[species])
        
        return None

    def refresh_table(self):
        if self.parent.session is None: return None
        self._check_for_spectral_models()
        self.proxy_spectral_models.reset()
        self.populate_filter_combo_box()
        # TODO
        #self.calculate_FeH()
        return None

    def refresh_plots(self):
        model, proxy_index, index = self._get_selected_model(True)
        start = time.time()
        self.update_spectrum_figure(redraw=False)
        self.update_selected_points_plot(redraw=False)
        self.update_line_strength_figure(redraw=True)
        #print("Time to refresh plots: {:.1f}s".format(time.time()-start))
        return None

    def fit_all(self):
        self._check_for_spectral_models()
        current_element_index = self.filter_combo_box.currentIndex()

        # Fit all acceptable
        num_unacceptable = 0
        for i,spectral_model in enumerate(self.all_spectral_models.spectral_models):
            if not spectral_model.is_acceptable:
                num_unacceptable += 1
                continue
            if isinstance(spectral_model, SpectralSynthesisModel):
                try:
                    res = spectral_model.fit()
                except (ValueError, RuntimeError, TypeError) as e:
                    logger.debug("Fitting error",spectral_model)
                    logger.debug(e)
            elif isinstance(spectral_model, ProfileFittingModel):
                try:
                    res = spectral_model.fit()
                except (ValueError, RuntimeError, TypeError) as e:
                    logger.debug("Fitting error",spectral_model)
                    logger.debug(e)
        # If none are acceptable, then fit all
        if num_unacceptable == self.all_spectral_models.rowCount(None):
            print("Found no acceptable spectral models, fitting all!")
            for i,spectral_model in enumerate(self.all_spectral_models.spectral_models):
                if isinstance(spectral_model, SpectralSynthesisModel):
                    try:
                        res = spectral_model.fit()
                    except (ValueError, RuntimeError, TypeError) as e:
                        logger.debug("Fitting error",spectral_model)
                        logger.debug(e)
                elif isinstance(spectral_model, ProfileFittingModel):
                    try:
                        res = spectral_model.fit()
                    except (ValueError, RuntimeError, TypeError) as e:
                        logger.debug("Fitting error",spectral_model)
                        logger.debug(e)

        self.proxy_spectral_models.reset()
        self.populate_filter_combo_box()
        self.refresh_cache()
        self.summarize_current_table()
        self.refresh_plots()
        self.filter_combo_box.setCurrentIndex(current_element_index)
        return None

    def measure_all(self):
        # Save this just to go back 
        current_element_index = self.filter_combo_box.currentIndex()

        # Gets abundances and uncertainties into session
        self.parent.session.measure_abundances()

        self.proxy_spectral_models.reset()
        self.populate_filter_combo_box()
        self.refresh_cache()
        self.summarize_current_table()
        self.refresh_plots()

        self.filter_combo_box.setCurrentIndex(current_element_index)
        return None

    def fit_one(self):
        spectral_model, proxy_index, index = self._get_selected_model(True)
        if spectral_model is None: return None
        try:
            res = spectral_model.fit()
        except (ValueError, RuntimeError) as e:
            logger.debug("Fitting error",spectral_model)
            logger.debug(e)
            return None
        self.table_view.update_row(proxy_index.row())
        self.update_cache(proxy_index)
        self.summarize_current_table()
        self.update_fitting_options()
        self.refresh_plots()
        return None

    def measure_one(self):
        spectral_model, proxy_index, index = self._get_selected_model(True)
        if spectral_model is None: return None

        self.parent.session.measure_abundances([spectral_model])

        self.table_view.update_row(proxy_index.row())
        self.update_cache(proxy_index)
        self.summarize_current_table()
        self.update_fitting_options()
        self.refresh_plots()
        return None

    def synthesize_current_model(self):
        spectral_model, proxy_index, index = self._get_selected_model(True)
        if spectral_model is None: return None
        spectral_model.update_fit_after_parameter_change()
        self.table_view.update_row(proxy_index.row())
        self.update_cache(proxy_index)
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
                dialog = TransitionsDialog(self.parent.session,
                    callbacks=[self.proxy_spectral_models.reset, 
                               self.refresh_table])
                dialog.exec_()

                # Do we even have any spectral models now?
                for sm in self.parent.session.metadata.get("spectral_models", []):
                    if sm.use_for_stellar_composition_inference: break
                else:
                    return False
            else:
                return False
        return True


    def figure_mouse_pick(self, event):
        """
        Trigger for when the mouse is used to select an item in the figure.

        :param event:
            The matplotlib event.
        """
        self.table_view.selectRow(event.ind[0])
        return None

    def figure_mouse_press(self, event):
        """
        Trigger for when the left mouse button is pressed in the figure.

        :param event:
            The matplotlib event.
        """
        print("figure_mouse_press",event)
        if event.button != 1: return None
        if event.inaxes in (self.ax_residual, self.ax_spectrum):
            self.spectrum_axis_mouse_press(event)
        return None


    def figure_mouse_release(self, event):
        """
        Trigger for when the left mouse button is released in the figure.

        :param event:
            The matplotlib event.
        """
        print("figure_mouse_release",event)
        if event.button != 1: return None
        if event.inaxes in (self.ax_residual, self.ax_spectrum):
            self.spectrum_axis_mouse_release(event)
        return None


    def key_press_zoom(self, event):
        print("key_press_zoom",event)
        if event.key not in "1234": return None
        if self.parent.session is None: return None
        ylim = self.parent.session.setting(["zoom_shortcuts",int(event.key)],
                                           default_return_value=[0.0,1.2])
        self.ax_spectrum.set_ylim(ylim)
        self.figure.draw()
        return None

    def key_press_check_uncheck(self, event):
        print("key_press_check_uncheck",event)
        if event.key not in ["u", "U", "a", "A"]: return None
        proxy_indices = self.table_view.selectionModel().selectedRows()
        if event.key in ["u", "U"]:
            print("Pressed",event.key,"marking as unacceptable")
            self.table_view.mark_selected_models_as_unacceptable(proxy_indices)
        elif event.key in ["a", "A"]:
            print("Pressed",event.key,"marking as acceptable")
            self.table_view.mark_selected_models_as_acceptable(proxy_indices)
        return None
            
    def spectrum_axis_mouse_press(self, event):
        """
        The mouse button was pressed in the spectrum axis.

        :param event:
            The matplotlib event.
        """

        if event.dblclick:

            # Double click.
            spectral_model, proxy_index, index = self._get_selected_model(True)

            for i, (s, e) in enumerate(spectral_model.metadata["mask"][::-1]):
                if e >= event.xdata >= s:
                    # Remove a mask
                    print("Removing mask")
                    # TODO this doesn't seem to work?
                    mask = spectral_model.metadata["mask"]
                    index = len(mask) - 1 - i
                    del mask[index]

                    # Re-fit the current spectral_model.
                    spectral_model.fit()

                    # Update the view for this row.
                    self.table_view.update_row(proxy_index.row())

                    # Update the view of the current model.
                    self.update_spectrum_figure(True)
                    break

            else:
                # No match with a masked region. 
                # TODO: Add a point that will be used for the continuum?
                # For the moment just refit the model.
                spectral_model.fit()

                # Update the view for this row.
                self.table_view.update_row(proxy_index.row())

                # Update the view of the current model.
                self.update_spectrum_figure(True)
                return None

        else:
            selected_model = self._get_selected_model()
            # HACK
            if "antimask_flag" not in selected_model.metadata:
                selected_model.metadata["antimask_flag"] = False
            # Clear all masks if shift key state is not same as antimask_flag
            # Also change the antimask state
            if selected_model.metadata["antimask_flag"] != self.figure.shift_key_pressed:
                selected_model.metadata["mask"] = []
                selected_model.metadata["antimask_flag"] = not selected_model.metadata["antimask_flag"]
                print("Switching antimask flag to",selected_model.metadata["antimask_flag"])
                # HACK
                self.update_fitting_options()

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
                patch.set_facecolor("g" if selected_model.metadata["antimask_flag"] else "r")

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
        print("update_mask_region",event)

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
            print("Fitting")
            spectral_model.fit()

            # Update the table view for this row.
            self.table_view.update_row(proxy_index.row())

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
        try:
            proxy_index = self.table_view.selectionModel().selectedRows()[-1]
        except IndexError:
            return (None, None, None) if full_output else None
        index = self.proxy_spectral_models.mapToSource(proxy_index).row()
        model = self.parent.session.metadata["spectral_models"][index]
        return (model, proxy_index, index) if full_output else model

    def selected_model_changed(self):
        self.update_fitting_options()
        self.refresh_plots()
        self.figure.reset_zoom_limits()
        return None

    def update_spectrum_figure(self, redraw=False):
        """
        TODO refactor all this plotting?
        Currently copied straight from stellar_parameters.py with minor changes
        """
        selected_model = self._get_selected_model()
        if selected_model is None:
            print("No model selected")
            return None
        transitions = selected_model.transitions
        try:
            window = float(self.edit_view_window.text())
        except:
            window = selected_model.metadata["window"]
        limits = [
            transitions["wavelength"][0] - window,
            transitions["wavelength"][-1] + window,
        ]

        if hasattr(self.parent, "session") \
        and hasattr(self.parent.session, "normalized_spectrum"):
            try:
                # Fix the memory leak!
                self.ax_spectrum.lines.remove(self._lines["spectrum"])
                self.ax_spectrum.collections.remove(self._lines["spectrum_fill"])
                self.ax_residual.collections.remove(self._lines["residual_fill"])
            except Exception as e:
                # TODO fail in a better way
                print(e)

            # Draw the spectrum.
            spectrum = self.parent.session.normalized_spectrum
            plot_ii = np.logical_and(spectrum.dispersion > limits[0]-10,
                                     spectrum.dispersion < limits[1]+10)
            if np.sum(plot_ii)==0: 
                # Can't plot, no points!
                return None
            self._lines["spectrum"] = self.ax_spectrum.plot(spectrum.dispersion[plot_ii],
                spectrum.flux[plot_ii], c="k", drawstyle="steps-mid")[0]

            sigma = 1.0/np.sqrt(spectrum.ivar[plot_ii])
            self._lines["spectrum_fill"] = \
            style_utils.fill_between_steps(self.ax_spectrum, spectrum.dispersion[plot_ii],
                spectrum.flux[plot_ii] - sigma, spectrum.flux[plot_ii] + sigma, 
                facecolor="#cccccc", edgecolor="#cccccc", alpha=1)

            self._lines["residual_fill"] = \
            style_utils.fill_between_steps(self.ax_residual, spectrum.dispersion[plot_ii],
                -sigma, +sigma, facecolor="#CCCCCC", edgecolor="none", alpha=1)

            self.ax_spectrum.set_ylim(0, 1.2)
            self.ax_spectrum.set_yticks([0, 0.5, 1])
            #self.ax_spectrum.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2])
            three_sigma = 3*np.median(sigma[np.isfinite(sigma)])
            self.ax_residual.set_ylim(-three_sigma, three_sigma)
        
        # Zoom to region.
        self.ax_spectrum.set_xlim(limits)
        self.ax_residual.set_xlim(limits)
        self.figure.reset_zoom_limits()
            
        # If this is a profile fitting line, show where the centroid is.
        x = transitions["wavelength"][0] \
            if isinstance(selected_model, ProfileFittingModel) else np.nan
        self._lines["transitions_center_main"].set_data([x, x], [0, 1.2])
        self._lines["transitions_center_residual"].set_data([x, x], [0, 1.2])
        # Model masks specified by the user.
        # (These should be shown regardless of whether there is a fit or not.)
        # HACK
        mask_color = "g" if "antimask_flag" in selected_model.metadata and \
            selected_model.metadata["antimask_flag"] else "r"
        for i, (start, end) in enumerate(selected_model.metadata["mask"]):
            try:
                patches = self._lines["model_masks"][i]

            except IndexError:
                self._lines["model_masks"].append([
                    self.ax_spectrum.axvspan(np.nan, np.nan,
                        facecolor=mask_color, edgecolor="none", alpha=0.25),
                    self.ax_residual.axvspan(np.nan, np.nan,
                        facecolor=mask_color, edgecolor="none", alpha=0.25)
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
                patch.set_facecolor(mask_color)
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

    def update_cache(self, proxy_index):
        """
        Update the point plotting cache at a single proxy table index.
        This is also used to compute the combo box summary.
        """
        if self.filter_combo_box.currentText() == "All":
            return None
        proxy_row = proxy_index.row()
        table_model = self.proxy_spectral_models
        try:
            if not table_model.data(table_model.createIndex(proxy_row, 0, None), QtCore.Qt.CheckStateRole):
                raise ValueError #to put in nan
            rew = float(table_model.data(table_model.createIndex(proxy_row, 4, None)))
            abund = float(table_model.data(table_model.createIndex(proxy_row, 2, None)))
            err = float(table_model.data(table_model.createIndex(proxy_row, 5, None)))
        except ValueError:
            self._rew_cache[proxy_row] = np.nan
            self._abund_cache[proxy_row] = np.nan
            self._err_cache[proxy_row] = np.nan
        else:
            self._rew_cache[proxy_row] = rew
            self._abund_cache[proxy_row] = abund
            self._err_cache[proxy_row] = err
        
    def refresh_cache(self):
        """
        Compute cache of REW and abundance from spectral model table model
        Wow this is the worst naming ever
        Note that we should use np.nan for REW for synthesis models to keep indices ok
        I believe that is correctly done in the table model
        """
        current_element =  self.filter_combo_box.currentText()
        self._currently_plotted_element = current_element
        if current_element == "All":
            print("Resetting cache for All")
            self._rew_cache = np.array([])
            self._abund_cache = np.array([])
            self._err_cache = np.array([])
            return None
        table_model = self.proxy_spectral_models
        rew_list = []
        abund_list = []
        err_list = []
        for row in range(table_model.rowCount()):
            try:
                if not table_model.data(table_model.createIndex(row, 0, None), QtCore.Qt.CheckStateRole):
                    raise ValueError #to put in nan
                rew = float(table_model.data(table_model.createIndex(row, 4, None)))
                abund = float(table_model.data(table_model.createIndex(row, 2, None)))
                err = float(table_model.data(table_model.createIndex(row, 5, None)))
            except ValueError:
                rew_list.append(np.nan); abund_list.append(np.nan); err_list.append(np.nan)
            else:
                rew_list.append(rew); abund_list.append(abund); err_list.append(err)
        self._rew_cache = np.array(rew_list)
        self._abund_cache = np.array(abund_list)
        self._err_cache = np.array(err_list)
        
    def update_selected_points_plot(self, redraw=False):
        """
        Plot selected points
        """
        if self.filter_combo_box.currentText() == "All":
            if redraw: self.figure.draw()
            return None
        # These are the proxy model indices
        indices = np.unique(np.array([index.row() for index in \
            self.table_view.selectionModel().selectedRows()]))
        if len(indices) == 0:
            self._lines["selected_point"][0].set_offsets(np.array([np.nan,np.nan]).T)
            if redraw: self.figure.draw()
            return None
        print("Update selected points plot: {} {} {}".format(indices, \
              self._rew_cache[indices],self._abund_cache[indices]))
        
        self._lines["selected_point"][0].set_offsets(\
            np.array([self._rew_cache[indices],self._abund_cache[indices]]).T)
        if redraw: self.figure.draw()
        return None

    def update_line_strength_figure(self, redraw=False, use_cache=True):
        current_element =  self.filter_combo_box.currentText()
        if current_element == "All":
            # This should remove all points
            self._points[0].set_offsets(np.array([self._rew_cache, self._abund_cache]).T)
            if redraw: self.figure.draw()
            return None
        # If new element or not using cache, refresh the cache
        if current_element != self._currently_plotted_element or not use_cache:
            self.refresh_cache()
            self.summarize_current_table()

        # Plot
        self._points[0].set_offsets(np.array([self._rew_cache, self._abund_cache]).T)
        style_utils.relim_axes(self.ax_line_strength)
        # TODO trend lines
        if redraw: self.figure.draw()
        return None

    def update_fitting_options(self):
        try:
            selected_model = self._get_selected_model()
        except IndexError:
            return None
        if selected_model is None: return None
    
        if isinstance(selected_model, ProfileFittingModel):
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
        else:
            self.opt_tabs.setTabEnabled(0, False)

        # Synthesis options.
        if isinstance(selected_model, SpectralSynthesisModel):
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
    def clicked_btn_update_abund_table(self):
        selected_model = self._get_selected_model()
        if selected_model is None: return None
        assert isinstance(selected_model, SpectralSynthesisModel), selected_model
        summary_dict = self.parent.session.summarize_spectral_models(organize_by_element=True)
        
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
            print("Run at least one fit before setting abundances of "
                  "fitted element {}!".format(elem))
        else:
            for i,elem in enumerate(selected_model.elements):
                abund = summary_dict[elem][1]
                key = "log_eps({})".format(elem)
                fitted_result[0][key] = abund
                fitted_result[2]["abundances"][i] = abund

        print(summary_dict)
        return None

class SpectralModelsTableView(SpectralModelsTableViewBase):
    def sizeHint(self):
        return QtCore.QSize(240,100)

    def minimumSizeHint(self):
        return QtCore.QSize(240,0)

    def refresh_gui(self):
        self.parent.summarize_current_table()
        self.parent.update_fitting_options()
        self.parent.refresh_plots()
        return None

    def fit_selected_models(self, proxy_indices):
        """ Fit the selected spectral models. """
        # Fit the models one by one
        for proxy_index in proxy_indices:
            index = self.model().mapToSource(proxy_index).row()
            self.parent.parent.session.metadata["spectral_models"][index].fit()

            # Update the data model, view, and cache
            self.update_row(proxy_index.row())
            self.parent.update_cache(proxy_index)

        self.refresh_gui()
        return None
    
    def measure_selected_models(self, proxy_indices):
        """ Fit the selected spectral models. """

        # Get list of selected spectral models
        spectral_models = []
        for proxy_index in proxy_indices:
            index = self.model().mapToSource(proxy_index).row()
            spectral_models.append(self.parent.parent.session.metadata["spectral_models"][index])

        # Fit abundances
        self.parent.parent.session.measure_abundances(spectral_models)

        # Update the data model and cache
        start = time.time()
        for proxy_index in proxy_indices:
            self.update_row(proxy_index.row())
            self.parent.update_cache(proxy_index)
        print("Time to update data model and cache: {:.1f}".format(time.time()-start))

        self.refresh_gui()
        return None

    def mark_selected_models_as_acceptable(self, proxy_indices):
        proxy_model = self.parent.proxy_spectral_models
        full_model = proxy_model.sourceModel()
        for proxy_index in proxy_indices:
            full_index = proxy_model.mapToSource(proxy_index)
            full_model.setData(full_index, 2, refresh_view=False)
        self.refresh_gui()
        return None
    def mark_selected_models_as_unacceptable(self, proxy_indices):
        proxy_model = self.parent.proxy_spectral_models
        full_model = proxy_model.sourceModel()
        for proxy_index in proxy_indices:
            full_index = proxy_model.mapToSource(proxy_index)
            full_model.setData(full_index, 0, refresh_view=False)
        self.refresh_gui()
        return None

    def set_fitting_option_value(self, proxy_indices, key, value,
                                 valid_for_profile=False,
                                 valid_for_synth=False):
        num_fit = 0
        num_unacceptable = 0
        num_profile_models = 0
        num_synthesis_models = 0
        for proxy_index in proxy_indices:
            idx = self.model().mapToSource(proxy_index).row()
            spectral_model \
                = self.parent.parent.session.metadata["spectral_models"][idx]
            run_fit = False
            if not spectral_model.is_acceptable: 
                num_unacceptable += 1
                continue
            if valid_for_profile and isinstance(spectral_model,ProfileFittingModel):
                num_profile_models += 1
                spectral_model.metadata[key] = value
                run_fit = True
            if valid_for_synth and isinstance(spectral_model,SpectralSynthesisModel):
                num_synthesis_models += 1
                spectral_model.metadata[key] = value
                run_fit = True
                
            if run_fit and "fitted_result" in spectral_model.metadata:
                num_fit += 1
                spectral_model.fit()
                self.update_row(proxy_index.row())
                self.parent.update_cache(proxy_index)
        print("Changed {0}={1}, fit {2} out of {3} models ({4} profile, {5} synth, skipped {6} unacceptable)".format(\
                key, value, num_fit, len(proxy_indices), num_profile_models, num_synthesis_models, num_unacceptable))
        self.refresh_gui()
        return None


class SpectralModelsTableModel(SpectralModelsTableModelBase):
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
        spectral_model = self.spectral_models[index.row()]

        if  column == 0 \
        and role in (QtCore.Qt.DisplayRole, QtCore.Qt.CheckStateRole):
            value = spectral_model.is_acceptable
            if role == QtCore.Qt.CheckStateRole:
                return QtCore.Qt.Checked if value else QtCore.Qt.Unchecked
            else:
                return None
        elif column == 1:
            value = spectral_model._repr_wavelength

        elif column == 2: #abundance
            try:
                abundances \
                    = spectral_model.metadata["fitted_result"][2]["abundances"]

            except (IndexError, KeyError):
                value = ""

            else:
                if len(abundances) == 1:
                    value = "{0:.2f}".format(abundances[0])
                else:
                    assert isinstance(spectral_model, SpectralSynthesisModel), spectral_model
                    current_element = self.parent.filter_combo_box.currentText().split()[0]
                    if current_element=="All":
                        try:
                            value = "; ".join(["{}".format(abund) \
                                               for abund in spectral_model.abundances])
                        except TypeError:
                            value = ""
                    else:
                        for i,elem in enumerate(spectral_model.elements):
                            if current_element == elem: break
                        else:
                            raise ValueError("{} not in {}".format(current_element, spectral_model.elements))
                        value = "{0:.2f}".format(abundances[i])
        elif column in [3, 4]: #EW, REW
            try:
                result = spectral_model.metadata["fitted_result"][2]
                equivalent_width = result["equivalent_width"][0]
            except:
                equivalent_width = np.nan

            if column == 3:
                value = "{0:.1f}".format(1000 * equivalent_width) \
                    if np.isfinite(equivalent_width) else ""
            if column == 4:
                value = "{:.2f}".format(np.log10(equivalent_width/float(spectral_model._repr_wavelength))) \
                    if np.isfinite(equivalent_width) else ""
        elif column == 5: #abundance err
            if isinstance(spectral_model, ProfileFittingModel):
                try:
                    result = spectral_model.metadata["fitted_result"][2]
                    err = result["abundance_uncertainties"][0]
                    value = "{:.2f}".format(err)
                except:
                    value = ""
            elif isinstance(spectral_model, SpectralSynthesisModel):
                current_element = self.parent.filter_combo_box.currentText().split()[0]
                if current_element=="All":
                    value = ""
                else:
                    for i,elem in enumerate(spectral_model.elements):
                        if current_element == elem: break
                    else:
                        raise ValueError("{} not in {}".format(current_element, spectral_model.elements))
                    try:
                        covar = spectral_model.metadata["fitted_result"][1]
                    except:
                        value = ""
                    else:
                        err = np.sqrt(covar[i,i])
                        value = "{:.2f}".format(err)
        elif column == 6: #EW err
            try:
                result = spectral_model.metadata["fitted_result"][2]
                err = 1000.*np.nanmax(np.abs(result["equivalent_width"][1:3]))
                value = "{:.2f}".format(err)
            except:
                value = ""
        elif column == 7:
            if isinstance(spectral_model, SpectralSynthesisModel):
                value = ""
            else:
                try:
                    loggf = spectral_model.transitions[0]['loggf']
                    value = "{:6.3f}".format(loggf)
                except:
                    value = ""
        elif column == 8:
            value = "; ".join(["{}".format(element) \
                      for element in spectral_model.elements])

        return value if role == QtCore.Qt.DisplayRole else None

    def setData(self, index, value, role=QtCore.Qt.DisplayRole, refresh_view=True):
        value = super(SpectralModelsTableModel, self).setData(index, value, role)
        if index.column() != 0: return False
        
        proxy_index = self.parent.table_view.model().mapFromSource(index)
        proxy_row = proxy_index.row()
        self.parent.table_view.rowMoved(proxy_row, proxy_row, proxy_row)

        self.parent.update_cache(proxy_index)
        if refresh_view:
            self.parent.summarize_current_table()
            self.parent.refresh_plots()

        return value

class SynthesisAbundanceTableView(QtGui.QTableView):
    """
    Make a small table view
    """
    def sizeHint(self):
        return QtCore.QSize(100,100)
    def minimumSizeHint(self):
        return QtCore.QSize(100,0)
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
            print(e)
            return 0
    def columnCount(self, parent):
        return 3
    def data(self, index, role):
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
                print("Run at least one fit before setting abundances!")
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
