#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The review tab in Spectroscopy Made Hard """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

#__all__ = ["ReviewTab"]

import logging
import numpy as np
import sys
from PySide2 import (QtCore, QtWidgets as QtGui)
import time

import smh
from smh.gui import base
from smh.gui.base import BaseTableView, MeasurementTableView, MeasurementTableDelegate
from smh.gui.base import MeasurementTableModelBase, MeasurementTableModelProxy
from smh.gui.base import MeasurementSummaryTableModel
from smh.gui.base import SMHScatterplot
from smh.gui import mpl, style_utils
from smh import utils

logger = logging.getLogger(__name__)
logger.addHandler(smh.handler)

class ReviewTab(QtGui.QWidget):
    def __init__(self, parent):
        super(ReviewTab, self).__init__(parent)
        self.parent = parent
        self.parent_layout = QtGui.QHBoxLayout(self)
        
        #######################################
        #### LHS: Summary and Measurement Table
        #######################################
        lhs_layout = QtGui.QVBoxLayout()
        self._init_summary_table()
        self._current_species = None
        lhs_layout.addWidget(self.summary_view)
        
        vbox = self._init_measurement_table()
        #lhs_layout.addWidget(self.measurement_view)
        lhs_layout.addLayout(vbox)
        self.parent_layout.addLayout(lhs_layout)
        
        ##################################
        #### RHS: Interactive Scatterplots
        ##################################
        rhs_layout = QtGui.QVBoxLayout()
        rhs_layout.setSpacing(0)
        self._init_scatterplots()
        rhs_layout.addWidget(self.plot1)
        rhs_layout.addWidget(self.plot2)
        rhs_layout.addWidget(self.plot3)
        self.parent_layout.addLayout(rhs_layout)
        return None
    
    def new_session_loaded(self):
        session = self.parent.session
        if session is None: return None
        self.measurement_model.beginResetModel()
        self.summary_model.new_session(session)
        self.full_measurement_model.new_session(session)
        self.measurement_model.reindex()
        self.measurement_model.endResetModel()
        self.measurement_view.update_session(session)
        self.refresh_plots()
        return None
        
    def selected_summary_changed(self):
        """
        Called when you select a new summary
        """
        start = time.time()
        # Get the species
        try:
            row = self.summary_view.selectionModel().selectedRows()[0].row()
        except IndexError:
            species = None
        else:
            species = self.summary_view.model().all_species[row]
        # Change measurement table filter
        logger.debug("species={}".format(species))
        self.change_measurement_table_species(species)
        # Change plots
        self.refresh_plots()
        logger.debug("selected_summary_changed to {}: {:.1f}".format(species, time.time()-start))
        return None
    def change_measurement_table_species(self, species):
        # Update the filter
        self.measurement_model.beginResetModel()
        if self._current_species is not None:
            _current_species_str = "{:.1f}".format(self._current_species)
            try:
                self.measurement_model.delete_filter_function(_current_species_str)
            except KeyError as e:
                logger.debug(e)
                logger.debug(self.measurement_model.filter_functions)
                raise
        if species is not None:
            filterfn = lambda model: _model_has_species(model, species)
            speciesstr = "{:.1f}".format(species)
            self.measurement_model.add_filter_function(speciesstr, filterfn)
        self._current_species = species
        # Reset the model (and its views)
        self.measurement_model.endResetModel()
    def refresh_plots(self):
        self.plot1.update_scatterplot(False)
        self.plot2.update_scatterplot(False)
        self.plot3.update_scatterplot(False)
        self.plot1.update_selected_points(True)
        self.plot2.update_selected_points(True)
        self.plot3.update_selected_points(True)
    def refresh_selected_points(self):
        self.plot1.update_selected_points(True)
        self.plot2.update_selected_points(True)
        self.plot3.update_selected_points(True)
    def selected_measurement_changed(self):
        self.refresh_selected_points()
        return None

    def _init_summary_table(self):
        # Create summary model
        self.summary_model = MeasurementSummaryTableModel(self, self.parent.session)
        # Create and link summary view
        self.summary_view = BaseTableView(self)
        self.summary_view.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        self.summary_view.setModel(self.summary_model)
        _ = self.summary_view.selectionModel()
        _.selectionChanged.connect(self.selected_summary_changed)

    def _init_measurement_table(self):
        self.full_measurement_model = MeasurementTableModelBase(self, self.parent.session, 
                                                    ["is_acceptable",
                                                     "wavelength","expot","loggf",
                                                     "equivalent_width","reduced_equivalent_width",
                                                     "abundances","abundances_to_solar",
                                                     "is_upper_limit","user_flag"])
        self.measurement_model = MeasurementTableModelProxy(self)
        self.measurement_model.setSourceModel(self.full_measurement_model)
        #self.measurement_view = MeasurementTableView(self)
        #self.measurement_view.setModel(self.measurement_model)
        vbox, measurement_view, btn_filter, btn_refresh = base.create_measurement_table_with_buttons(
            self, self.measurement_model, self.parent.session,
            callbacks_after_menu=[self.new_session_loaded],
            display_fitting_options=False)
        self.measurement_view = measurement_view
        self.measurement_model.add_view_to_update(self.measurement_view)
        _ = self.measurement_view.selectionModel()
        _.selectionChanged.connect(self.selected_measurement_changed)
        ## Sorting crashes :(
        #self.measurement_view.setSortingEnabled(True)
        #self.measurement_view.setItemDelegate(MeasurementTableDelegate(self,self.measurement_view))
        self.measurement_view.setSizePolicy(QtGui.QSizePolicy(
            QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.MinimumExpanding))

        self.btn_filter_acceptable = btn_filter
        self.btn_filter_acceptable.clicked.connect(self.refresh_plots)
        self.btn_refresh = btn_refresh
        self.btn_refresh.clicked.connect(self.new_session_loaded)

        return vbox

    def _init_scatterplots(self):
        filters = [lambda x: (x.is_acceptable) and (not x.user_flag) and (not x.is_upper_limit),
                   lambda x: not x.is_acceptable,
                   lambda x: x.user_flag,
                   lambda x: x.is_upper_limit]
        point_styles = [{"s":30,"facecolor":"k","edgecolor":"k","alpha":0.5},
                        {"s":30,"c":"c","marker":"x","linewidths":3},
                        {"s":30,"edgecolor":"r","facecolor":"k","alpha":0.5,"linewidths":3},
                        {"s":30,"edgecolor":"r","facecolor":"none","marker":"v"}
        ]
        linefit_styles = [{"color":"k","linestyle":"--","zorder":-99},None,None,None]
        linemean_styles = [{"color":"k","linestyle":":","zorder":-999},None,None,None]
        self.plot1 = SMHScatterplot(None, "expot", "abundances",
                                    tableview=self.measurement_view,
                                    filters=filters, point_styles=point_styles,
                                    linefit_styles=linefit_styles,linemean_styles=linemean_styles)
        self.plot2 = SMHScatterplot(None, "reduced_equivalent_width", "abundances",
                                    tableview=self.measurement_view,
                                    filters=filters, point_styles=point_styles,
                                    linefit_styles=linefit_styles,linemean_styles=linemean_styles)
        self.plot3 = SMHScatterplot(None, "wavelength", "abundances",
                                    tableview=self.measurement_view,
                                    filters=filters, point_styles=point_styles,
                                    linefit_styles=linefit_styles,linemean_styles=linemean_styles)
        
        sp = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, 
                               QtGui.QSizePolicy.MinimumExpanding)
        for fig in [self.plot1, self.plot2, self.plot3]:
            fig.setSizePolicy(sp)
        return None
        
def _model_has_species(model, species):
    for s in model.species:
        if isinstance(s, list) and species in s: return True
        elif isinstance(s, float) and species == s: return True
    return False
