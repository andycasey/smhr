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
from matplotlib.colors import ColorConverter
from matplotlib.ticker import MaxNLocator
from time import time

import mpl, style_utils
import astropy.table
import smh
from smh.gui.base import SMHSpecDisplay
from smh.gui.base import BaseTableView, MeasurementTableView, MeasurementTableDelegate
from smh.gui.base import MeasurementTableModelBase, MeasurementTableModelProxy
from smh.gui.base import SMHScatterplot
from smh.photospheres import available as available_photospheres
from smh.photospheres.abundances import asplund_2009 as solar_composition
from smh.spectral_models import (ProfileFittingModel, SpectralSynthesisModel)
from smh import utils
from linelist_manager import TransitionsDialog
from smh.optimize_stellar_params import optimize_stellar_parameters

from spectral_models_table import SpectralModelsTableViewBase, SpectralModelsFilterProxyModel, SpectralModelsTableModelBase
from quality_control import QualityControlDialog
from sp_solve_options import SolveOptionsDialog

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

class StellarParametersTab(QtGui.QWidget):


    def __init__(self, parent):
        """
        Create a tab for the determination of stellar parameters by excitation
        and ionization equalibrium.

        :param parent:
            The parent widget.
        """

        super(StellarParametersTab, self).__init__(parent)
        self.parent = parent

        sp = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, 
            QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)

        self.parent_layout = QtGui.QHBoxLayout(self)
        self.parent_layout.setContentsMargins(20, 20, 20, 0)

        ###########################################
        ############ Create LHS Layout ############
        ###########################################
        lhs_layout = QtGui.QVBoxLayout()
        lhs_layout.setSpacing(0)
        lhs_layout.setContentsMargins(0,0,0,0)
        
        ## Add RT options
        grid_layout = self._init_rt_options(parent)
        lhs_layout.addLayout(grid_layout)
        
        ## Add RT buttons
        hbox_layout = self._init_rt_buttons(parent)
        lhs_layout.addLayout(hbox_layout)

        ## Add state table (slopes)
        self._init_state_table(parent)
        lhs_layout.addWidget(self.state_table_view)

        ## Add measurement table
        self._init_measurement_table(parent)
        lhs_layout.addWidget(self.measurement_view)

        ## Add table buttons
        hbox_layout = self._init_other_buttons(parent)
        lhs_layout.addLayout(hbox_layout)
        self.parent_layout.addLayout(lhs_layout)
        ###########################################
        ############ Finish LHS Layout ############
        ###########################################

        ###########################################
        ############ Create RHS Layout ############
        ###########################################
        rhs_layout = QtGui.QVBoxLayout()
        rhs_layout.setSpacing(0)
        rhs_layout.setContentsMargins(0,0,0,0)
        self._init_scatterplots(parent)
        self._init_specfig(parent)
        rhs_layout.addWidget(self.expotfig)
        rhs_layout.addWidget(self.rewfig)
        rhs_layout.addWidget(self.specfig)
        self.parent_layout.addLayout(rhs_layout)
        ###########################################
        ############ Finish RHS Layout ############
        ###########################################

    def measure_abundances(self):
        """ 
        The measure abundances button has been clicked.
        - if no acceptable measurements, fit all
        - update session stellar parameters 
        - calculate abundances
        - update table and plots
        - update plots
        """
        if self.parent.session is None or not self._check_for_spectral_models():
            return None
        # If no acceptable measurements, fit all
        for model in self.parent.session.metadata["spectral_models"]:
             if model.use_for_stellar_parameter_inference \
            and "fitted_result" in model.metadata:
                break # do not fit if any stellar parameter lines are already fit
        else:
            for model in self.parent.session.metadata["spectral_models"]:
                if isinstance(model, SpectralSynthesisModel): continue
                try:
                    model.fit()
                except:
                    logger.exception(
                        "Exception in fitting spectral model {}".format(model))
                    continue
        self.update_stellar_parameter_session()
        
        ## Loop through the spectral model table and measure relevant abundances
        spectral_models = []
        for i,spectral_model in enumerate(self.full_measurement_model.spectral_models):
            if isinstance(spectral_model, ProfileFittingModel) and spectral_model.is_acceptable \
               and spectral_model.use_for_stellar_parameter_inference and (not spectral_model.is_upper_limit):
                spectral_models.append(spectral_model)
        self.parent.session.measure_abundances(spectral_models=spectral_models)
        self.measurement_model.reset()
        self.update_stellar_parameter_state_table()
        self.refresh_plots()
        return None
    
    def _init_rt_options(self, parent):
        grid_layout = QtGui.QGridLayout()
        # Effective temperature.
        label = QtGui.QLabel(self)
        label.setText("Teff")
        label.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum))
        grid_layout.addWidget(label, 0, 0, 1, 1)
        self.edit_teff = QtGui.QLineEdit(self)
        self.edit_teff.setMinimumSize(QtCore.QSize(40, 0))
        self.edit_teff.setMaximumSize(QtCore.QSize(50, 16777215))
        self.edit_teff.setAlignment(QtCore.Qt.AlignCenter)
        self.edit_teff.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum))
        self.edit_teff.setValidator(
            QtGui.QDoubleValidator(3000, 8000, 0, self.edit_teff))
        self.edit_teff.textChanged.connect(self._check_lineedit_state)
        grid_layout.addWidget(self.edit_teff, 0, 1)
        
        # Surface gravity.
        label = QtGui.QLabel(self)
        label.setText("logg")
        label.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum))

        grid_layout.addWidget(label, 1, 0, 1, 1)
        self.edit_logg = QtGui.QLineEdit(self)
        self.edit_logg.setMinimumSize(QtCore.QSize(40, 0))
        self.edit_logg.setMaximumSize(QtCore.QSize(50, 16777215))
        self.edit_logg.setAlignment(QtCore.Qt.AlignCenter)
        self.edit_logg.setValidator(
            QtGui.QDoubleValidator(-1, 6, 3, self.edit_logg))
        self.edit_logg.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum))
        self.edit_logg.textChanged.connect(self._check_lineedit_state)
        grid_layout.addWidget(self.edit_logg, 1, 1)

        # Metallicity.
        label = QtGui.QLabel(self)
        label.setText("[M/H]")
        label.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum))

        grid_layout.addWidget(label, 2, 0, 1, 1)
        self.edit_metallicity = QtGui.QLineEdit(self)
        self.edit_metallicity.setMinimumSize(QtCore.QSize(40, 0))
        self.edit_metallicity.setMaximumSize(QtCore.QSize(50, 16777215))
        self.edit_metallicity.setAlignment(QtCore.Qt.AlignCenter)
        self.edit_metallicity.setValidator(
            QtGui.QDoubleValidator(-5, 1, 3, self.edit_metallicity))
        self.edit_metallicity.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum))
        self.edit_metallicity.textChanged.connect(self._check_lineedit_state)
        grid_layout.addWidget(self.edit_metallicity, 2, 1)


        # Microturbulence.
        label = QtGui.QLabel(self)
        label.setText("vt")
        label.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum))

        grid_layout.addWidget(label, 3, 0, 1, 1)
        self.edit_xi = QtGui.QLineEdit(self)
        self.edit_xi.setMinimumSize(QtCore.QSize(40, 0))
        self.edit_xi.setMaximumSize(QtCore.QSize(50, 16777215))
        self.edit_xi.setAlignment(QtCore.Qt.AlignCenter)
        self.edit_xi.setValidator(QtGui.QDoubleValidator(0, 5, 3, self.edit_xi))
        self.edit_xi.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum))
        self.edit_xi.textChanged.connect(self._check_lineedit_state)
        grid_layout.addWidget(self.edit_xi, 3, 1)

        # Alpha-enhancement.
        label = QtGui.QLabel(self)
        label.setText("alpha")
        label.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum))
        
        grid_layout.addWidget(label, 4, 0, 1, 1)
        self.edit_alpha = QtGui.QLineEdit(self)
        self.edit_alpha.setMinimumSize(QtCore.QSize(40, 0))
        self.edit_alpha.setMaximumSize(QtCore.QSize(50, 16777215))
        self.edit_alpha.setAlignment(QtCore.Qt.AlignCenter)
        self.edit_alpha.setValidator(QtGui.QDoubleValidator(-1, 1, 3, self.edit_alpha))
        #self.edit_alpha.setValidator(QtGui.QDoubleValidator(0, 0.4, 3, self.edit_alpha))
        self.edit_alpha.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum))
        self.edit_alpha.textChanged.connect(self._check_lineedit_state)
        grid_layout.addWidget(self.edit_alpha, 4, 1)

        return grid_layout

    def new_session_loaded(self):
        if not hasattr(self.parent, "session") or self.parent.session is None:
            return None
        self.update_stellar_parameter_labels()
        self.full_measurement_model.new_session(session)
        self.measurement_model.reset()
        self.measurement_view.update_session(session)
        self.specfig.new_session(self.parent.session)
        self.update_stellar_parameter_state_table()
        self.refresh_plots()
        return None
    def update_stellar_parameter_session(self):
        """ Update the stellar parameters with the values in the GUI. """
        self.parent.session.metadata["stellar_parameters"].update({
            "effective_temperature": float(self.edit_teff.text()),
            "surface_gravity": float(self.edit_logg.text()),
            "metallicity": float(self.edit_metallicity.text()),
            "microturbulence": float(self.edit_xi.text()),
            "alpha": float(self.edit_alpha.text())
        })
        return True
    def update_stellar_parameter_labels(self):
        if not hasattr(self.parent, "session") or self.parent.session is None:
            return None
        widget_info = [
            (self.edit_teff, "{0:.0f}", "effective_temperature"),
            (self.edit_logg, "{0:.2f}", "surface_gravity"),
            (self.edit_metallicity, "{0:+.2f}", "metallicity"),
            (self.edit_xi, "{0:.2f}", "microturbulence"),
            (self.edit_alpha, "{0:.2f}", "alpha")
        ]
        for widget, fmt, key in widget_info:
            widget.setText(fmt.format(self.parent.session.metadata["stellar_parameters"][key]))
        return None
        
    def selected_measurement_changed(self):
        ta = time()
        try:
            selected_model = self._get_selected_model()
        except IndexError:
            self.update_selected_points(redraw=True)
            logger.debug("Time taken B: {}".format(time() - ta))
            return None
        if selected_model is None:
            logger.debug("No selected model: {}".format(time() - ta))
            return None
        logger.debug("selected model is at {}".format(selected_model._repr_wavelength))
        self.refresh_plots()
        logger.debug("Time taken: {}".format(time() - ta))
        return None
    def update_stellar_parameter_state_table(self):
        raise NotImplementedError
    def refresh_plots(self):
        self.expotfig.update_scatterplot(False)
        self.rewfig.update_scatterplot(False)
        self.expotfig.update_selected_points(True)
        self.rewfig.update_selected_points(True)
        self.specfig.update_spectrum_figure(True)
        return None
    def refresh_selected_points(self):
        self.expotfig.update_selected_points(True)
        self.rewfig.update_selected_points(True)
        return None
    def _get_selected_model(self, full_output=False):
        # Map the first selected row back to the source model index.
        try:
            proxy_index = self.measurement_view.selectionModel().selectedIndexes()[-1]
        except IndexError:
            return (None, None, None) if full_output else None
        index = self.measurement_model.mapToSource(proxy_index).row()
        model = self.parent.session.metadata["spectral_models"][index]
        return (model, proxy_index, index) if full_output else model

    def options(self):
        """ Open a GUI for the radiative transfer and solver options. """
        return SolveOptionsDialog(self.parent.session).exec_()
    def sperrors_dialog(self):
        """ Open a GUI for the stellar parameter errors. """
        return StellarParameterUncertaintiesDialog(self.parent.session).exec_()

    def solve_parameters(self):
        """ Solve the stellar parameters. """
        if self.parent.session is None or not self._check_for_spectral_models():
            return None

        ## use current state as initial guess
        logger.info("Setting [alpha/Fe]=0.4 to solve")
        self.update_stellar_parameters()
        sp = self.parent.session.metadata["stellar_parameters"]
        sp["alpha"] = 0.4
        initial_guess = [sp["effective_temperature"], sp["microturbulence"],
                         sp["surface_gravity"], sp["metallicity"]]

        ## grab transitions for stellar parameters
        transitions = []
        EWs = []
        for i, model in enumerate(self.parent.session.metadata["spectral_models"]):
            if model.use_for_stellar_parameter_inference and model.is_acceptable and not model.is_upper_limit:
                transitions.append(model.transitions[0])
                EWs.append(1e3 * model.metadata["fitted_result"][-1]["equivalent_width"][0])
        transitions = smh.LineList.vstack(transitions)
        transitions["equivalent_width"] = EWs
        
        ## TODO the optimization does not use the error weights yet
        ## TODO allow specification of tolerances and other params in QDialog widget
        out = optimize_stellar_parameters(initial_guess, transitions, \
                                              max_attempts=5, total_tolerance=1e-4, \
                                              individual_tolerances=None, \
                                              maxfev=30, use_nlte_grid=None)
        tolerance_achieved, initial_guess, num_moog_iterations, i, \
            t_elapsed, final_parameters, final_parameters_result, \
            all_sampled_points = out
        logger.info("Optimization took {:.1f}s".format(t_elapsed))
        if not tolerance_achieved:
            logger.warn("Did not converge in {}/{}!!! {:.5f} > {:.5f}".format( \
                    num_moog_iterations, i, 
                    np.sum(np.array(final_parameters_result)**2), 1e-4))
        
        ## update stellar params in metadata and gui
        sp["effective_temperature"] = final_parameters[0]
        sp["microturbulence"] = final_parameters[1]
        sp["surface_gravity"] = final_parameters[2]
        sp["metallicity"] = final_parameters[3]
        self.update_stellar_parameter_labels()
        ## refresh plot
        self.measure_abundances()


    def quality_control(self):
        """
        Show a dialog to specify quality control constraints for spectral models
        used in the determination of stellar parameters.
        """
        dialog = QualityControlDialog(self.parent.session,
            filter_spectral_models=lambda m: m.use_for_stellar_parameter_inference)
        dialog.exec_()
        if len(dialog.affected_indices) > 0:
            self.measurement_model.reset()
            self.refresh_plots()
        return None




    def _init_rt_buttons(self, parent):
        # Buttons for solving/measuring.        
        hbox = QtGui.QHBoxLayout()
        self.btn_measure = QtGui.QPushButton(self)
        self.btn_measure.setAutoDefault(True)
        self.btn_measure.setDefault(True)
        self.btn_measure.setText("Measure abundances")
        hbox.addWidget(self.btn_measure)

        self.btn_options = QtGui.QPushButton(self)
        self.btn_options.setText("Options..")
        hbox.addWidget(self.btn_options)

        self.btn_solve = QtGui.QPushButton(self)
        self.btn_solve.setText("Solve")
        hbox.addWidget(self.btn_solve)

        return hbox
    def _init_state_table(self, parent):
        self.state_table_view = QtGui.QTableView(self)
        self.state_table_view.setModel(StateTableModel(self))
        self.state_table_view.setSortingEnabled(False)
        self.state_table_view.verticalHeader().setResizeMode(QtGui.QHeaderView.Fixed)
        self.state_table_view.verticalHeader().setDefaultSectionSize(_ROWHEIGHT)
        self.state_table_view.setMaximumSize(QtCore.QSize(400, 3*(_ROWHEIGHT+1))) # MAGIC
        self.state_table_view.setSizePolicy(QtGui.QSizePolicy(
            QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.MinimumExpanding))
        self.state_table_view.setSelectionBehavior(
            QtGui.QAbstractItemView.SelectRows)

        self.state_table_view.horizontalHeader().setResizeMode(
            QtGui.QHeaderView.Stretch)

        self.state_table_view.horizontalHeader().setResizeMode(
            1, QtGui.QHeaderView.Fixed)
        self.state_table_view.horizontalHeader().resizeSection(1, 35) # MAGIC
        return None
    def _init_measurement_table(self, parent):
        self.full_measurement_model = MeasurementTableModelBase(self, self.parent.session, 
                                                    ["is_acceptable",
                                                     "wavelength","species","equivalent_width",
                                                     "equivalent_width_uncertainty",
                                                     "abundances","abundance_uncertainties",
                                                     "user_flag"])
        self.measurement_model = MeasurementTableModelProxy(self)
        self.measurement_model.setSourceModel(self.full_measurement_model)
        vbox, measurement_view, btn_filter, btn_refresh = base.create_measurement_table_with_buttons(
            self, self.measurement_model, self.parent.session,
            callbacks_after_menu=[self.new_session_loaded],
            display_fitting_options=False)
        self.measurement_view = measurement_view
        self.measurement_model.add_view_to_update(self.measurement_view)
        _ = self.measurement_view.selectionModel()
        _.selectionChanged.connect(self.selected_measurement_changed)
        self.measurement_view.setSizePolicy(QtGui.QSizePolicy(
            QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.MinimumExpanding))
        self.measurement_view.add_filter_function(
            "use_for_stellar_parameter_inference",
            lambda model: model.use_for_stellar_parameter_inference)
        
        self.btn_filter_acceptable = btn_filter
        self.btn_filter_acceptable.clicked.connect(self.refresh_plots)
        self.btn_refresh = btn_refresh
        self.btn_refresh.clicked.connect(self.new_session_loaded)
        
        return vbox
    def _init_other_buttons(self, parent):
        hbox = QtGui.QHBoxLayout()
        self.btn_quality_control = QtGui.QPushButton(self)
        self.btn_quality_control.setText("Quality control..")
        self.btn_sperrors = QtGui.QPushButton(self)
        self.btn_sperrors.setText("Stellar Parameter Uncertainties..")
        hbox.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Preferred,
            QtGui.QSizePolicy.Minimum))
        hbox.addWidget(self.btn_quality_control)
        return hbox
    def _init_scatterplots(self, parent):
        filters = [lambda x: (x.is_acceptable) and (x.species[0]==26.0) and (not x.is_upper_limit),
                   lambda x: (x.is_acceptable) and (x.species[0]==26.1) and (not x.is_upper_limit),
                   lambda x: (x.is_acceptable) and (x.species[0]==26.0) and (x.user_flag) and (not x.is_upper_limit),
                   lambda x: (x.is_acceptable) and (x.species[0]==26.1) and (x.user_flag) and (not x.is_upper_limit),
        ]
        point_styles = [{"s":40,"facecolor":"#FFFFFF","edgecolor":"k","linewidths":1,zorder=1},
                        {"s":40,"facecolor":"r","edgecolor":"k","linewidths":1,zorder=9},
                        {"s":70,"facecolor":"none","edgecolor":"red","linewidths":3,zorder=10},
                        {"s":70,"facecolor":"none","edgecolor":"red","linewidths":3,zorder=19},
        ]
        linefit_styles = [{"color":"k","linestyle":"--","zorder":-99},
                          {"color":"r","linestyle":"--","zorder":-99},
                          None,None]
        linemean_styles = [{"color":"k","linestyle":":","zorder":-999},
                           {"color":"r","linestyle":":","zorder":-999},
                           None,None]
        self.expotfig = SMHScatterplot(None, "expot", "abundances",
                                       tableview=self.measurement_view,
                                       filters=filters, point_styles=point_styles,
                                       linefit_styles=linefit_styles,linemean_styles=linemean_styles)
        self.rewfig = SMHScatterplot(None, "reduced_equivalent_width", "abundances",
                                     tableview=self.measurement_view,
                                     filters=filters, point_styles=point_styles,
                                     linefit_styles=linefit_styles,linemean_styles=linemean_styles)
        
        sp = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, 
                               QtGui.QSizePolicy.MinimumExpanding)
        for fig in [self.expotfig, self.rewfig]:
            fig.setSizePolicy(sp)
        return None
    def _init_specfig(self, parent):
        self.specfig = SMHSpecDisplay(None, self.parent.session, enable_masks=True,
                                      get_selected_model=self._get_selected_model)
        self.ax_spectrum = self.specfig.ax_spectrum
        self.ax_residual = self.specfig.ax_residual
    
    def _check_lineedit_state(self, *args, **kwargs):
        """
        Update the background color of a QLineEdit object based on whether the
        input is valid.
        """
        # TODO: Implement from
        # http://stackoverflow.com/questions/27159575/pyside-modifying-widget-colour-at-runtime-without-overwriting-stylesheet
        sender = self.sender()
        validator = sender.validator()
        state = validator.validate(sender.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            color = 'none' # normal background color
        elif state == QtGui.QValidator.Intermediate:
            color = '#fff79a' # yellow
        else:
            color = '#f6989d' # red
        sender.setStyleSheet('QLineEdit { background-color: %s }' % color)
        return None
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
                dialog = TransitionsDialog(self.parent.session,
                    callbacks=[
                        self.parent.transition_dialog_callback
                    ])
                dialog.exec_()

                # Do we even have any spectral models now?
                for sm in self.parent.session.metadata.get("spectral_models", []):
                    if sm.use_for_stellar_parameter_inference: break
                else:
                    return False
            else:
                return False
        return True




        
class StateTableModel(QtCore.QAbstractTableModel):
    header = ["Species", "N", u"〈[X/H]〉", u"σ", 
        u"∂A/∂χ", u"∂A/∂REW"]
    def __init__(self, parent, *args):
        super(StateTableModel, self).__init__(parent, *args)
        self.parent = parent
    def rowCount(self, parent):
        raise NotImplementedError
    def columnCount(self, parent):
        return len(self.header)
    def data(self, index, role):
        raise NotImplementedError
    def setData(self, index, value, role):
        return False
    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal \
        and role == QtCore.Qt.DisplayRole:
            return self.header[col]
        if role==QtCore.Qt.FontRole:
            return _QFONT
        return None
    def flags(self, index):
        if not index.isValid(): return
        return QtCore.Qt.ItemIsSelectable
    
class StellarParameterUncertaintiesDialog(QtGui.QDialog):
    def __init__(self, session,
                 default_tols=[5,0.01,0.01],
                 default_syserrs=[150,0.3,0.2],**kwargs):
        """
        A widget to calculate stellar parameter uncertainties.
        default_tols = [Teff, logg, vt] (5K, 0.01, 0.01)
        default_syserrs = [Teff, logg, vt] (150K, 0.3, 0.2)
        """
        super(StellarParameterUncertaintiesDialog, self).__init__(**kwargs)
        
        self.session = session
        Teff, logg, vt, MH = session.stellar_parameters
        stat_Teff, stat_logg, stat_vt, stat_MH = session.stellar_parameters_staterr
        sys_Teff, sys_logg, sys_vt, sys_MH = session.stellar_parameters_syserr
        tot_Teff, tot_logg, tot_vt, tot_MH = session.stellar_parameters_syserr
        
        # Display dialog in center and set size policy.
        self.setGeometry(320, 160, 320, 160)
        self.move(QtGui.QApplication.desktop().screen().rect().center() \
            - self.rect().center())
        self.setWindowTitle("Stellar parameter uncertainty analysis")
        
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, 
            QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)
        
        vbox = QtGui.QVBoxLayout(self)
        ### Create table grid layout
        grid = QtGui.QGridLayout()
        collabel, coltol, colstaterr, colsyserr, coltoterr = 0,1,2,3,4
        ## Column Labels
        grid.addWidget(QtGui.QLabel("Tolerance",self), 0, coltol)
        grid.addWidget(QtGui.QLabel("Stat Error",self), 0, colstaterr)
        grid.addWidget(QtGui.QLabel("Sys Error",self), 0, colsyserr)
        grid.addWidget(QtGui.QLabel("Total Error",self), 0, coltoterr)
        ## Row Labels
        self.label_Teff = QtGui.QLabel(self)
        self.label_logg = QtGui.QLabel(self)
        self.label_MH = QtGui.QLabel(self)
        self.label_vt = QtGui.QLabel(self)
        for i, label in enumerate([self.label_Teff, self.label_logg, self.label_MH, self.label_vt]):
            label.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum))
            grid.addWidget(label, i+1, collabel)
        ## Tolerances
        self.edit_tol_Teff = QtGui.QLineEdit(self)
        self.edit_tol_logg = QtGui.QLineEdit(self)
        self.edit_tol_vt = QtGui.QLineEdit(self)
        self.edit_tol_Teff.setText(str(default_tols[0]))
        self.edit_tol_logg.setText(str(default_tols[1]))
        self.edit_tol_vt.setText(str(default_tols[2]))
        for i, edit in zip([0,1,3],[self.edit_tol_Teff, self.edit_tol_logg, self.edit_tol_vt]):
            edit.setMinimumSize(QtCore.QSize(40, 0))
            edit.setMaximumSize(QtCore.QSize(50, 16777215))
            edit.setAlignment(QtCore.Qt.AlignCenter)
            edit.setValidator(QtGui.QDoubleValidator(0.001, 50, 3, edit))
            edit.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum))
            grid.addWidget(edit, i+1, coltol)
        ## StatErrs
        self.label_staterr_Teff = QtGui.QLabel(self)
        self.label_staterr_logg = QtGui.QLabel(self)
        self.label_staterr_MH = QtGui.QLabel(self)
        self.label_staterr_vt = QtGui.QLabel(self)
        for i, label in enumerate([self.label_staterr_Teff, self.label_staterr_logg,
                                   self.label_staterr_MH, self.label_staterr_vt]):
            label.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum))
            grid.addWidget(label, i+1, colstaterr)
        ## SysErrs
        self.edit_syserr_Teff = QtGui.QLineEdit(self)
        self.edit_syserr_logg = QtGui.QLineEdit(self)
        self.edit_syserr_MH = QtGui.QLineEdit(self)
        self.edit_syserr_vt = QtGui.QLineEdit(self)
        for i, edit in enumerate([self.edit_syserr_Teff, self.edit_syserr_logg,
                                  self.edit_syserr_MH, self.edit_syserr_vt]):
            edit.setMinimumSize(QtCore.QSize(40, 0))
            edit.setMaximumSize(QtCore.QSize(50, 16777215))
            edit.setAlignment(QtCore.Qt.AlignCenter)
            edit.setValidator(QtGui.QDoubleValidator(0.01, 1000, 2, edit))
            edit.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum))
            grid.addWidget(edit, i+1, colsyserr)
        ## TotErrs
        self.label_toterr_Teff = QtGui.QLabel(self)
        self.label_toterr_logg = QtGui.QLabel(self)
        self.label_toterr_MH = QtGui.QLabel(self)
        self.label_toterr_vt = QtGui.QLabel(self)
        for i, label in enumerate([self.label_toterr_Teff, self.label_toterr_logg,
                                   self.label_toterr_MH, self.label_toterr_vt]):
            label.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum))
            grid.addWidget(label, i+1, coltoterr)
        
        self.refresh_table()
        vbox.addLayout(grid)
        
        ### Add buttons to bottom
        hbox = QtGui.QHBoxLayout()
        self.btn_staterr = QtGui.QPushButton(self)
        self.btn_staterr.setText("Calculate Stat Err")
        self.btn_staterr.clicked.connect(self.calc_stat_err)
        hbox.addWidget(self.btn_staterr)
        self.btn_exit = QtGui.QPushButton(self)
        self.btn_exit.setText("Save to Session and Exit")
        self.btn_exit.setDefault(True)
        self.btn_exit.clicked.connect(self.save_and_exit)
        hbox.addWidget(self.btn_exit)
        vbox.addLayout(hbox)

    def refresh_table(self):
        Teff, logg, vt, MH = self.session.stellar_parameters
        stat_Teff, stat_logg, stat_vt, stat_MH = self.session.stellar_parameters_staterr
        sys_Teff, sys_logg, sys_vt, sys_MH = self.session.stellar_parameters_syserr
        tot_Teff, tot_logg, tot_vt, tot_MH = self.session.stellar_parameters_err
        
        self.label_Teff.setText("Teff={:.0f}".format(Teff))
        self.label_logg.setText("logg={:.2f}".format(logg))
        self.label_MH.setText("[M/H]={:.2f}".format(MH))
        self.label_vt.setText("vt={:.2f}".format(vt))
        
        self.label_staterr_Teff.setText("{:.0f}".format(stat_Teff))
        self.label_staterr_logg.setText("{:.2f}".format(stat_logg))
        self.label_staterr_MH.setText("{:.2f}".format(stat_MH))
        self.label_staterr_vt.setText("{:.2f}".format(stat_vt))
        
        self.edit_syserr_Teff.setText("{:.0f}".format(sys_Teff))
        self.edit_syserr_logg.setText("{:.2f}".format(sys_logg))
        self.edit_syserr_MH.setText("{:.2f}".format(sys_MH))
        self.edit_syserr_vt.setText("{:.2f}".format(sys_vt))
        
        self.label_toterr_Teff.setText("{:.0f}".format(tot_Teff))
        self.label_toterr_logg.setText("{:.2f}".format(tot_logg))
        self.label_toterr_MH.setText("{:.2f}".format(tot_MH))
        self.label_toterr_vt.setText("{:.2f}".format(tot_vt))
        
    def calc_stat_err(self):
        self.btn_staterr.setText("Calculating Stat Err...(may take a while)")
        tols = [float(x.text()) for x in [self.edit_tol_Teff, self.edit_tol_logg, self.edit_tol_vt]]
        syserrs = [float(x.text()) for x in [self.edit_syserr_Teff, self.edit_syserr_logg,
                                             self.edit_syserr_vt, self.edit_syserr_MH]]
        staterr_Teff, staterr_logg, staterr_vt, staterr_MH = self.session.stellar_parameter_uncertainty_analysis(
            tolerances=tols, systematic_errors=syserrs)
        self.label_staterr_Teff.setText("{:.0f}".format(staterr_Teff))
        self.label_staterr_logg.setText("{:.2f}".format(staterr_logg))
        self.label_staterr_MH.setText("{:.2f}".format(staterr_MH))
        self.label_staterr_vt.setText("{:.2f}".format(staterr_vt))
        self.btn_staterr.setText("Calculate Stat Err")
        self.save_to_session()
    def save_to_session(self):
        self.session.set_stellar_parameters_errors("stat",
            float(self.label_staterr_Teff.text()),
            float(self.label_staterr_logg.text()),
            float(self.label_staterr_vt.text()),
            float(self.label_staterr_MH.text()))
        self.session.set_stellar_parameters_errors("sys",
            float(self.edit_syserr_Teff.text()),
            float(self.edit_syserr_logg.text()),
            float(self.edit_syserr_vt.text()),
            float(self.edit_syserr_MH.text()))
        self.refresh_table()
    
    def save_and_exit(self):
        """
        Save to session metadata
        """
        self.save_to_session()
        self.close()
        
