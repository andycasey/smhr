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
from PyQt5 import (QtCore, QtGui as QtGui2, QtWidgets as QtGui)
from matplotlib.colors import ColorConverter
from matplotlib.ticker import MaxNLocator
from time import time

import mpl, style_utils
import astropy.table
import smh
from smh.gui.base import SMHSpecDisplay
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
        QtGui2.QFont.insertSubstitution(*substitute)


_QFONT = QtGui2.QFont("Helvetica Neue", 10)
_ROWHEIGHT = 20
DOUBLE_CLICK_INTERVAL = 0.1 # MAGIC HACK

class StateTableModel(QtCore.QAbstractTableModel):

    #header = ["Species", "N", u"〈[X/H]〉\n[dex]", u"σ\n[dex]", 
    #    u"∂A/∂χ\n[dex/eV]", u"∂A/∂REW\n[-]"]

    header = ["Species", "N", u"〈[X/H]", u"σ", 
        u"∂A/∂χ", u"∂A/∂REW"]


    def __init__(self, parent, *args):
        super(StateTableModel, self).__init__(parent, *args)
        self.parent = parent


    def rowCount(self, parent):
        try:
            state = self.parent._state_transitions
            finite_species = state["species"][np.isfinite(state["abundance"])]
            return len(set(finite_species))

        except AttributeError:
            return 0

    def columnCount(self, parent):
        return len(self.header)

    def data(self, index, role):
        if role==QtCore.Qt.FontRole:
            return _QFONT

        if not index.isValid() or role != QtCore.Qt.DisplayRole:
            return None

        try:
            state = self.parent._state_transitions

        except AttributeError:
            return None

        column = index.column()
        finite_abundances = np.isfinite(state["abundance"])
        finite_species = np.unique(state["species"][finite_abundances])

        if column == 0:
            return utils.species_to_element(finite_species[index.row()])

        elif column == 1:
            mask = finite_abundances \
                 * (state["species"] == finite_species[index.row()])
            return "{:.0f}".format(mask.sum())

        elif column == 2:
            species = finite_species[index.row()]
            mask = finite_abundances * (state["species"] == species)

            return "{0:.2f}".format(np.mean(state["abundance"][mask]) \
                - solar_composition(species))

        elif column == 3:
            mask = finite_abundances \
                 * (state["species"] == finite_species[index.row()])

            return "{0:.2f}".format(np.std(state["abundance"][mask]))

        elif column in (4, 5):

            key = ["expot", "reduced_equivalent_width"][column - 4]
            species = finite_species[index.row()]
            try:
                slope = self.parent._state_slopes[species][key][0]

            except (AttributeError, KeyError):
                return None

            else:
                return "{0:+.3f}".format(slope)

        return None


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
        lhs_layout.addWidget(self.table_view)

        _ = self.table_view.selectionModel()
        _.selectionChanged.connect(self.selected_model_changed)

        ## Add table buttons
        hbox_layout = self._init_table_buttons(parent)
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

        self._init_mpl_figure(parent)
        rhs_layout.addWidget(self.figure)
        rhs_layout.addWidget(self.specfig)
        self.parent_layout.addLayout(rhs_layout)

        # Some empty figure objects that we will use later.
        self._lines = {
            "excitation_slope_text": {},
            "line_strength_slope_text": {},
            "abundance_text": {},

            "excitation_trends": {},
            "line_strength_trends": {},
            "excitation_medians": {},
            "line_strength_medians": {},
            "scatter_points": {},
            "scatter_point_errors": {},
            "selected_point": [
                self.ax_excitation.scatter([], [],
                    edgecolor="b", facecolor="none", s=150, linewidth=3, zorder=1e4),
                self.ax_line_strength.scatter([], [],
                    edgecolor="b", facecolor="none", s=150, linewidth=3, zorder=1e4)
            ],
        }


        # Connect buttons.
        self.btn_measure.clicked.connect(self.measure_abundances)
        self.btn_options.clicked.connect(self.options)
        self.btn_solve.clicked.connect(self.solve_parameters)
        self.btn_filter.clicked.connect(self.filter_models)
        self.btn_quality_control.clicked.connect(self.quality_control)
        self.edit_teff.returnPressed.connect(self.measure_abundances)
        self.edit_logg.returnPressed.connect(self.measure_abundances)
        self.edit_metallicity.returnPressed.connect(self.measure_abundances)
        self.edit_xi.returnPressed.connect(self.measure_abundances)
        self.edit_alpha.returnPressed.connect(self.measure_abundances)

        # Connect matplotlib.
        self.figure.mpl_connect("button_press_event", self.figure_mouse_press)
        # Zoom box
        self.figure.enable_interactive_zoom()
        self.figure.setFocusPolicy(QtCore.Qt.ClickFocus)

        return None


    def populate_widgets(self):
        """ Update the stellar parameter edit boxes from the session. """

        if not hasattr(self.parent, "session") or self.parent.session is None:
            return None

        widget_info = [
            (self.edit_teff, "{0:.0f}", "effective_temperature"),
            (self.edit_logg, "{0:.2f}", "surface_gravity"),
            (self.edit_metallicity, "{0:+.2f}", "metallicity"),
            (self.edit_xi, "{0:.2f}", "microturbulence"),
            (self.edit_alpha, "{0:.2f}", "alpha")
        ]
        metadata = self.parent.session.metadata["stellar_parameters"]

        for widget, format, key in widget_info:
            widget.setText(format.format(metadata[key]))

        self.specfig.new_session(self.parent.session)

        return None


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
        if state == QtGui2.QValidator.Acceptable:
            color = 'none' # normal background color
        elif state == QtGui2.QValidator.Intermediate:
            color = '#fff79a' # yellow
        else:
            color = '#f6989d' # red
        sender.setStyleSheet('QLineEdit { background-color: %s }' % color)
    
        return None


    def filter_models(self):
        """
        Filter the view of the models used in the determination of stellar
        parameters. 
        """

        hide = self.btn_filter.text().startswith("Hide")

        if hide:
            self.proxy_spectral_models.add_filter_function(
                "is_acceptable", lambda model: model.is_acceptable)
        else:
            self.proxy_spectral_models.delete_filter_function("is_acceptable")

        text = "{} unacceptable models".format(("Hide", "Show")[hide])
        self.btn_filter.setText(text)
        return None


    def quality_control(self):
        """
        Show a dialog to specify quality control constraints for spectral models
        used in the determination of stellar parameters.
        """

        dialog = QualityControlDialog(self.parent.session,
            filter_spectral_models=lambda m: m.use_for_stellar_parameter_inference)
        dialog.exec_()

        # Update the state.
        if len(dialog.affected_indices) > 0:

            if hasattr(self, "_state_transitions"):
                indices = np.array(dialog.affected_indices)
                self._state_transitions["abundance"][indices] = np.nan
                self._state_transitions["reduced_equivalent_width"][indices] = np.nan

            # Update table and view.
            self.proxy_spectral_models.reset()
            self.update_scatter_plots()
            self.update_trend_lines(redraw=True)

        return None


    def figure_mouse_pick(self, event):
        """
        Trigger for when the mouse is used to select an item in the figure.

        :param event:
            The matplotlib event.
        """
        
        ycol = "abundance"
        xcol = {
            self.ax_excitation_twin: "expot",
            self.ax_line_strength_twin: "reduced_equivalent_width"
        }[event.inaxes]

        xscale = np.ptp(event.inaxes.get_xlim())
        yscale = np.ptp(event.inaxes.get_ylim())
        try:
            distance = np.sqrt(
                    ((self._state_transitions[ycol] - event.ydata)/yscale)**2 \
                +   ((self._state_transitions[xcol] - event.xdata)/xscale)**2)
        except AttributeError:
            # Stellar parameters have not been measured yet
            return None

        index = np.nanargmin(distance)

        # Because the state transitions are linked to the parent source model of
        # the table view, we will have to get the proxy index.
        proxy_index = self.table_view.model().mapFromSource(
            self.proxy_spectral_models.sourceModel().createIndex(index, 0)).row()

        self.table_view.selectRow(proxy_index)
        return None


    def figure_mouse_press(self, event):
        """
        Trigger for when the left mouse button is pressed in the figure.

        :param event:
            The matplotlib event.
        """

        if event.button != 1: return None

        if event.inaxes \
        in (self.ax_excitation, self.ax_excitation_twin,
            self.ax_line_strength, self.ax_line_strength_twin):
            self.figure_mouse_pick(event)

        return None


    def update_selected_scatter_point(self):
        """ The current-selected spectral model has been re-fit. """

        spectral_model, proxy_index, index = self._get_selected_model(True)

        # Update the reduced equivalent width for this model in the state.
        if hasattr(self, "_state_transitions"):

            meta = spectral_model.metadata["fitted_result"][-1]
            self._state_transitions["reduced_equivalent_width"][index] \
                = meta["reduced_equivalent_width"][0]
            self._state_transitions["abundance"][index] \
                = meta.get("abundances", [np.nan])[0]
            self._state_transitions["abundance_uncertainty"][index] \
                = meta.get("abundance_uncertainties", [np.nan])[0]

            self.update_scatter_plots()
            self.update_selected_points(redraw=True)

        return None


    def update_stellar_parameters(self):
        """ Update the stellar parameters with the values in the GUI. """

        self.parent.session.metadata["stellar_parameters"].update({
            "effective_temperature": float(self.edit_teff.text()),
            "surface_gravity": float(self.edit_logg.text()),
            "metallicity": float(self.edit_metallicity.text()),
            "microturbulence": float(self.edit_xi.text()),
            "alpha": float(self.edit_alpha.text())
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


    def update_scatter_plots(self, redraw=False):
        """
        Update the axes showing the abundances with respect to excitation
        potential and abundances with respect to reduced equivalent width.

        :param redraw: [optional]
            Force a redraw of the figure.
        """

        try:
            state = self._state_transitions

        except AttributeError:
            if redraw:
                self.figure.draw()
            return None


        # Group by species with finite abundances.
        finite_species = set(state["species"][np.isfinite(state["abundance"])])
        for species in finite_species:

            # We don't use setdefault because it would need to create the
            # object value before checking whether the key exists in the
            # dictionary, and we don't want that because matplotlib won't
            # clean it up.
            if species not in self._lines["scatter_points"]:
                zorder = self._zorders.get(species, 1)
                facecolor = self._colors.get(species, "#FFFFFF")

                self._lines["scatter_points"][species] = [
                    self.ax_excitation.scatter([], [], 
                        s=40, zorder=zorder, facecolor=facecolor,
                        linewidths=1),
                    self.ax_line_strength.scatter([], [], 
                        s=40, zorder=zorder, facecolor=facecolor,
                        linewidths=1)
                ]

                self._lines["scatter_point_errors"][species] = [
                    self.ax_excitation.errorbar(
                        np.nan * np.ones(2), np.nan * np.ones(2), 
                        yerr=np.nan * np.ones((2, 2)), fmt=None, 
                        ecolor="#666666", elinewidth=2, zorder=-10,
                        capsize=3),
                    self.ax_line_strength.errorbar(
                        np.nan * np.ones(2), np.nan * np.ones(2), 
                        xerr=np.nan * np.ones((2, 2)),
                        yerr=np.nan * np.ones((2, 2)), fmt=None,
                        ecolor="#666666", elinewidth=2,  zorder=-10,
                        capsize=3)
                ]

            ex_collection, ls_collection \
                = self._lines["scatter_points"][species]

            mask = state["species"] == species
            y = state["abundance"][mask]
            
            x = state["expot"][mask]
            ex_collection.set_offsets(np.array([x, y]).T)

            x = state["reduced_equivalent_width"][mask]
            ls_collection.set_offsets(np.array([x, y]).T)


            ex_container, ls_container \
                = self._lines["scatter_point_errors"][species]

            # No error in x-direction.
            _, (ex_yerr_top, ex_yerr_bot), (ybars, ) = ex_container

            x = state["expot"][mask]
            ex_collection.set_offsets(np.array([x, y]).T)

            # Update errors.
            yerr = state["abundance_uncertainty"][mask]
            yerr_top = y + yerr
            yerr_bot = y - yerr

            ex_yerr_top.set_xdata(x)
            ex_yerr_bot.set_xdata(x)
            ex_yerr_top.set_ydata(yerr_top)
            ex_yerr_bot.set_ydata(yerr_bot)

            ex_ysegments = [np.array([[xi, yt], [xi, yb]]) for xi, yt, yb in \
                zip(x, yerr_top, yerr_bot)]
            ybars.set_segments(ex_ysegments)


            _,  (ls_xerr_lt, ls_xerr_rt, ls_yerr_top, ls_yerr_bot), \
                (xbars, ybars) = ls_container
            
            # Update errors (y-errors the same as above).
            x = state["reduced_equivalent_width"][mask]
            xerr = np.nan * np.ones(mask.sum())
            # TODO: Show REW errors, even if they are small?
            xerr_rt = x + xerr
            xerr_lt = x - xerr

            ls_yerr_top.set_xdata(x)
            ls_yerr_bot.set_xdata(x)
            ls_yerr_top.set_ydata(yerr_top)
            ls_yerr_bot.set_ydata(yerr_bot)

            ls_xerr_lt.set_ydata(y)
            ls_xerr_rt.set_ydata(y)
            ls_xerr_lt.set_xdata(xerr_lt)
            ls_xerr_rt.set_xdata(xerr_rt)

            ybars.set_segments([np.array([[xi, yt], [xi, yb]]) \
                for xi, yt, yb in zip(x, yerr_top, yerr_bot)])

            xbars.set_segments([np.array([[xt, yi], [xb, yi]]) \
                for xt, xb, yi in zip(xerr_rt, xerr_lt, y)])


        # Update limits on the excitation and line strength figures.
        style_utils.relim_axes(self.ax_excitation)
        style_utils.relim_axes(self.ax_line_strength)

        self.ax_excitation_twin.set_ylim(self.ax_excitation.get_ylim())
        self.ax_line_strength_twin.set_ylim(self.ax_line_strength.get_ylim())
        self.ax_excitation_twin.set_ylabel(r"$\log_\epsilon({\rm X})$")
        self.ax_line_strength_twin.set_ylabel(r"$\log_\epsilon({\rm X})$")

        # Scale the left hand ticks to [X/H] or [X/M]

        # How many atomic number?
        Z = set([int(species) for species in finite_species])
        if len(Z) == 1:

            scaled_ticks = np.array(
                self.ax_excitation.get_yticks()) - solar_composition(Z)[0]

            self.ax_excitation.set_yticklabels(scaled_ticks)
            self.ax_line_strength.set_yticklabels(scaled_ticks)

            label = "[{}/H]".format(state["element"][mask][0].split()[0])
            self.ax_excitation.set_ylabel(label)
            self.ax_line_strength.set_ylabel(label)

        else:
            raise NotImplementedError
        
        # Update trend lines.
        self.update_trend_lines()

        if redraw:
            self.figure.draw()
        return None


    def measure_abundances(self):
        """ The measure abundances button has been clicked. """

        if self.parent.session is None or not self._check_for_spectral_models():
            return None

        # If no acceptable measurements, fit all
        for model in self.parent.session.metadata["spectral_models"]:
             if model.use_for_stellar_parameter_inference \
            and "fitted_result" in model.metadata:
                break # do not fit if any stellar parameter lines are already fit
        else:
            for model in self.parent.session.metadata["spectral_models"]:
                ## Only fit models for stellar parameter inference
                ##if not model.use_for_stellar_parameter_inference: continue
                # Actually just fit all profile models
                if isinstance(model, SpectralSynthesisModel): continue
                try:
                    model.fit()
                except:
                    logger.exception(
                        "Exception in fitting spectral model {}".format(model))
                    continue

        # Update the session with the stellar parameters in the GUI, and then
        # calculate abundances.
        self.update_stellar_parameters()

        filtering = lambda model: model.use_for_stellar_parameter_inference
        try:
            self._state_transitions, state, \
                = self.parent.session.stellar_parameter_state(full_output=True,
                    filtering=filtering)

        except ValueError as e:
            logger.warn(e)
            logger.warn("No measured transitions to calculate abundances for.")
            return None

        # The order of transitions may differ from the order in the table view.
        # We need to re-order the transitions by hashes.
        """
        print("STATE")
        print(self._state_transitions)

        print("MODELS")
        print([each.transition["wavelength"][0] for each in self.parent.session.metadata["spectral_models"]])

        # The number of transitions should match what is shown in the view.
        assert len(self._state_transitions) == self.table_view.model().rowCount(
            QtCore.QModelIndex())
        """

        # Otherwise we're fucked:
        expected_hashes = np.array([smh.LineList.hash(each.transitions[0]) for each in \
            self.parent.session.metadata["spectral_models"]]) 

        assert np.all(expected_hashes == self._state_transitions.compute_hashes())

        self.update_scatter_plots(redraw=True)

        # Draw trend lines based on the data already there.
        self.update_trend_lines()

        # Update selected entries.
        self.selected_model_changed()

        # Update abundance column for all rows.
        data_model = self.proxy_spectral_models.sourceModel()
        # 3 is the abundance column
        data_model.dataChanged.emit(
            data_model.createIndex(0, 3),
            data_model.createIndex(
                data_model.rowCount(QtCore.QModelIndex()), 3))

        # It ought to be enough just to emit the dataChanged signal, but
        # there is a bug when using proxy models where the data table is
        # updated but the view is not, so we do this hack to make it
        # work:
        self.table_view.columnMoved(3, 3, 3)

        self.proxy_spectral_models.reset()

        return None


    def update_trend_lines(self, redraw=False):
        """
        Update the trend lines in the figures.
        """

        
        if not hasattr(self, "_state_transitions"):
            if redraw:
                self.figure.draw()
            return None

        # Use abundance errors in fit?
        if self.parent.session.setting(("stellar_parameter_inference", 
            "use_abundance_uncertainties_in_line_fits"), True):
            yerr_column = "abundance_uncertainty"

        else:
            yerr_column = None

        self.state_table_view.model().beginResetModel()
        states = utils.equilibrium_state(self._state_transitions,
            columns=("expot", "reduced_equivalent_width"), ycolumn="abundance",
            yerr_column=yerr_column)

        self._state_slopes = states
        #self.state_table_view.model().reset()
        self.state_table_view.model().endResetModel()


        """
        # Offsets from the edge of axes.
        x_offset = 0.0125
        y_offset = 0.10
        y_space = 0.15
        """

        no_state = (np.nan, np.nan, np.nan, np.nan, 0)
        for i, (species, state) in enumerate(states.items()):
            if not state: continue

            color = self._colors[species]

            # Create defaults.
            if species not in self._lines["excitation_medians"]:
                self._lines["excitation_medians"][species] \
                    = self.ax_excitation.plot([], [], c=color, linestyle=":")[0]
            if species not in self._lines["line_strength_medians"]:
                self._lines["line_strength_medians"][species] \
                    = self.ax_line_strength.plot([], [], c=color, linestyle=":")[0]

            if species not in self._lines["excitation_trends"]:
                self._lines["excitation_trends"][species] \
                    = self.ax_excitation.plot([], [], c=color)[0]
            if species not in self._lines["line_strength_trends"]:
                self._lines["line_strength_trends"][species] \
                    = self.ax_line_strength.plot([], [], c=color)[0]

            # Do actual updates.
            #(m, b, np.median(y), np.std(y), len(x))
            m, b, median, sigma, N = state.get("expot", no_state)

            x = np.array(self.ax_excitation.get_xlim())
            self._lines["excitation_medians"][species].set_data(x, median)
            self._lines["excitation_trends"][species].set_data(x, m * x + b)

            m, b, median, sigma, N = state.get("reduced_equivalent_width", no_state)
            x = np.array(self.ax_line_strength.get_xlim())
            self._lines["line_strength_medians"][species].set_data(x, median)
            self._lines["line_strength_trends"][species].set_data(x, m * x + b)

            """
            # Show text.
            # TECH DEBT:
            # If a new species is added during stellar parameter determination
            # and some text is already shown, new text could appear on top of
            # that due to the way dictionaries (the state dictionary) is hashed.
            if species not in self._lines["abundance_text"]:
                self._lines["abundance_text"][species] \
                    = self.ax_excitation.text(
                        x_offset, 1 - y_offset - i * y_space, "",
                        color=color, transform=self.ax_excitation.transAxes,
                        horizontalalignment="left", verticalalignment="center")
            
            # Only show useful text.
            text = ""   if N == 0 \
                        else r"$\log_\epsilon{{\rm ({0})}} = {1:.2f} \pm {2:.2f}$"\
                             r" $(N = {3:.0f})$".format(
                                utils.species_to_element(species).replace(" ", "\,"),
                                median, sigma, N)
            self._lines["abundance_text"][species].set_text(text)


            m, b, median, sigma, N = state.get("expot", no_state)
            if species not in self._lines["excitation_slope_text"]:
                self._lines["excitation_slope_text"][species] \
                    = self.ax_excitation.text(
                        1 - x_offset, 1 - y_offset - i * y_space, "",
                        color=color, transform=self.ax_excitation.transAxes,
                        horizontalalignment="right", verticalalignment="center")

            # Only show useful text.
            text = ""   if not np.isfinite(m) \
                        else r"${0:+.3f}$ ${{\rm dex\,eV}}^{{-1}}$".format(m)
            self._lines["excitation_slope_text"][species].set_text(text)


            m, b, median, sigma, N = state.get("reduced_equivalent_width",
                no_state)
            if species not in self._lines["line_strength_slope_text"]:
                self._lines["line_strength_slope_text"][species] \
                    = self.ax_line_strength.text(
                        1 - x_offset, 1 - y_offset - i * y_space, "",
                        color=color, transform=self.ax_line_strength.transAxes,
                        horizontalalignment="right", verticalalignment="center")

            # Only show useful text.
            text = ""   if not np.isfinite(m) else r"${0:+.3f}$".format(m)
            self._lines["line_strength_slope_text"][species].set_text(text)
            """

        if redraw:
            self.figure.draw()

        return None



    def _get_selected_model(self, full_output=False):

        # Map the first selected row back to the source model index.
        try:
            proxy_index = self.table_view.selectionModel().selectedIndexes()[-1]
        except IndexError:
            return (None, None, None) if full_output else None
        index = self.proxy_spectral_models.mapToSource(proxy_index).row()
        model = self.parent.session.metadata["spectral_models"][index]
        return (model, proxy_index, index) if full_output else model


    def update_selected_points(self, redraw=False):
        # Show selected points.


        proxy_indices = np.unique(np.array([index for index in \
            self.table_view.selectionModel().selectedIndexes()]))

        if 1 > proxy_indices.size:
            for collection in self._lines["selected_point"]:
                collection.set_offsets(np.array([np.nan, np.nan]).T)
            if redraw:
                self.figure.draw()
            return None


        # These indices are proxy indices, which must be mapped back.
        indices = np.unique([self.table_view.model().mapToSource(index).row() \
            for index in proxy_indices])

        try:
            x_excitation = self._state_transitions["expot"][indices]
            x_strength = self._state_transitions["reduced_equivalent_width"][indices]
            y = self._state_transitions["abundance"][indices]

        except:
            x_excitation, x_strength, y = (np.nan, np.nan, np.nan)

        point_excitation, point_strength = self._lines["selected_point"]
        point_excitation.set_offsets(np.array([x_excitation, y]).T)
        point_strength.set_offsets(np.array([x_strength, y]).T)

        if redraw:
            self.figure.draw()

        return None


    def selected_model_changed(self):
        """
        The selected model was changed.
        """
        ta = time()

        # Show point on excitation/line strength plot.
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
        
        self.update_selected_points()

        # Show spectrum.
        self.update_spectrum_figure(redraw=True)
        self.figure.reset_zoom_limits()

        logger.debug("Time taken: {}".format(time() - ta))
        return None


    def update_spectrum_figure(self, redraw=True):
        """ Update the spectrum figure. """
        self.specfig.update_spectrum_figure(redraw=redraw)


    def options(self):
        """ Open a GUI for the radiative transfer and solver options. """

        return SolveOptionsDialog(self.parent.session).exec_()


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
        self.populate_widgets()
        ## refresh plot
        self.measure_abundances()


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
            QtGui2.QDoubleValidator(3000, 8000, 0, self.edit_teff))
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
            QtGui2.QDoubleValidator(-1, 6, 3, self.edit_logg))
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
            QtGui2.QDoubleValidator(-5, 1, 3, self.edit_metallicity))
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
        self.edit_xi.setValidator(QtGui2.QDoubleValidator(0, 5, 3, self.edit_xi))
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
        self.edit_alpha.setValidator(QtGui2.QDoubleValidator(-1, 1, 3, self.edit_alpha))
        #self.edit_alpha.setValidator(QtGui.QDoubleValidator(0, 0.4, 3, self.edit_alpha))
        self.edit_alpha.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum))

        self.edit_alpha.textChanged.connect(self._check_lineedit_state)
        grid_layout.addWidget(self.edit_alpha, 4, 1)

        return grid_layout

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
        self.state_table_view.verticalHeader().setSectionResizeMode(QtGui.QHeaderView.Fixed)
        self.state_table_view.verticalHeader().setDefaultSectionSize(_ROWHEIGHT)
        self.state_table_view.setMaximumSize(QtCore.QSize(400, 3*(_ROWHEIGHT+1))) # MAGIC
        self.state_table_view.setSizePolicy(QtGui.QSizePolicy(
            QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.MinimumExpanding))
        self.state_table_view.setSelectionBehavior(
            QtGui.QAbstractItemView.SelectRows)

        self.state_table_view.horizontalHeader().setSectionResizeMode(
            QtGui.QHeaderView.Stretch)

        self.state_table_view.horizontalHeader().setSectionResizeMode(
            1, QtGui.QHeaderView.Fixed)
        self.state_table_view.horizontalHeader().resizeSection(1, 35) # MAGIC
        return None

    def _init_measurement_table(self, parent):
        #header = ["", u"λ\n[Å]", "Element", u"EW\n[mÅ]", u"σ(EW)\n[mÅ]",
        #          "log ε\n[dex]", "σ(log ε)\n[dex]"]
        header = ["", u"λ", "Element", u"EW", u"σ(EW)",
                  "log ε", "σ(log ε)", "ul"]
        attrs = ("is_acceptable", "_repr_wavelength", "_repr_element", 
                 "equivalent_width", "err_equivalent_width", "abundance", "err_abundance",
                 "is_upper_limit")

        self.table_view = SpectralModelsTableView(self)
        self.table_view.verticalHeader().setSectionResizeMode(QtGui.QHeaderView.Fixed)
        self.table_view.verticalHeader().setDefaultSectionSize(_ROWHEIGHT)
        
        # Set up a proxymodel.
        self.proxy_spectral_models = SpectralModelsFilterProxyModel(self)
        self.proxy_spectral_models.add_filter_function(
            "use_for_stellar_parameter_inference",
            lambda model: model.use_for_stellar_parameter_inference)

        self.proxy_spectral_models.setDynamicSortFilter(True)
        self.proxy_spectral_models.setSourceModel(SpectralModelsTableModel(self, header, attrs))

        self.table_view.setModel(self.proxy_spectral_models)
        self.table_view.setSelectionBehavior(
            QtGui.QAbstractItemView.SelectRows)

        # TODO: Re-enable sorting.
        self.table_view.setSortingEnabled(False)
        self.table_view.setMaximumSize(QtCore.QSize(400, 16777215))        
        self.table_view.setSizePolicy(QtGui.QSizePolicy(
            QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.MinimumExpanding))
        
        # Keep the first colum to a fixed width, but the rest stretched.
        self.table_view.horizontalHeader().setSectionResizeMode(
            0, QtGui.QHeaderView.Fixed)
        self.table_view.horizontalHeader().resizeSection(0, 30) # MAGIC

        for i in range(1, len(header)):
            self.table_view.horizontalHeader().setSectionResizeMode(
                i, QtGui.QHeaderView.Stretch)
            
        return None

    def _init_table_buttons(self, parent):
        hbox = QtGui.QHBoxLayout()
        self.btn_filter = QtGui.QPushButton(self)
        self.btn_filter.setText("Hide unacceptable models")
        self.btn_quality_control = QtGui.QPushButton(self)
        self.btn_quality_control.setText("Quality control..")
        hbox.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Preferred,
            QtGui.QSizePolicy.Minimum))
        hbox.addWidget(self.btn_filter)
        hbox.addWidget(self.btn_quality_control)
        return hbox

    def _init_mpl_figure(self, parent):
        # Matplotlib figure.
        self.figure = mpl.MPLWidget(None, tight_layout=True, autofocus=True)
        self.figure.setMinimumSize(QtCore.QSize(10, 10))
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Expanding)
        #sp.setHorizontalStretch(0)
        #sp.setVerticalStretch(0)
        #sp.setHeightForWidth(self.figure.sizePolicy().hasHeightForWidth())
        self.figure.setSizePolicy(sp)
        #self.figure.setFocusPolicy(QtCore.Qt.StrongFocus)

        gs = matplotlib.gridspec.GridSpec(2, 1)
        gs.update(hspace=0.40)

        self._colors = {
            26.0: "#666666",
            26.1: "r"
        }
        self._zorders = {
            26.1: 10,
        }

        self.ax_excitation = self.figure.figure.add_subplot(gs[0])
        self.ax_excitation.xaxis.get_major_formatter().set_useOffset(False)
        self.ax_excitation.yaxis.set_major_locator(MaxNLocator(4))
        self.ax_excitation.set_xlabel(u"Excitation potential, χ (eV)")
        self.ax_excitation.set_ylabel("[X/H]")

        self.ax_excitation_twin = self.ax_excitation.twinx()
        self.ax_excitation_twin.yaxis.set_major_locator(MaxNLocator(4))
        self.ax_excitation_twin.set_ylabel(r"$\log_\epsilon({\rm X})$")

        self.ax_line_strength = self.figure.figure.add_subplot(gs[1])
        self.ax_line_strength.xaxis.get_major_formatter().set_useOffset(False)
        self.ax_line_strength.yaxis.set_major_locator(MaxNLocator(4))
        self.ax_line_strength.set_xlabel(
            r"Reduced equivalent width (REW), $\log({\rm EW}/\lambda)$")
        self.ax_line_strength.set_ylabel("[X/H]")

        self.ax_line_strength_twin = self.ax_line_strength.twinx()
        self.ax_line_strength_twin.yaxis.set_major_locator(MaxNLocator(4))
        self.ax_line_strength_twin.set_ylabel(r"$\log_\epsilon({\rm X})$")

        self.specfig = SMHSpecDisplay(None, self.parent.session, enable_masks=True,
                                      get_selected_model=self._get_selected_model)
        self.ax_spectrum = self.specfig.ax_spectrum
        self.ax_residual = self.specfig.ax_residual


class SpectralModelsTableView(SpectralModelsTableViewBase):

    def contextMenuEvent(self, event):
        """
        Provide a context (right-click) menu for the table containing the
        spectral models to use for stellar parameter inference.

        :param event:
            The mouse event that triggered the menu.
        """

        proxy_indices = self.selectionModel().selectedRows()
        indices = np.unique([self.model().mapToSource(index).row() \
            for index in proxy_indices])


        N = len(indices)

        menu = QtGui.QMenu(self)
        fit_models = menu.addAction(
            "Fit selected model{}..".format(["", "s"][N != 1]))
        menu.addSeparator()
        mark_as_acceptable = menu.addAction("Mark as acceptable")
        mark_as_unacceptable = menu.addAction("Mark as unacceptable")
        menu.addSeparator()

        # Common options.
        set_fitting_window = menu.addAction("Set fitting window..")

        continuum_menu = menu.addMenu("Set continuum")
        set_no_continuum = continuum_menu.addAction("No continuum",
            checkable=True)
        continuum_menu.addSeparator()
        set_continuum_order = [continuum_menu.addAction(
            "Order {}".format(i), checkable=True) for i in range(0, 10)]

        menu.addSeparator()
        menu_profile_type = menu.addMenu("Set profile type")

        set_gaussian = menu_profile_type.addAction("Gaussian")
        set_lorentzian = menu_profile_type.addAction("Lorentzian")
        set_voigt = menu_profile_type.addAction("Voigt")

        menu.addSeparator()

        enable_central_weighting = menu.addAction("Enable central weighting")
        disable_central_weighting = menu.addAction("Disable central weighting")

        menu.addSeparator()

        set_detection_sigma = menu.addAction("Set detection sigma..")
        set_detection_pixels = menu.addAction("Set detection pixels..")

        menu.addSeparator()

        set_rv_tolerance = menu.addAction("Set RV tolerance..")
        set_wl_tolerance = menu.addAction("Set WL tolerance..")

        if N == 0:
            menu.setEnabled(False)

        action = menu.exec_(self.mapToGlobal(event.pos()))

        update_spectrum_figure = False

        if action == fit_models:
            for idx in indices:
                spectral_model \
                    = self.parent.parent.session.metadata["spectral_models"][idx]
                spectral_model.fit()

            update_spectrum_figure = True


        elif action == mark_as_unacceptable:
            for idx in indices:
                spectral_model \
                    = self.parent.parent.session.metadata["spectral_models"][idx]
                spectral_model.metadata["is_acceptable"] = False

            # Update trend lines and scatter points.
            if len(indices) > 0:
                self.update_scatter_plots()
                self.update_selected_scatter_point()
                self.update_trend_lines(redraw=True)


        elif action == mark_as_acceptable:
            for idx in indices:
                spectral_model \
                    = self.parent.parent.session.metadata["spectral_models"][idx]
                if "fitted_result" in spectral_model.metadata:
                    spectral_model.metadata["is_acceptable"] = True

            # Update trend lines and scatter points.
            if len(indices) > 0:
                self.update_scatter_plots()
                self.update_selected_scatter_point()
                self.update_trend_lines(redraw=True)


        elif action == set_fitting_window:

            first_spectral_model = \
                self.parent.parent.session.metadata["spectral_models"][indices[0]]

            window, is_ok = QtGui.QInputDialog.getDouble(
                None, "Set fitting window", u"Fitting window (Å):", 
                value=first_spectral_model.metadata["window"],
                minValue=0.1, maxValue=1000)
            if not is_ok: return None

            for idx, proxy_index in zip(indices, proxy_indices):
                spectral_model \
                    = self.parent.parent.session.metadata["spectral_models"][idx]
                spectral_model.metadata["window"] = window

                if "fitted_result" in spectral_model.metadata:
                    spectral_model.fit()
                    update_spectrum_figure = True
                    self.parent.table_view.rowMoved(proxy_index.row(), 
                        proxy_index.row(), proxy_index.row())


        elif action == set_no_continuum:
            for idx, proxy_index in zip(indices, proxy_indices):
                spectral_model \
                    = self.parent.parent.session.metadata["spectral_models"][idx]
                spectral_model.metadata["continuum_order"] = -1

                # Re-fit if it already had a result.
                if "fitted_result" in spectral_model.metadata:
                    spectral_model.fit()
                    update_spectrum_figure = True
                    self.parent.table_view.rowMoved(proxy_index.row(), 
                        proxy_index.row(), proxy_index.row())


        elif action in set_continuum_order:            
            order = set_continuum_order.index(action)
            for idx, proxy_index in zip(indices, proxy_indices):
                spectral_model \
                    = self.parent.parent.session.metadata["spectral_models"][idx]
                spectral_model.metadata["continuum_order"] = order

                # Re-fit if it already had a result.
                if "fitted_result" in spectral_model.metadata:
                    spectral_model.fit()
                    update_spectrum_figure = True
                    self.parent.table_view.rowMoved(proxy_index.row(), 
                        proxy_index.row(), proxy_index.row())


        elif action in (set_gaussian, set_lorentzian, set_voigt):
            kind = {
                set_gaussian: "gaussian",
                set_lorentzian: "lorentzian",
                set_voigt: "voigt"
            }[action]

            for idx, proxy_index in zip(indices, proxy_indices):
                spectral_model \
                    = self.parent.parent.session.metadata["spectral_models"][idx]
                if isinstance(spectral_model, ProfileFittingModel):
                    spectral_model.metadata["profile"] = kind

                    if "fitted_result" in spectral_model.metadata:
                        spectral_model.fit()
                        update_spectrum_figure = True
                        self.parent.table_view.rowMoved(proxy_index.row(), 
                            proxy_index.row(), proxy_index.row())


        elif action in (enable_central_weighting, disable_central_weighting):
            toggle = (action == enable_central_weighting)
            for idx, proxy_index in zip(indices, proxy_indices):
                spectral_model \
                    = self.parent.parent.session.metadata["spectral_models"][idx]

                if isinstance(spectral_model, ProfileFittingModel):
                    spectral_model.metadata["central_weighting"] = toggle

                    if "fitted_result" in spectral_model.metadata:
                        spectral_model.fit()
                        update_spectrum_figure = True
                        self.parent.table_view.rowMoved(proxy_index.row(), 
                            proxy_index.row(), proxy_index.row())


        elif action == set_detection_sigma:

            # Get the first profile model
            for idx in indices:
                detection_sigma = self.parent.parent.session\
                    .metadata["spectral_models"][idx]\
                    .metadata.get("detection_sigma", None)
                if detection_sigma is not None:
                    break
            else:
                detection_sigma = 0.5

            detection_sigma, is_ok = QtGui.QInputDialog.getDouble(
                None, "Set detection sigma", u"Detection sigma:", 
                value=detection_sigma, minValue=0.1, maxValue=1000)
            if not is_ok: return None

            for idx, proxy_index in zip(indices, proxy_indices):
                spectral_model \
                    = self.parent.parent.session.metadata["spectral_models"][idx]

                if isinstance(spectral_model, ProfileFittingModel):
                    spectral_model.metadata["detection_sigma"] = detection_sigma

                    if "fitted_result" in spectral_model.metadata:
                        spectral_model.fit()
                        update_spectrum_figure = True
                        self.parent.table_view.rowMoved(proxy_index.row(), 
                            proxy_index.row(), proxy_index.row())



        elif action == set_detection_pixels:

            # Get the first profile model
            for idx in indices:
                detection_pixels = self.parent.parent.session\
                    .metadata["spectral_models"][idx]\
                    .metadata.get("detection_pixels", None)
                if detection_pixels is not None:
                    break
            else:
                detection_pixels = 3

            detection_pixels, is_ok = QtGui.QInputDialog.getInt(
                None, "Set detection pixel", u"Detection pixels:", 
                value=detection_pixels, minValue=1, maxValue=1000)
            if not is_ok: return None

            for idx, proxy_index in zip(indices, proxy_indices):
                spectral_model \
                    = self.parent.parent.session.metadata["spectral_models"][idx]

                if isinstance(spectral_model, ProfileFittingModel):
                    spectral_model.metadata["detection_pixels"] = detection_pixels

                    if "fitted_result" in spectral_model.metadata:
                        spectral_model.fit()
                        update_spectrum_figure = True
                        self.parent.table_view.rowMoved(proxy_index.row(), 
                            proxy_index.row(), proxy_index.row())


        elif action == set_rv_tolerance:

            # Get the first profile model
            for idx in indices:
                velocity_tolerance = self.parent.parent.session\
                    .metadata["spectral_models"][idx]\
                    .metadata.get("velocity_tolerance", None)
                if velocity_tolerance is not None:
                    break
            else:
                velocity_tolerance = 5

            velocity_tolerance, is_ok = QtGui.QInputDialog.getInt(
                None, "Set velocity tolerance", u"Velocity tolerance:", 
                value=velocity_tolerance, minValue=0.01, maxValue=100)
            if not is_ok: return None

            for idx, proxy_index in zip(indices, proxy_indices):
                spectral_model \
                    = self.parent.parent.session.metadata["spectral_models"][idx]

                if isinstance(spectral_model, ProfileFittingModel):
                    spectral_model.metadata["velocity_tolerance"] = velocity_tolerance

                    if "fitted_result" in spectral_model.metadata:
                        spectral_model.fit()
                        update_spectrum_figure = True
                        self.parent.table_view.rowMoved(proxy_index.row(), 
                            proxy_index.row(), proxy_index.row())



        elif action == set_wl_tolerance:

            # Get the first profile model
            for idx in indices:
                wavelength_tolerance = self.parent.parent.session\
                    .metadata["spectral_models"][idx]\
                    .metadata.get("wavelength_tolerance", None)
                if wavelength_tolerance is not None:
                    break
            else:
                wavelength_tolerance = 0.1

            wavelength_tolerance, is_ok = QtGui.QInputDialog.getInt(
                None, "Set wavelength tolerance", u"Wavelength tolerance:", 
                value=wavelength_tolerance, minValue=0.01, maxValue=1)
            if not is_ok: return None

            for idx, proxy_index in zip(indices, proxy_indices):
                spectral_model \
                    = self.parent.parent.session.metadata["spectral_models"][idx]

                if isinstance(spectral_model, ProfileFittingModel):
                    spectral_model.metadata["wavelength_tolerance"] = wavelength_tolerance

                    if "fitted_result" in spectral_model.metadata:
                        spectral_model.fit()
                        update_spectrum_figure = True
                        self.parent.table_view.rowMoved(proxy_index.row(), 
                            proxy_index.row(), proxy_index.row())


        if update_spectrum_figure:
            self.parent.update_spectrum_figure(redraw=True)
            

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

        if role==QtCore.Qt.FontRole:
            return _QFONT

        column = index.column()
        # TODO this is broken
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

        elif column == 2:
            value = spectral_model._repr_element

        elif column == 3:
            try:
                result = spectral_model.metadata["fitted_result"][2]
                equivalent_width = result["equivalent_width"][0]

            except:
                equivalent_width = np.nan

            value = "{0:.1f}".format(1000 * equivalent_width) \
                if np.isfinite(equivalent_width) else ""

        elif column == 4:

            try:
                result = spectral_model.metadata["fitted_result"][2]
                u_equivalent_width \
                    = 1000 * np.max(np.abs(result["equivalent_width"][1:]))
                value = "{0:.1f}".format(u_equivalent_width)

            except:
                value = ""


        elif column == 5:
            try:
                abundances \
                    = spectral_model.metadata["fitted_result"][2]["abundances"]

            except (IndexError, KeyError):
                value = ""

            else:
                # How many elements were measured?
                value = "; ".join(["{0:.2f}".format(abundance) \
                    for abundance in abundances])


        elif column == 6:

            try:
                result = spectral_model.metadata["fitted_result"][2]
                u_abundance = result["abundance_uncertainties"][0]
                value = "{0:.2f}".format(u_abundance)

            except:
                value = ""

        elif column == 7 \
        and role in (QtCore.Qt.DisplayRole, QtCore.Qt.CheckStateRole):
            value = spectral_model.is_upper_limit
            if role == QtCore.Qt.CheckStateRole:
                return QtCore.Qt.Checked if value else QtCore.Qt.Unchecked
            else:
                return None
            
        return value if role == QtCore.Qt.DisplayRole else None
    

    def setData(self, index, value, role=QtCore.Qt.DisplayRole):
        
        # We only allow the checkbox to be ticked or unticked here. The other
        # columns cannot be edited. 
        # This is handled by the SpectralModelsTableModelBase class.

        # If the checkbox has just been ticked by the user but the selected
        # spectral model does not have a result, we should not allow the
        # spectral model to be marked as acceptable.
        if index.column() == 0 and value and \
        "fitted_result" not in self.spectral_models[index.row()].metadata:
            return False
        
        value = super(SpectralModelsTableModel, self).setData(index, value, role)

        # If we have a cache of the state transitions, update the entries.
        if hasattr(self.parent, "_state_transitions"):
            cols = ("equivalent_width", "reduced_equivalent_width", "abundance")
            for col in cols:
                self.parent._state_transitions[col][index.row()] = np.nan
            self.parent.update_scatter_plots(redraw=False)
            self.parent.update_selected_points(redraw=False)

        # TODO: Any cheaper way to update this?
        #       layoutAboutToBeChanged() and layoutChanged() didn't work
        #       neither did rowCountChanged or rowMoved()
        self.parent.proxy_spectral_models.reset()

        # Update figures.
        self.parent.update_scatter_plots(redraw=False)
        self.parent.update_selected_points(redraw=False)
        self.parent.update_trend_lines(redraw=True)
        
        return value

