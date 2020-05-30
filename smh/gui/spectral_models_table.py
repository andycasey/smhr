#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["SpectralModelsTableViewBase", "SpectralModelsFilterProxyModel", "SpectralModelsTableModelBase"]

import logging
import numpy as np
import sys
from PySide2 import (QtCore, QtGui as QtGui2, QtWidgets as QtGui)
import time

from smh.photospheres import available as available_photospheres
from smh.spectral_models import (ProfileFittingModel, SpectralSynthesisModel)
from smh import utils
from linelist_manager import TransitionsDialog

logger = logging.getLogger(__name__)

if sys.platform == "darwin":
        
    # See http://successfulsoftware.net/2013/10/23/fixing-qt-4-for-mac-os-x-10-9-mavericks/
    substitutes = [
        (".Lucida Grande UI", "Lucida Grande"),
        (".Helvetica Neue DeskInterface", "Helvetica Neue")
    ]
    for substitute in substitutes:
        QtGui2.QFont.insertSubstitution(*substitute)


DOUBLE_CLICK_INTERVAL = 0.1 # MAGIC HACK
PICKER_TOLERANCE = 10 # MAGIC HACK



class SpectralModelsTableViewBase(QtGui.QTableView):

    def __init__(self, parent, *args):
        super(SpectralModelsTableViewBase, self).__init__(parent, *args)
        self.parent = parent

    def update_row(self,row):
        """
        Row is proxy_index.row()
        """
        #data_model = self.model().sourceModel()
        #start = time.time()
        #data_model.dataChanged.emit(
        #    data_model.createIndex(row, 0),
        #    data_model.createIndex(row,data_model.columnCount(None)))
        #print("Time to emit update_row: {:.1f}".format(time.time()-start))
        # It ought to be enough just to emit the dataChanged signal, but
        # there is a bug when using proxy models where the data table is
        # updated but the view is not, so we do this hack to make it
        # work:
        self.rowMoved(row, row, row)
        return None

    def contextMenuEvent(self, event):
        """
        Provide a context (right-click) menu for the line list table.

        :param event:
            The mouse event that triggered the menu.
        """

        proxy_indices = self.selectionModel().selectedRows()
        N = len(proxy_indices)
        
        menu = QtGui.QMenu(self)
        fit_models = menu.addAction("Fit selected model{}..".format(["", "s"][N != 1]))
        menu.addSeparator()
        measure_models = menu.addAction("Measure selected model{}..".format(["", "s"][N != 1]))
        menu.addSeparator()
        mark_as_acceptable = menu.addAction("Mark as acceptable")
        mark_as_unacceptable = menu.addAction("Mark as unacceptable")
        menu.addSeparator()

        # Fitting Options
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
        if action == fit_models:
            self.fit_selected_models(proxy_indices)
        elif action == measure_models:
            self.measure_selected_models(proxy_indices)
        elif action == mark_as_acceptable:
            self.mark_selected_models_as_acceptable(proxy_indices)
        elif action == mark_as_unacceptable:
            self.mark_selected_models_as_unacceptable(proxy_indices)
        elif action == set_fitting_window:
            self.set_fitting_window(proxy_indices)
        elif action == set_no_continuum:
            self.set_no_continuum(proxy_indices)
        elif action in set_continuum_order:
            order = set_continuum_order.index(action)
            self.set_continuum_order(proxy_indices, order)
        elif action in (set_gaussian, set_lorentzian, set_voigt):
            kind = {
                set_gaussian: "gaussian",
                set_lorentzian: "lorentzian",
                set_voigt: "voigt"
            }[action]
            self.set_profile(proxy_indices, kind)
        elif action in (enable_central_weighting, disable_central_weighting):
            toggle = action==enable_central_weighting
            self.set_central_weighting(proxy_indices, toggle)
        elif action == set_detection_sigma:
            self.set_detection_sigma(proxy_indices)
        elif action == set_detection_pixels:
            self.set_detection_pixels(proxy_indices)
        elif action == set_rv_tolerance:
            self.set_rv_tolerance(proxy_indices)
        elif action == set_wl_tolerance:
            self.set_wl_tolerance(proxy_indices)

        return None

    def fit_selected_models(self,indices):
        raise NotImplementedError("Base must be subclassed")
    def measure_selected_models(self,indices):
        raise NotImplementedError("Base must be subclassed")
    def mark_selected_models_as_acceptable(self,indices):
        raise NotImplementedError("Base must be subclassed")
    def mark_selected_models_as_unacceptable(self,indices):
        raise NotImplementedError("Base must be subclassed")
    def set_fitting_option_value(self, proxy_indices, key, value,
                                 valid_for_profile=False,
                                 valid_for_synth=False):
        """
        This should be called to update all fitting options
        and keep the GUI up to date.
        """
        raise NotImplementedError("Base must be subclassed")
    def set_fitting_window(self, proxy_indices):
        window, is_ok = QtGui.QInputDialog.getDouble(
            None, "Set fitting window", u"Fitting window (Å):", 
            value=5, minValue=0.1, maxValue=1000)
        if not is_ok: return None
        self.set_fitting_option_value(proxy_indices, 
                                      "window", window,
                                      valid_for_profile=True,
                                      valid_for_synth=True)
        return None
    def set_no_continuum(self, proxy_indices):
        self.set_fitting_option_value(proxy_indices, 
                                      "continuum_order", -1,
                                      valid_for_profile=True,
                                      valid_for_synth=True)
        return None
    def set_continuum_order(self, proxy_indices, order):
        self.set_fitting_option_value(proxy_indices, 
                                      "continuum_order", order,
                                      valid_for_profile=True,
                                      valid_for_synth=True)
        return None
    def set_profile(self, proxy_indices, kind):
        self.set_fitting_option_value(proxy_indices, 
                                      "profile", kind,
                                      valid_for_profile=True,
                                      valid_for_synth=False)
        return None
    def set_central_weighting(self, proxy_indices, toggle):
        self.set_fitting_option_value(proxy_indices, 
                                      "central_weighting", toggle,
                                      valid_for_profile=True,
                                      valid_for_synth=False)
        return None
    def set_detection_sigma(self, proxy_indices):
        detection_sigma, is_ok = QtGui.QInputDialog.getDouble(
            None, "Set detection sigma", u"Detection sigma:", 
            value=0.5, minValue=0.1, maxValue=1000)
        if not is_ok: return None
        self.set_fitting_option_value(proxy_indices, 
                                      "detection_sigma", detection_sigma,
                                      valid_for_profile=True,
                                      valid_for_synth=False)
    def set_detection_pixels(self, proxy_indices):
        detection_pixels, is_ok = QtGui.QInputDialog.getInt(
            None, "Set detection pixel", u"Detection pixels:", 
            value=3, minValue=1, maxValue=1000)
        if not is_ok: return None
        self.set_fitting_option_value(proxy_indices,
                                      "detection_pixels", detection_pixels,
                                      valid_for_profile=True,
                                      valid_for_synth=False)
    def set_rv_tolerance(self, proxy_indices):
        velocity_tolerance, is_ok = QtGui.QInputDialog.getDouble(
            None, "Set velocity tolerance", u"Velocity tolerance:", 
            value=5, minValue=0.01, maxValue=100)
        if not is_ok: return None
        # TODO cannot turn it back to None right now
        self.set_fitting_option_value(proxy_indices,
                                      "velocity_tolerance", velocity_tolerance,
                                      valid_for_profile=True,
                                      valid_for_synth=True)
    def set_wl_tolerance(self, proxy_indices):
        wavelength_tolerance, is_ok = QtGui.QInputDialog.getDouble(
            None, "Set wavelength tolerance", u"Wavelength tolerance:", 
            value=0.1, minValue=0.01, maxValue=1)
        if not is_ok: return None
        # TODO cannot turn it back to None right now
        self.set_fitting_option_value(proxy_indices,
                                      "wavelength_tolerance", wavelength_tolerance,
                                      valid_for_profile=True,
                                      valid_for_synth=False)

class SpectralModelsFilterProxyModel(QtCore.QSortFilterProxyModel):

    def __init__(self, parent=None):
        super(SpectralModelsFilterProxyModel, self).__init__(parent)
        self.filter_functions = {}
        self.filter_indices = []
        return None


    def add_filter_function(self, name, filter_function):
        """
        Add a filtering function to the proxy model.

        :param name:
            The name of the filtering function.

        :param filter_function:
            A function that accepts a single argument (the spectral model) and
            returns True or False whether to display this row in the table.
        """

        self.filter_functions[name] = filter_function
        self.invalidateFilter()
        self.reindex()
        return None


    def delete_filter_function(self, name):
        """
        Delete a filtering function from the proxy model.

        :param name:
            The name of the filtering function:
        """


        try:
            del self.filter_functions[name]
            self.invalidateFilter()
            self.reindex()

        except KeyError:
            raise

        else:
            return None

    def delete_all_filter_functions(self):
        self.filter_functions = {}
        self.filter_indices = []
        self.invalidateFilter()
        self.reindex()
        return None

    def reset(self, *args):
        #super(SpectralModelsFilterProxyModel, self).reset(*args)
        super().beginResetModel()
        self.reindex()
        super().endResetModel()
        return None


    def reindex(self):

        try: 
            self.sourceModel().spectral_models

        except AttributeError:
            return None

        lookup_indices = []
        for i, model in enumerate(self.sourceModel().spectral_models):
            for name, filter_function in self.filter_functions.items():
                if not filter_function(model):
                    break
            else:
                # No problems with any filter functions.
                lookup_indices.append(i)

        self.lookup_indices = np.array(lookup_indices)
        return None


    def filterAcceptsRow(self, row, parent):
        """
        Return whether all of the filters for this proxy model agree that this
        row should be shown.

        :param row:
            The row to check.

        :param parent:
            The parent widget.
        """

        # Check if we need to update the filter indices for mapping.
        model = self.sourceModel().spectral_models[row]
        for filter_name, filter_function in self.filter_functions.items():
            if not filter_function(model): break
        else:
            # No problems.
            return True

        # We broke out of the for loop.
        return False


    def mapFromSource(self, data_index):
        """
        Map a data table index to a proxy table index.

        :param data_index:
            The index of the item in the data table.
        """
        if not data_index.isValid():
            return data_index

        # TODO is this necessary every time?
        #self.reindex()

        return self.createIndex(
            np.where(self.lookup_indices == data_index.row())[0],
            data_index.column())


    def mapToSource(self, proxy_index):
        """
        Map a proxy data table index back to the source data table indices.

        :param proxy_index:
            The index of the item in the table.
        """

        if not proxy_index.isValid():
            return proxy_index

        # TODO is this necessary every time?
        #self.reindex()

        try:
            return self.createIndex(self.lookup_indices[proxy_index.row()],
                proxy_index.column())
            
        except AttributeError:
            return proxy_index


class SpectralModelsTableModelBase(QtCore.QAbstractTableModel):

    def __init__(self, parent, header, attrs, *args):
        """
        An abstract table model for spectral models.
        Need to subclass and specify .data() method function!

        :param parent:
            The parent. This *must* have an attribute of `parent.parent.session`.
        """

        super(SpectralModelsTableModelBase, self).__init__(parent, *args)

        # Normally you should never do this, but here I know "better". See:
        #http://stackoverflow.com/questions/867938/qabstractitemmodel-parent-why
        self.parent = parent 
        self.header = header
        self.attrs = attrs
        return None

    @property
    def spectral_models(self):
        try:
            return self.parent.parent.session.metadata.get("spectral_models", [])
        except AttributeError:
            # session is None
            return []

    def rowCount(self, parent):
        """ Return the number of rows in the table. """
        return len(self.spectral_models)

    def columnCount(self, parent):
        """ Return the number of columns in the table. """
        return len(self.header)

    def data(self, index, role):
        raise NotImplementedError("Must subclass this model")

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal \
        and role == QtCore.Qt.DisplayRole:
            return self.header[col]
        return None

    def setData(self, index, value, role=QtCore.Qt.DisplayRole):
        if index.column() != 0:
            return False
        else:
            # value appears to be 0 or 2. Set it to True or False
            row = index.row()
            value = (value != 0)
            model = self.spectral_models[row]
            model.metadata["is_acceptable"] = value
            
            # Emit data change for this row.
            # TODO this is slow.
            #_start = time.time()
            #self.dataChanged.emit(self.createIndex(row, 0),
            #                      self.createIndex(row, 
            #                      self.columnCount(None)))
            #print("Time to emit setData: {:.1f}s".format(time.time()-_start))
            return value
    """
    def sort(self, column, order):
        print("NO SORTING")
        return None

        self.emit(QtCore.SIGNAL("layoutAboutToBeChanged()"))

        def get_equivalent_width(model):
            try:
                return model.metadata["fitted_result"][-1]["equivalent_width"][0]
            except (IndexError, KeyError):
                return np.nan

        def get_abundance(model):
            try:
                return model.metadata["fitted_result"][-1]["abundances"][0]
            except (IndexError, KeyError):
                return np.nan

        sorter = {
            0: lambda model: model.is_acceptable,
            1: lambda model: model._repr_wavelength,
            2: lambda model: model._repr_element,
            3: get_equivalent_width,
            4: get_abundance,
        }

        self._parent.parent.session.metadata["spectral_models"] \
            = sorted(self._parent.parent.session.metadata["spectral_models"], key=sorter[column])
        
        if order == QtCore.Qt.DescendingOrder:
            self._parent.parent.session.metadata["spectral_models"].reverse()

        self.dataChanged.emit(self.createIndex(0, 0),
            self.createIndex(self.rowCount(0), self.columnCount(0)))
        self.emit(QtCore.SIGNAL("layoutChanged()"))
        return None

    """

    def flags(self, index):
        if not index.isValid(): return
        return  QtCore.Qt.ItemIsSelectable|\
                QtCore.Qt.ItemIsEnabled|\
                QtCore.Qt.ItemIsUserCheckable


