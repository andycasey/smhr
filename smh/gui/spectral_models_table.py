#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["SpectralModelsTableView", "SpectralModelsFilterProxyModel", "SpectralModelsTableModelBase"]

import logging
import numpy as np
import sys
from PySide import QtCore, QtGui
from time import time

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
        QtGui.QFont.insertSubstitution(*substitute)


DOUBLE_CLICK_INTERVAL = 0.1 # MAGIC HACK
PICKER_TOLERANCE = 10 # MAGIC HACK



class SpectralModelsTableView(QtGui.QTableView):

    def __init__(self, parent, *args):
        super(SpectralModelsTableView, self).__init__(parent, *args)
        self.parent = parent


    def contextMenuEvent(self, event):
        """
        Provide a context (right-click) menu for the line list table.

        :param event:
            The mouse event that triggered the menu.
        """

        N = len(self.selectionModel().selectedRows())
        
        menu = QtGui.QMenu(self)
        fit_models = menu.addAction("Fit model{}..".format(["", "s"][N != 1]))
        menu.addSeparator()
        option = menu.addAction("Option A")
        option.setEnabled(False)
        option = menu.addAction("Option B")
        option.setEnabled(False)
        option = menu.addAction("Option C")
        option.setEnabled(False)
        option = menu.addAction("Option D")
        option.setEnabled(False)
        menu.addSeparator()
        option = menu.addAction("Option E")
        option.setEnabled(False)

        if N == 0:
            fit_models.setEnabled(False)

        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == fit_models:
            self.fit_selected_models()

        return None


    def fit_selected_models(self):
        """ Fit the selected spectral models. """

        data_model = self.model().sourceModel()
        for i, proxy_index in enumerate(self.selectionModel().selectedIndexes()):
            # Fit the models.

            index = self.model().mapToSource(proxy_index).row()
            self.parent.parent.session.metadata["spectral_models"][index].fit()

            # Update the view if this is the first one.
            if i == 0:
                # TODO parent has update_spectrum_figure()
                self.parent.update_spectrum_figure()

            # Update the data model.
            data_model.dataChanged.emit(
                data_model.createIndex(proxy_index.row(), 0),
                data_model.createIndex(proxy_index.row(),
                    data_model.columnCount(QtCore.QModelIndex())))

            # It ought to be enough just to emit the dataChanged signal, but
            # there is a bug when using proxy models where the data table is
            # updated but the view is not, so we do this hack to make it
            # work:
            # TODO parent has table_view
            self.parent.table_view.rowMoved(
                proxy_index.row(), proxy_index.row(), proxy_index.row())

        return None



class SpectralModelsFilterProxyModel(QtGui.QSortFilterProxyModel):

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


    def reset(self, *args):
        super(SpectralModelsFilterProxyModel, self).reset(*args)
        self.reindex()
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

        return self.createIndex(self.lookup_indices[proxy_index.row()],
            proxy_index.column())



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
        if index.column() != 0 or value:
            return False

        model = self.spectral_models[index.row()]
        print("test test",index.row(),value,model.metadata["is_acceptable"])
        model.metadata["is_acceptable"] = value

        try:
            del model.metadata["fitted_result"]
        except KeyError:
            None

        # Emit data change for this row.
        self.dataChanged.emit(self.createIndex(index.row(), 0),
            self.createIndex(index.row(), 
                self.columnCount(QtCore.QModelIndex())))

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


