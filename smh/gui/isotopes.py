#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The radial velocity tab view for the Spectroscopy Made Hard GUI. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import os
import sys
import yaml
from PySide import QtCore, QtGui

from matplotlib import (gridspec, pyplot as plt)

from smh import (Session, isoutils)
from smh.linelists import LineList

__all__ = ["IsotopeWidget"]

class IsotopeWidget(QtGui.QWidget):
    """
    Widget to edit isotope ratios.

    For simplicity, implemented by copying all isotopes and modifying internally.
    Then it just returns a new (validated) set of isotope ratios.

    Importantly, you MUST close the widget window before it updates the session.
    """

    def __init__(self, orig_isotopes, needed_isotopes=None, linelist=None, parent=None):
        """
        Get the current isotopes and the needed isotopes.
        Setup an editable table.
        """

        super(IsotopeWidget, self).__init__(parent)
        self.parent = parent
        self.setGeometry(300, 200, 570, 450)

        # TODO get orig_isotopes from 
        # parent.session.metadata['isotopes']

        # Ensure that orig_isotopes is actually isotopes
        if not isinstance(orig_isotopes, dict):
            raise Error(str(orig_isotopes))
        for key in orig_isotopes:
            if not isinstance(orig_isotopes[key], dict):
                raise Error(str(orig_isotopes))
        
        # Figure out where to get needed isotopes
        if isinstance(needed_isotopes, dict):
            if isinstance(linelist, LineList):
                raise Error("Cannot specify both needed_isotopes AND linelist")
        else:
            assert needed_isotopes == None, needed_isotopes
            if isinstance(linelist, LineList):
                needed_isotopes = isoutils.identify_needed_isotopes(linelist)
            else:
                assert linelist == None, linelist
                # No needed isotopes!
                needed_isotopes = {}

        # This is the data we will work with
        current_isotopes = orig_isotopes.copy()
        self.orig_isotopes = orig_isotopes
        self.needed_isotopes = needed_isotopes
        self.current_isotopes = current_isotopes

        # Establish the GUI for this tab.
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)

        layout = QtGui.QHBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)

        # Create editable table
        tab = isoutils.convert_isodict_to_array(current_isotopes)
        self.tab = tab

        table_model = IsotopeModel(self, tab)
        table_view  = QtGui.QTableView(self)
        table_view.setModel(table_model)
        table_view.resizeRowsToContents()
        table_view.resizeColumnsToContents()

        # Create and link buttons
        button_widget = QtGui.QWidget(self)

        button_add_row = QtGui.QPushButton("Add Row To Bottom")
        button_add_row.clicked.connect(self.add_row)
        button_ins_row = QtGui.QPushButton("Insert Row")
        button_ins_row.clicked.connect(self.insert_row)
        button_del_row = QtGui.QPushButton("Delete Row")
        button_del_row.clicked.connect(self.delete_row)
        button_rproc  = QtGui.QPushButton("Use r-process")
        button_sproc  = QtGui.QPushButton("Use s-process")
        button_sneden = QtGui.QPushButton("Use Sneden solar")
        button_asplund= QtGui.QPushButton("Use Asplund solar")
        button_finish = QtGui.QPushButton("Finish setting isotopes")
        button_cancel = QtGui.QPushButton("Cancel without saving")

        all_buttons = [button_add_row, button_ins_row, button_del_row, 
                       button_rproc, button_sproc, button_sneden, button_asplund,
                       button_finish, button_cancel]
        button_layout = QtGui.QVBoxLayout(button_widget)
        for button in all_buttons:
            button_layout.addWidget(button)

        # Final widget
        layout.addWidget(table_view)
        layout.addWidget(button_widget)

        self.table_view = table_view
        self.table_model = table_model
        self.button_widget = button_widget
        self.all_buttons = all_buttons

        return None

    def add_row(self):
        # Insert a blank row at the bottom
        self.table_model.insertRow(self.table_model.rowCount(self))
    def insert_row(self):
        index = self.table_view.selectionModel().currentIndex()
        row = index.row()
        self.table_model.insertRow(row)
    def delete_row(self):
        index = self.table_view.selectionModel().currentIndex()
        row = index.row()
        self.table_model.removeRow(row)

    def finish(self):
        """
        Return a new isotope ratio and close the window?
        Also allow canceling?
        """
        new_isotopes = isoutils.convert_array_to_isodict(self.tab)
        # Modify the session's isotope ratios
        raise NotImplementedError

class IsotopeModel(QtCore.QAbstractTableModel):
    def __init__(self, parent, tab, *args):
        QtCore.QAbstractTableModel.__init__(self, parent, *args)
        # Expand dictionary to a table
        self.tab = tab
    def rowCount(self, parent):
        return len(self.tab)
    def columnCount(self, parent):
        return self.tab.shape[1]
    def data(self, index, role):
        if not index.isValid():
            return None
        elif role != QtCore.Qt.DisplayRole:
            return None
        return self.tab[index.row(),index.column()]
    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            if col == 0: return 'Elem'
            if col == 1: return 'A'
            if col == 2: return 'Frac'
        return None

    def setData(self,index,value,role):
        self.tab[index.row(),index.column()] = value
        #QtCore.QAbstractItemModel.dataChanged(self)
        return True
    def flags(self,index):
        return QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsEditable
    
    def insertRow(self,row,parent=QtCore.QModelIndex()):
        QtCore.QAbstractItemModel.beginInsertRows(self,parent,row,row)
        self.tab = np.insert(self.tab,row,['??','xxx','1.0'],axis=0)
        QtCore.QAbstractItemModel.endInsertRows(self)
    def removeRow(self,row,parent=QtCore.QModelIndex()):
        QtCore.QAbstractItemModel.beginRemoveRows(self,parent,row,row)
        self.tab = np.delete(self.tab,row,axis=0)
        QtCore.QAbstractItemModel.endRemoveRows(self)
