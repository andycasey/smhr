#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The radial velocity tab view for the Spectroscopy Made Hard GUI. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import os
import sys
from PySide2 import (QtCore, QtWidgets as QtGui)

from matplotlib import (gridspec, pyplot as plt)

from smh import (Session, isoutils)
from smh.linelists import LineList

import logging
logger = logging.getLogger(__name__)

__all__ = ["IsotopeDialog","IsotopeError","IsotopeLog"]

class IsotopeDialog(QtGui.QDialog):
    """
    Widget to edit isotope ratios.

    For simplicity, implemented by copying all isotopes and modifying internally.
    Then it just returns a new (validated) set of isotope ratios.

    Importantly, you MUST push the "Finish" button to update the session
    """

    def __init__(self, session=None, orig_isotopes=None, needed_isotopes=None, linelist=None):
        """
        Get the current isotopes and the needed isotopes.
        Setup an editable table.
        """

        super(IsotopeDialog, self).__init__()
        self.setGeometry(300, 200, 570, 450)

        if isinstance(session, Session):
            if isinstance(orig_isotopes, dict):
                logger.warn("IsotopeDialog: Session already has isotopes, clobbering...")
            orig_isotopes = session.metadata['isotopes']
            if isinstance(linelist, LineList):
                logger.warn("IsotopeDialog: Session already has isotopes, clobbering...")
            if isinstance(needed_isotopes, dict):
                logger.warn("IsotopeDialog: Using session, clobbering needed_isotopes...")
            needed_isotopes = None
            if 'line_list' in session.metadata.keys():
                linelist = session.metadata['line_list']
        self.session = session

        # Ensure that orig_isotopes is actually isotopes
        # Doesn't have to be validated though
        if not isinstance(orig_isotopes, dict):
            raise Error("IsotopeDialog: orig_isotopes is not a dict \n"+str(orig_isotopes))
        for key in orig_isotopes:
            if not isinstance(orig_isotopes[key], dict):
                raise Error("IsotopeDialog: orig_isotopes values are not dicts \n"+str(orig_isotopes))
        
        # Figure out where to get needed isotopes
        if isinstance(needed_isotopes, dict):
            if isinstance(linelist, LineList):
                raise Error("IsotopeDialog: Cannot specify both needed_isotopes and linelist")
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

        self.default_isotopes = [isoutils.load_isotope_data(x) for x in ['rproc','sproc','sneden','asplund']]

        # Establish the GUI for this tab.
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)

        layout = QtGui.QHBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)

        # Create editable table
        tab = isoutils.convert_isodict_to_array(current_isotopes)

        table_model = IsotopeModel(self, tab)
        table_view  = IsotopeTableView(self)
        table_view.setModel(table_model)
        table_view.resizeRowsToContents()
        table_view.resizeColumnsToContents()

        # Create and link buttons
        button_widget = QtGui.QWidget(self)
        text_box = IsotopeLog(self)
        text_box.setReadOnly(True)
        
        button_needed = QtGui.QPushButton("Print Needed Isotopes")
        button_needed.clicked.connect(self.print_needed_isotopes)
        button_add_row = QtGui.QPushButton("Add Row To Bottom")
        button_add_row.clicked.connect(self.add_row)
        button_ins_row = QtGui.QPushButton("Insert Row")
        button_ins_row.clicked.connect(self.insert_row)
        button_del_row = QtGui.QPushButton("Delete Row")
        button_del_row.clicked.connect(self.delete_row)
        button_rproc  = QtGui.QPushButton("Use r-process (Sneden08)")
        button_rproc.clicked.connect(self.use_rproc)
        button_sproc  = QtGui.QPushButton("Use s-process (Sneden08)")
        button_sproc.clicked.connect(self.use_sproc)
        button_sneden = QtGui.QPushButton("Use Sneden08 solar")
        button_sneden.clicked.connect(self.use_sneden)
        button_asplund= QtGui.QPushButton("Use Asplund09 solar")
        button_asplund.clicked.connect(self.use_asplund)
        button_finish = QtGui.QPushButton("Save isotopes to session and exit")
        button_finish.clicked.connect(self.finish)
        if not isinstance(session, Session):
            button_finish.setEnabled(False)
        button_cancel = QtGui.QPushButton("Restore original isotopes and exit (TODO)")
        button_cancel.clicked.connect(self.reject) #built-in finish button

        all_buttons = [button_needed, button_add_row, button_ins_row, button_del_row, 
                       button_rproc, button_sproc, button_sneden, button_asplund,
                       button_finish, button_cancel]
        button_layout = QtGui.QVBoxLayout(button_widget)
        button_layout.addWidget(text_box)
        for button in all_buttons:
            button_layout.addWidget(button)

        # Final widget
        layout.addWidget(table_view)
        layout.addWidget(button_widget)

        self.table_view = table_view
        self.table_model = table_model
        self.text_box = text_box
        self.button_widget = button_widget
        self.all_buttons = all_buttons

        self.print_needed_isotopes()

        return None

    def print_needed_isotopes(self):
        # TODO figure out new needed isotopes
        self.log("Needed isotopes:\n"+isoutils.pretty_print_isotopes(self.needed_isotopes))

    def _use_data(self,i):
        row = self.table_view.selectionModel().currentIndex().row()
        elem = self.table_model.tab[row,0]
        try:
            self.table_model.set_element_with_isotopes(elem,self.default_isotopes[i])
        except ValueError as e:
            self.log(str(e))
        
    #def set_colors(self):
    #    self.table_model.setData(self.table_model.index(0,0), QtCore.Qt.red, QtCore.Qt.BackgroundRole)

    def use_rproc(self):
        self._use_data(0)
        self.log('(using rproc)')
    def use_sproc(self):
        self._use_data(1)
        self.log('(using sproc)')
    def use_sneden(self):
        self._use_data(2)
        self.log('(using sneden)')
    def use_asplund(self):
        self._use_data(3)
        self.log('(using asplund)')

    def add_row(self):
        # Insert a blank row at the bottom
        self.table_model.insertRow(self.table_model.rowCount(self))
        self.log("Added row at bottom")
    def insert_row(self):
        index = self.table_view.selectionModel().currentIndex()
        row = index.row()
        self.table_model.insertRow(row)
        self.log("Insert row at row {}".format(row))
    def delete_row(self):
        index = self.table_view.selectionModel().currentIndex()
        row = index.row()
        elem,A,frac = self.table_model.tab[row,:]
        self.log("Deleting row ({},{},{})".format(elem,A,frac))
        self.table_model.removeRow(row)

    def validate(self):
        ## TODO add this
        self.log("Validating isotopes...")
        try:
            new_isotopes = isoutils.convert_array_to_isodict(self.table_model.tab)
        except Exception as e:
            self.log("ERROR: Can't create isotope dict (invalid values?)")
            self.log(str(e))
            return None
        try:
            isoutils.validate_isotopes(new_isotopes)
        except isoutils.IsotopeError as e:
            self.log("ERROR: Invalid isotopes!")
            self.log("Bad isotopes:")
            self.log(isoutils.pretty_print_isotopes(e.bad_isotopes))
        else:
            self.log("Isotopes validated!")
        return new_isotopes
    def finish(self):
        """
        Return a new isotope ratio and close the window
        """
        try:
            new_isotopes = isoutils.convert_array_to_isodict(self.table_model.tab)
        except Exception as e:
            self.log(str(e))
            return None

        try:
            isoutils.validate_isotopes(new_isotopes)
        except isoutils.IsotopeError as e:
            self.log("ERROR: Invalid isotopes! Not saving.")
            self.log("Bad isotopes:")
            self.log(isoutils.pretty_print_isotopes(e.bad_isotopes))
        else:
            self.log("Isotopes validated!")
            self.log("Saved isotopes to session")
            self.session.metadata['isotopes'] = new_isotopes
            self.accept()
    def log(self,msg):
        self.text_box.appendMessage(msg)

class IsotopeTableView(QtGui.QTableView):
    def __init__(self, parent, *args):
        super(IsotopeTableView, self).__init__(parent, *args)
        self.parent = parent
    def contextMenuEvent(self, event):
        """
        Right-click menu
        """
        menu = QtGui.QMenu(self)
        insert_action = menu.addAction("Insert Line")
        delete_action = menu.addAction("Delete Line")
        
        insert_action.setEnabled(False)
        delete_action.setEnabled(False)
#        action = menu.exec_(self.mapToGlobal(event.pos()))
#        if action == delete_action:
#            self.parent.delete_row()
#        elif action == insert_action:
#            self.parent.insert_row()

class IsotopeModel(QtCore.QAbstractTableModel):
    def __init__(self, parent, tab, *args):
        QtCore.QAbstractTableModel.__init__(self, parent, *args)
        self.parent = parent
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

    def set_element_with_isotopes(self,elem,isotopes):
        if elem in isotopes:
            self.parent.text_box.appendMessage("Updating {}".format(elem))
            self.beginResetModel()
            isodict = isoutils.convert_array_to_isodict(self.tab)
            isodict[elem] = isotopes[elem]
            self.tab = isoutils.convert_isodict_to_array(isodict)
            self.endResetModel()
        else:
            self.parent.text_box.appendMessage("{} not in isotopes".format(elem))

    def setData(self,index,value,role):
        oldval = self.tab[index.row(),index.column()]
        self.tab[index.row(),index.column()] = value
        self.parent.text_box.appendMessage("Changed {} to {}".format(oldval,value))
        return True
    def flags(self,index):
        return QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsEditable
    
    def insertRow(self,row,parent=QtCore.QModelIndex()):
        self.beginInsertRows(parent,row,row)
        self.tab = np.insert(self.tab,row,['??','0','1.0'],axis=0)
        self.endInsertRows()
    def removeRow(self,row,parent=QtCore.QModelIndex()):
        self.beginRemoveRows(parent,row,row)
        self.tab = np.delete(self.tab,row,axis=0)
        self.endRemoveRows()

class IsotopeLog(QtGui.QPlainTextEdit):
    def log(self,msg):
        self.appendMessage(msg)
    def appendMessage(self,msg):
        self.appendPlainText(msg)
        self.verticalScrollBar().setValue(self.verticalScrollBar().maximum())
        
if __name__=="__main__":
    # For development
    datadir = os.path.dirname(os.path.abspath(__file__))+'/../tests/test_data'
    session = Session([datadir+"/spectra/hd122563.fits"])

    orig_isotopes = isoutils.load_isotope_data('rproc')
    ll = LineList.read(os.path.dirname(os.path.abspath(__file__))+'/../tests/test_data/linelists/lin4554new')
    session.metadata['isotopes'] = orig_isotopes
    session.metadata['line_list'] = ll

    app = QtGui.QApplication(sys.argv)
    window = IsotopeDialog(session) #orig_isotopes=orig_isotopes, linelist=ll)
    window.exec_()
    #print(session.metadata['isotopes'])
    #print('U in isotopes?', 'U' in session.metadata['isotopes'])
