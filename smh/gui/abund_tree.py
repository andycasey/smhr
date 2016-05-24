from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import os
import sys
from PySide import QtCore, QtGui
from six import string_types

from smh import (Session, specutils)
from smh.linelists import LineList

from astropy import table

import logging
logger = logging.getLogger(__name__)

_treecols = ['is_selected','wavelength','expot','loggf','element','A(X)','e(X)','equivalent_width','A(X)','A(X)']
_treecolmap = dict(zip(range(len(_treecols)),_treecols))

def summarize_abundances_species(ttab):
    element = ttab['element'][0]
    ttab = ttab[ttab['is_selected']]
    N = len(ttab)
    if N==0: return [element, N, np.nan, np.nan, np.nan, np.nan]
    weights = 1./ttab['e(X)']**2
    total_weights = np.sum(weights)
    abund = np.mean(ttab['A(X)']*weights)/total_weights
    stdev = np.std(ttab['A(X)']) #TODO weight
    XH = np.nan
    XFe = np.nan
    return [element,N,abund,stdev,XH,XFe]
def summarize_abundances(tab):
    tab = tab.group_by('species')
    summary = []
    for i,species in enumerate(tab.groups.keys):
        ttab = tab.groups[i]
        summary.append(summarize_abundances_species(ttab))
    return summary

class AbundTreeView(QtGui.QTreeView):
    pass
class AbundTreeItem(object):
    def __init__(self, parent, row):
        self.parent = parent
        self.row = row
        self.subnodes = self._getChildren()
    def _getChildren(self):
        raise NotImplementedError()
    def child(self, row):
        return self.subnodes[row]
    def childCount(self):
        return len(self.subnodes)
class AbundTreeElementSummaryItem(AbundTreeItem):
    def __init__(self, parent, tab, index, model):
        self.index = index
        self.model = model
        super(AbundTreeElementSummaryItem, self).__init__(parent, index)
        
        self.compute_summary()
    def _getChildren(self):
        N = len(self.model.tab.groups[self.index])
        return [AbundTreeMeasurementItem(row,self) for row in range(N)]
    def columnCount(self):
        return 6
    def data(self, column):
        if column >= self.columnCount(): return None
        return str(self.summary[column])
    def compute_summary(self):
        ttab = self.model.tab.groups[self.index]
        self.summary = summarize_abundances_species(ttab)
        return None

class AbundTreeMeasurementItem(AbundTreeItem):
    def __init__(self,row,parent):
        self.row = row
        assert isinstance(parent, AbundTreeElementSummaryItem)
        self.parent = parent
        super(AbundTreeMeasurementItem, self).__init__(parent, row)
    def _getChildren(self):
        return []
    def columnCount(self):
        return 10
    def data(self, column):
        col = _treecolmap[column]
        data = self.parent.model.tab.groups[self.parent.index][self.row][col]
        if col=='is_selected': return data
        return str(data)
    
class AbundTreeModel(QtCore.QAbstractItemModel):
    def __init__(self, session=None, parent=None):
        super(AbundTreeModel, self).__init__(parent)
        self.session = session
        self.tab = self.obtain_measurements_from_session(session)
        # TODO set up a map from spectral models to items in the tree
        # That way you can selectively update the internal table here
        self.summaries = self._getSummaries()

    def session_updated(self,spectral_model=None):
        # TODO call this whenever the session updates measurements in any way
        if spectral_model is None:
            self.obtain_measurements_from_session()
        else:
            # TODO get the item from the measurement
            # TODO update just those items with the new measurements
            raise NotImplementedError
    def obtain_measurements_from_session(self):
        # TODO actually do this from the session
        session = self.session
        ll = LineList.read("/Users/alexji/smhr/smh/tests/test_data/linelists/complete.list")
        col1 = table.Column(np.ones(len(ll)),name='A(X)')
        col2 = table.Column(np.ones(len(ll)),name='e(X)')
        col3 = table.Column(np.ones(len(ll), dtype=bool),name='is_selected')
        ll.add_columns([col1,col2,col3])
        ll['equivalent_width'] = 2.0
        ll = table.Table(ll)
        cols = ['wavelength','expot','loggf','element','species','A(X)','e(X)','is_selected','equivalent_width']
        ll = ll[cols]
        tab = ll.group_by('species')
        return tab

    def _getSummaries(self):
        summaries = []
        # Initialize summary objects for each element
        for index,ttab in enumerate(self.tab.groups):
            # Initialize the abundance summary 
            elem_summary = AbundTreeElementSummaryItem(None, self.tab, index, self)
            summaries.append(elem_summary)
        return summaries

    def index(self, row, column, parent):
        if not parent.isValid(): # Root node
            return self.createIndex(row, column, self.summaries[row])
        parentNode = parent.internalPointer()
        return self.createIndex(row, column, parentNode.subnodes[row])
    def parent(self, index):
        if not index.isValid():
            return QtCore.QModelIndex()
        node = index.internalPointer()
        if node.parent is None:
            return QtCore.QModelIndex()
        else:
            return self.createIndex(node.parent.row, 0, node.parent)
    def reset(self):
        self.summaries = self._getSummaries()
        QtCore.QAbstractItemModel.reset(self)
    def rowCount(self, parent):
        if not parent.isValid():
            return len(self.summaries)
        node = parent.internalPointer()
        return len(node.subnodes)
    def columnCount(self, parent):
        if isinstance(parent, AbundTreeElementSummaryItem):
            return 6
        return 10

    def data(self, index, role):
        if not index.isValid(): return None
        if role == QtCore.Qt.CheckStateRole and index.column()==0:
            item = index.internalPointer()
            if isinstance(item, AbundTreeMeasurementItem):
                checked = item.data(0)# == 'True' #grrr
                if checked: return QtCore.Qt.Checked
                else: return QtCore.Qt.Unchecked
            else: return None
        if role != QtCore.Qt.DisplayRole: return None
        item = index.internalPointer()
        return item.data(index.column())

    def setData(self, index, value, role):
        print("Setting data")
        item = index.internalPointer()
        if isinstance(item, AbundTreeMeasurementItem) and index.column()==0:
            # Check or uncheck a box
            if role == QtCore.Qt.EditRole: return False
            if role == QtCore.Qt.CheckStateRole:
                self.tab.groups[item.parent.index][item.row]['is_selected'] = value
                # Also edit the parent summary here!
                item.parent.compute_summary()
                topLeft = self.createIndex(0, 0, item.parent)
                botRight = self.createIndex(0, item.parent.columnCount(), item.parent)
                self.dataChanged.emit(index,index)
                self.dataChanged.emit(topLeft,botRight)
                return True
        return False
            
    def flags(self, index):
        if not index.isValid(): return None
        item = index.internalPointer()
        if isinstance(item, AbundTreeElementSummaryItem):
            pass
        if isinstance(item, AbundTreeMeasurementItem):
            if index.column() == 0:
                return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsUserCheckable
        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable

    def headerData(self, section, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            cols = ['','wl','EP','loggf','elem','A(X)','e(X)','EW','A(X)','A(X)']
            if section < len(cols):
                return cols[section]
        return None
        
if __name__=="__main__":
    app = QtGui.QApplication(sys.argv)
    abundtree = AbundTreeView()
    model = AbundTreeModel()
    abundtree.setModel(model)
    
    abundtree.show()
    sys.exit(app.exec_())
