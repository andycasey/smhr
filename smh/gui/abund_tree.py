from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import os
import sys
from PySide import QtCore, QtGui
from six import string_types

from smh import (Session, specutils)
from smh.linelists import LineList

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
    def __init__(self, parent, row, tab, index, model):
        self.tab = tab
        self.index = index
        self.model = model
        super(AbundTreeElementSummaryItem, self).__init__(parent, row)
        
        self.summary = []
        self.compute_summary()
    def _getChildren(self):
        N = len(self.tab.groups[self.index])
        return [AbundTreeMeasurementItem(row,self) for row in range(N)]
    def columnCount(self):
        return 6
    def data(self, column):
        return str(self.summary[column])

    # TODO: if my measurements changed, recompute summary
    def compute_summary(self):
        ttab = self.tab.groups[self.index]
        element = ttab[0]['element']
        N = len(ttab)
        self.summary = summarize_abundances_species(ttab)
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
        return str(self.parent.tab.groups[self.parent.index][self.row][_treecolmap[column]])
    
class AbundTreeModel(QtCore.QAbstractItemModel):
    """
    """
    def __init__(self, session=None, parent=None):
        self.session = session # TODO

        ll = LineList.read("/Users/alexji/smhr/smh/tests/test_data/linelists/complete.list")
        from astropy import table
        col1 = table.Column(np.ones(len(ll)),name='A(X)')
        col2 = table.Column(np.ones(len(ll)),name='e(X)')
        col3 = table.Column(np.ones(len(ll), dtype=bool),name='is_selected')
        ll.add_columns([col1,col2,col3])
        ll['equivalent_width'] = 2.0
        ll = table.Table(ll)
        cols = ['wavelength','expot','loggf','element','species','A(X)','e(X)','is_selected','equivalent_width']
        ll = ll[cols]
        tab = ll.group_by('species')
        self.tab = tab
        
        self.summaries = self._getSummaries()
        super(AbundTreeModel, self).__init__(parent)

    def _getSummaries(self):
        summaries = []
        # Initialize summary objects for each element
        for index,ttab in enumerate(self.tab.groups):
            # Initialize the abundance summary 
            parent = None
            elem_summary = AbundTreeElementSummaryItem(parent, index, self.tab, index, self)
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
        self.rootNodes = self._getRootNodes()
        QtCore.QAbstractItemModel.reset(self)
    def rowCount(self, parent):
        if not parent.isValid():
            return len(self.summaries)
        node = parent.internalPointer()
        return len(node.subnodes)
    def columnCount(self, parent):
        return 6
        #if isinstance(parent, AbundTreeElementSummaryItem):
        #    return 6
        #return 10

    def data(self, index, role):
        if not index.isValid(): return None
        if role == QtCore.Qt.CheckStateRole and index.column()==0:
            item = index.internalPointer()
            if isinstance(item, AbundTreeMeasurementItem):
                checked = item.data(0) == 'True'
                if checked: return QtCore.Qt.Checked
                else: return QtCore.Qt.Unchecked
            else: return None
        if role != QtCore.Qt.DisplayRole: return None
            
        item = index.internalPointer()
        return item.data(index.column())
    def setData(self, index, value, role):
        # TODO there is still a problem with the checkbox not changing state!
        item = index.internalPointer()
        if isinstance(item, AbundTreeMeasurementItem) and index.column()==0:
            if role == QtCore.Qt.EditRole: return False
            if role == QtCore.Qt.CheckStateRole:
                self.tab.groups[item.parent.index][item.row]['is_selected'] = value
                self.dataChanged.emit(index,index)
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

#if __name__=="__main__":
def test():
    app = QtGui.QApplication(sys.argv)
    abundtree = QtGui.QTreeWidget()

    ll = LineList.read("/Users/alexji/smhr/smh/tests/test_data/linelists/complete.list")
    from astropy import table
    col1 = table.Column(np.ones(len(ll)),name='A(X)')
    col2 = table.Column(np.ones(len(ll)),name='e(X)')
    col3 = table.Column(np.ones(len(ll),dtype=bool),name='is_selected')
    ll.add_columns([col1,col2,col3])
    ll['equivalent_width'] = 2.0
    ll = table.Table(ll)
    cols = ['wavelength','expot','loggf','element','species','A(X)','e(X)','is_selected','equivalent_width']
    ll = ll[cols]
    tab = ll.group_by('species')
    for ttab in tab.groups:
        summary = summarize_abundances_species(ttab)
        elem_root = QtGui.QTreeWidgetItem(abundtree, summary)
        for row in ttab:
            data = [row['wavelength'],row['expot'],row['loggf'],row['A(X)'],row['e(X)'],row['equivalent_width']]
            line = QtGui.QTreeWidgetItem(elem_root, [str(x) for x in data])
            line.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
    
#    header = QtGui.QTreeWidgetItem(["Virtual folder tree","Comment"])
#    abundtree.setHeaderItem(header)   #Another alternative is setHeaderLabels(["Tree","First",...])
#    
#    root = QtGui.QTreeWidgetItem(abundtree, ["Untagged files"])
#    root.setData(2, QtCore.Qt.EditRole, 'Some hidden data here')# Data set to column 2, which is not visible
#    
#    folder1 = QtGui.QTreeWidgetItem(root, ["Interiors"])
#    folder1.setData(2, QtCore.Qt.EditRole, 'Some hidden data here')# Data set to column 2, which is not visible
#    
#    folder2 = QtGui.QTreeWidgetItem(root, ["Exteriors"])
#    folder2.setData(2, QtCore.Qt.EditRole, 'Some hidden data here')# Data set to column 2, which is not visible
#    
#    folder1_1 = QtGui.QTreeWidgetItem(folder1, ["Bathroom", "Seg was here"])
#    folder1_1.setData(2, QtCore.Qt.EditRole, 'Some hidden data here')# Data set to column 2, which is not visible
#    
#    folder1_2 = QtGui.QTreeWidgetItem(folder1, ["Living room", "Approved by client"])
#    folder1_2.setData(2, QtCore.Qt.EditRole, 'Some hidden data here')# Data set to column 2, which is not visible
#    
#    
#    
#    def printer( treeItem ):
#        foldername = treeItem.text(0)
#        comment = treeItem.text(1)
#        data = treeItem.text(2)
#        print(foldername + ': ' + comment + ' (' + data + ')')
#    
#    
#    abundtree.itemClicked.connect( lambda : printer( abundtree.currentItem() ) )
    
    
    abundtree.show()
    sys.exit(app.exec_())

def test():
#if __name__=="__main__":
    app = QtGui.QApplication(sys.argv)
    


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # init widgets
    view = QtGui.QTreeView()
    view.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
    model = QtGui.QStandardItemModel()
    model.setHorizontalHeaderLabels(['col1', 'col2', 'col3'])
    view.setModel(model)
    view.setUniformRowHeights(True)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # populate data
    for i in range(3):
        parent1 = QtGui.QStandardItem('Family {}. Some long status text for sp'.format(i))
        for j in range(3):
            child0 = QtGui.QStandardItem('Checkbox')
            child0.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
            child1 = QtGui.QStandardItem('Child {}'.format(i*3+j))
            child2 = QtGui.QStandardItem('row: {}, col: {}'.format(i, j+1))
            child3 = QtGui.QStandardItem('row: {}, col: {}'.format(i, j+2))
            parent1.appendRow([child0, child1, child2, child3])
        model.appendRow(parent1)
        # span container columns
        view.setFirstColumnSpanned(i, view.rootIndex(), True)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # expand third container
    index = model.indexFromItem(parent1)
    view.expand(index)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # select last row
    selmod = view.selectionModel()
    index2 = model.indexFromItem(child3)
    selmod.select(index2, QtGui.QItemSelectionModel.Select|QtGui.QItemSelectionModel.Rows)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    view.show()
    sys.exit(app.exec_())
