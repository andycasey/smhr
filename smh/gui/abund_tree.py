from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import os
import sys
from PySide import QtCore, QtGui
from six import string_types

from smh import (Session, specutils)
import smh.spectral_models
import smh.radiative_transfer as rt
from smh.linelists import LineList

from astropy import table

import time
import logging
logger = logging.getLogger(__name__)

_treecols = ['is_selected','wavelength','expot','loggf','element','A(X)','e(X)','equivalent_width','A(X)','A(X)']
_treecolmap = dict(zip(range(len(_treecols)),_treecols))

def summarize_abundances_species(ttab,use_weights=True):
    element = ttab['element'][0]
    ttab = ttab[ttab['is_selected']]
    N = len(ttab)
    if N==0: return [element, N, np.nan, np.nan, np.nan, np.nan]
    if use_weights:
        weights = 1./ttab['e(X)']**2
    else:
        weights = np.ones(N)
    total_weights = np.sum(weights)
    abund = np.sum(ttab['A(X)']*weights)/total_weights
    stdev = np.sum(ttab['e(X)']*weights**2)/(total_weights**2)
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
    def __init__(self, parent, *args):
        super(AbundTreeView, self).__init__(parent, *args)
        self.parent = parent
        font = QtGui.QFont("Monospace")
        font.setStyleHint(QtGui.QFont.TypeWriter)
        self.setFont(font)
    def span_cols(self):
        """
        Have to call after connecting a model
        """
        for i in range(10):
            if i==0: self.setColumnWidth(i,70)
            else: self.setColumnWidth(i,70)
        for i in range(len(self.model().summaries)):
            self.setFirstColumnSpanned(i,self.rootIndex(), True)
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
        self.fmts = ["{:5}", "N={:3}", "A(X)={:5.2f}", "e(X)={:5.2f}", "[X/H]={:5.2f}", "[X/Fe]={:5.2f}"]
        self.precols = ['','N=','A(X)=','e(X)=','[X/H]=','[X/Fe]=']
    def _getChildren(self):
        N = len(self.model.tab.groups[self.index])
        return [AbundTreeMeasurementItem(row,self) for row in range(N)]
    def columnCount(self):
        return 6 #1
    def data(self, column):
        if column >= self.columnCount(): return None
        return self.fmts[column].format(self.summary[column])
        #return self.print_summary()
    def print_summary(self):
        return "{0:5} N={1:3} A(X)={2:5.2f} e(X)={3:5.2f} [X/H]={4:5.2f} [X/Fe]={5:5.2f}".format(*self.summary)
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
    def __init__(self, parenttab, *args):
        super(AbundTreeModel, self).__init__(parenttab, *args)
        #there were some display bugs when this was named "parent"!
        self.parenttab = parenttab 
        self.tab = None #self.obtain_measurements_from_session()
        # TODO set up a map from spectral models to items in the tree
        # That way you can selectively update the internal table here
        # TODO alternatively, let's just hook all the data into a spectral model list 
        self.summaries = self._getSummaries()

    def session_updated(self,spectral_model=None):
        # TODO call this whenever the session updates measurements in any way
        if spectral_model is None:
            self.obtain_measurements_from_parent_tab()
        else:
            # TODO get the item from the measurement
            # TODO update just those items with the new measurements
            raise NotImplementedError
    def obtain_measurements_from_parent_tab(self):
        print("Summarizing measurements"); start = time.time()
        measurements = self.parenttab.spectral_models
        wl = []
        EP = []
        loggf = []
        element = []
        species = []
        abund = []
        err = []
        is_selected = []
        EW = []
        for m in measurements:
            if isinstance(m,smh.spectral_models.ProfileFittingModel):
                try:
                    ab = m.abundances[0]
                except KeyError:
                    abund.append(np.nan)
                    err.append(np.nan)
                    EW.append(np.nan)
                    is_selected.append(False)
                except rt.RTError:
                    abund.append(np.nan)
                    err.append(np.nan)
                    EW.append(1000.*m.metadata["fitted_result"][2]["equivalent_width"][0])
                    is_selected.append(False)
                else:
                    abund.append(ab)
                    err.append(0.1) #TODO
                    EW.append(1000.*m.metadata["fitted_result"][2]["equivalent_width"][0])
                    is_selected.append(m.metadata['is_acceptable'])
                line = m.transitions[0]
                wl.append(line['wavelength'])
                EP.append(line['expot'])
                loggf.append(line['loggf'])
                element.append(line['element'])
                species.append(line['species'])
            if isinstance(m,smh.spectral_models.SpectralSynthesisModel):
                raise NotImplementedError
        tab = table.Table([wl,EP,loggf,element,species,abund,err,is_selected,EW],
                          names=['wavelength','expot','loggf','element','species','A(X)','e(X)','is_selected','equivalent_width'])
        tab = tab.group_by('species')
        print("Computed! {:.1f}s".format(time.time()-start))
        self.tab = tab

    def _getSummaries(self):
        summaries = []
        if self.tab is None: return summaries
        # Initialize summary objects for each element
        for index,ttab in enumerate(self.tab.groups):
            # Initialize the abundance summary 
            elem_summary = AbundTreeElementSummaryItem(None, self.tab, index, self)
            summaries.append(elem_summary)
        return summaries

    def index(self, row, column, parent):
        if not parent.isValid(): # Root node
            #print("Index (root node)",row,column,parent)
            return self.createIndex(row, column, self.summaries[row])
        parentNode = parent.internalPointer()
        #print("Index (sub node)",row,column,parent)
        return self.createIndex(row, column, parentNode.subnodes[row])
    def parent(self, index):
        if not index.isValid():
            #print("Creating parent (bad index)")
            return QtCore.QModelIndex()
        node = index.internalPointer()
        if node.parent is None:
            #print("Creating parent (no parent)",index,node)
            return QtCore.QModelIndex()
        else:
            #print("Creating parent (node has parent)",index,node)
            return self.createIndex(node.parent.row, 0, node.parent)
    def reset(self):
        #self.beginResetModel()
        self.session_updated()
        self.summaries = self._getSummaries()
        QtCore.QAbstractItemModel.reset(self)
        #self.endResetModel()
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
    import main_ui
    app.window = main_ui.Ui_MainWindow()
    for i in range(app.window.tabs.count()):
        app.window.tabs.setTabEnabled(i, True)

    app.window.show()

    import yaml
    with open(Session._default_settings_path, "rb") as fp:
        defaults = yaml.load(fp)
    datadir = os.path.dirname(os.path.abspath(__file__))+'/../tests/test_data'
    session = Session([datadir+"/spectra/hd122563.fits"])
    session.metadata.update(defaults)
    ll = LineList.read(os.path.dirname(os.path.abspath(__file__))+'/../tests/test_data/linelists/complete.list')
    session.metadata['line_list'] = ll
    import cPickle as pickle
    print("Loading pre-saved spectral models"); start = time.time()
    with open(datadir+'/ewtest.pkl','rb') as fp:
        session.metadata['spectral_models'] = pickle.load(fp)
    print("Done!",time.time()-start)

    app.window.session = session
    chemical_abundances_tab = app.window.chemical_abundances_tab
    abundtree = chemical_abundances_tab.abundtree
    model = abundtree.model()

    #abundtree = AbundTreeView(None)
    #model = AbundTreeModel(None)
    #abundtree.setModel(model)
    #abundtree.span_cols()
    #abundtree.setGeometry(900, 400, 900, 400)
    #abundtree.move(QtGui.QApplication.desktop().screen().rect().center() \
    #              - abundtree.rect().center())
    #sp = QtGui.QSizePolicy(
    #    QtGui.QSizePolicy.MinimumExpanding, 
    #    QtGui.QSizePolicy.MinimumExpanding)
    #sp.setHeightForWidth(abundtree.sizePolicy().hasHeightForWidth())
    #abundtree.setSizePolicy(sp)
    #abundtree.show()
    
    #chemical_abundances_tab.show()
    
    sys.exit(app.exec_())
