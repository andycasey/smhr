from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import os
import sys
from PySide2 import (QtCore, QtWidgets as QtGui)
from six import string_types, iteritems

from smh import (Session, specutils)
from smh.spectral_models import ProfileFittingModel, SpectralSynthesisModel
import smh.radiative_transfer as rt
from smh.linelists import LineList

from astropy import table

import time
import logging
logger = logging.getLogger(__name__)

_treecols = ["is_selected","wavelength","expot","loggf","element","A(X)","e(X)","equivalent_width","A(X)","A(X)"]
_treecolmap = dict(zip(range(len(_treecols)),_treecols))

def summarize_abundances_species(ttab,use_weights=True):
    element = ttab["element"][0]
    ttab = ttab[ttab["is_selected"]]
    N = len(ttab)
    if N==0: return [element, N, np.nan, np.nan, np.nan, np.nan]
    if use_weights:
        weights = 1./ttab["e(X)"]**2
    else:
        weights = np.ones(N)
    total_weights = np.sum(weights)
    abund = np.sum(ttab["A(X)"]*weights)/total_weights
    stdev = np.sum(ttab["e(X)"]*weights**2)/(total_weights**2)
    XH = np.nan
    XFe = np.nan
    return [element,N,abund,stdev,XH,XFe]
def summarize_abundances(tab):
    tab = tab.group_by("species")
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
        self.setUniformRowHeights(True)
        self.setItemDelegate(AbundTreeViewDelegate(self))
    def span_cols(self):
        """
        Have to call after connecting a model
        """
        for i in range(10):
            if i==0: self.setColumnWidth(i,70)
            else: self.setColumnWidth(i,80)
        for i in range(len(self.model().summaries)):
            self.setFirstColumnSpanned(i,None, True)
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
    def __init__(self, parent, sm_indices, ab_indices, index, model):
        self.sm_indices = sm_indices
        self.ab_indices = ab_indices
        self.index = index
        self.model = model
        super(AbundTreeElementSummaryItem, self).__init__(parent, index)
        
        self.compute_summary()
        self.fmts = ["{:5}", "N={:3}", "A(X)={:5.2f}", "e(X)={:5.2f}", "[X/H]={:5.2f}", "[X/Fe]={:5.2f}"]
        self.precols = ["","N=","A(X)=","e(X)=","[X/H]=","[X/Fe]="]
    def _getChildren(self):
        N = len(self.sm_indices)
        return [AbundTreeMeasurementItem(row,self.sm_indices[row],self.ab_indices[row],self) for row in range(N)]
    def columnCount(self):
        return 6 #1
    def data(self, column):
        if column >= self.columnCount(): return None
        return self.fmts[column].format(self.summary[column])
        #return self.print_summary()
    def print_summary(self):
        return "{0:5} N={1:3} A(X)={2:5.2f} e(X)={3:5.2f} [X/H]={4:5.2f} [X/Fe]={5:5.2f}".format(*self.summary)
    def compute_summary(self):
        assert len(self.subnodes) >= 1
        elem = self.subnodes[0].data(8)
        N = len(self.subnodes)
        abunds = np.zeros(N)*np.nan
        errs = np.zeros(N)*np.nan
        for i,node in enumerate(self.subnodes):
            if node.data(0):
                errs[i] = float(node.data(5))
                abunds[i] = float(node.data(2))
        weights = 1./errs**2
        # TODO errors must be > 0 right now
        ii = np.logical_and(~np.isnan(abunds), ~np.isnan(weights))
        abunds = abunds[ii]
        errs = errs[ii]
        weights = weights[ii]
        N = len(abunds)
        total_weights = np.sum(weights)
        abund = np.sum(abunds*weights)/total_weights
        stdev = np.sum(errs*weights**2)/(total_weights**2)
        # TODO [X/H]
        # TODO [X/Fe]
        self.summary = [elem,N,abund,stdev,np.nan,np.nan]
        return None
class AbundTreeMeasurementItem(AbundTreeItem):
    def __init__(self,row,sm_ix,ab_ix,parent):
        self.row = row
        self.sm_ix = sm_ix # which spectral model
        self.ab_ix = ab_ix # which element abundance in that spectral model
        assert isinstance(parent, AbundTreeElementSummaryItem)
        self.parent = parent
        super(AbundTreeMeasurementItem, self).__init__(parent, row)
    def _getChildren(self):
        return []
    def columnCount(self):
        return 10
    def data(self, column):
        m = self.parent.model.parenttab.spectral_models[self.sm_ix]
        if column==0:
            return m.is_acceptable
        if isinstance(m,ProfileFittingModel):
            if column==1: #wl
                return "{:6.1f}".format(m.transitions["wavelength"][0])
            elif column==2: #A(X)
                try:
                    return "{:5.2f}".format(m.metadata["fitted_result"][2]["abundances"][0])
                except:
                    return str(np.nan)
            elif column==3: #equivalent_width
                try:
                    return "{:6.2f}".format(m.metadata["fitted_result"][2]["equivalent_width"][0]*1000.)
                except:
                    return str(np.nan)
            elif column==4: #REW
                try:
                    return "{:6.2f}".format(np.log10(m.metadata["fitted_result"][2]["equivalent_width"][0]/m.transitions["wavelength"][0]))
                except:
                    return str(np.nan)
            elif column==5: #e(X)
                #TODO
                return "{:4.2f}".format(0.1)
            elif column==6: #EP
                return "{:4.2f}".format(m.transitions["expot"][0])
            elif column==7: #loggf
                return "{:7.3f}".format(m.transitions["loggf"][0])
            elif column==8: #element
                return m.transitions["element"][0]
            else:
                return None
        else: #SpectralSynthesisModel
            raise NotImplementedError
    
class AbundTreeModel(QtCore.QAbstractItemModel):
    def __init__(self, parenttab, *args):
        super(AbundTreeModel, self).__init__(parenttab, *args)
        #there were some display bugs when this was named "parent"!
        self.parenttab = parenttab 
        # Summaries and all_species should be the same size
        self.summaries = []
        self.all_species = []
        self.all_sm_indices = {}
        self.all_ab_indices = {}

    def session_updated(self):
        # TODO call this whenever the session updates measurements in any way
        self.obtain_measurements_from_parent_tab()
    def obtain_measurements_from_parent_tab(self):
        logger.debug("Sorting measurements..."); start = time.time()
        measurements = self.parenttab.spectral_models
        all_species = []
        all_sm_indices = {} #species to index list
        all_ab_indices = {} #species to index list
        for i,m in enumerate(measurements):
            if isinstance(m,ProfileFittingModel):
                line = m.transitions[0]
                species = line["species"]
                if species in all_sm_indices: 
                    all_sm_indices[species].append(i)
                    all_ab_indices[species].append(0)
                else: 
                    all_species.append(species)
                    all_sm_indices[species] = [i]
                    all_ab_indices[species] = [0]
            if isinstance(m,SpectralSynthesisModel):
                raise NotImplementedError
        self.all_species = np.sort(all_species)
        assert len(self.all_species) == len(np.unique(self.all_species))
        self.all_sm_indices = all_sm_indices
        self.all_ab_indices = all_ab_indices

    def _getSummaries(self):
        summaries = []
        for index, species in enumerate(self.all_species):
            # Initialize the abundance summary 
            sm_indices = self.all_sm_indices[species]
            ab_indices = self.all_ab_indices[species]
            elem_summary = AbundTreeElementSummaryItem(None, sm_indices, ab_indices, index, self)
            summaries.append(elem_summary)
        return summaries

    def index(self, row, column, parent):
        if not parent.isValid(): # Root node, get a summary
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
        # TODO figure out current open/closed/selected items in tree
        # repopulate after resetting the model
        
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
                checked = item.data(0)# == "True" #grrr
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
                #assert isinstance(value, bool),value
                m = self.parenttab.spectral_models[item.sm_ix]
                m.metadata["is_acceptable"] = value
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
            cols = ["","wl","A(X)","EW","REW","e(X)","EP","loggf","element",""]
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
    datadir = os.path.dirname(os.path.abspath(__file__))+"/../tests/test_data"
    session = Session([datadir+"/spectra/hd122563.fits"])
    session.metadata.update(defaults)
    ll = LineList.read(os.path.dirname(os.path.abspath(__file__))+"/../tests/test_data/linelists/complete.list")
    session.metadata["line_list"] = ll
    import cPickle as pickle
    print("Loading pre-saved spectral models"); start = time.time()
    with open(datadir+"/ewtest.pkl","rb") as fp:
        session.metadata["spectral_models"] = pickle.load(fp)
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

class AbundTreeViewDelegate(QtGui.QStyledItemDelegate):
    def paint(self, painter, option, index):
        if isinstance(index.internalPointer(), AbundTreeElementSummaryItem):
            option.font.setWeight(QtGui.QFont.Bold)
        super(AbundTreeViewDelegate, self).paint(painter, option, index)

