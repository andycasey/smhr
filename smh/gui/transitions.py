#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A widget for fitting spectral models. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)


import operator
from PySide import QtCore, QtGui

import mpl


class SpectralModelsTableModel(QtCore.QAbstractTableModel):

    header = [" ", "Wavelength", "Element"]
    attrs = ("is_acceptable", "_repr_wavelength", "_repr_element")

    def __init__(self, parent, spectral_models, *args):
        super(SpectralModelsTableModel, self).__init__(parent, *args)
        self.spectral_models = spectral_models



    def row_selected(self, *args):
        print("selected ", args)

    def rowCount(self, parent):
        return len(self.spectral_models)

    def columnCount(self, parent):
        return len(self.header)

    def data(self, index, role):
        if not index.isValid():
            return None

        value = getattr(
            self.spectral_models[index.row()],
            self.attrs[index.column()])

        if index.column() == 0:
            if role == QtCore.Qt.CheckStateRole:
                return QtCore.Qt.Checked if value else QtCore.Qt.Unchecked
            else:
                return None

        elif role != QtCore.Qt.DisplayRole:
            return None

        return getattr(self.spectral_models[index.row()],
            self.attrs[index.column()])

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal \
        and role == QtCore.Qt.DisplayRole:
            return self.header[col]
        return None


    def setData(self, index, value, role=QtCore.Qt.DisplayRole):
        attr = self.attrs[index.column()]
        print(attr)
        if attr != "is_acceptable":
            return False

        self.spectral_models[index.row()].metadata[attr] = value
        self.dataChanged.emit(index, index)
        return value


    def sort(self, column, order):

        self.emit(QtCore.SIGNAL("layoutAboutToBeChanged()"))
        self.spectral_models = sorted(self.spectral_models,
            key=lambda sm: getattr(sm, self.attrs[column]))
        
        if order == QtCore.Qt.DescendingOrder:
            self.spectral_models.reverse()

        self.dataChanged.emit(self.createIndex(0, 0), self.createIndex(self.rowCount(0), self.columnCount(0)))
        self.emit(QtCore.SIGNAL("layoutChanged()"))


    def flags(self, index):
        if not index.isValid(): return
        return QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsUserCheckable




if __name__ == "__main__":

    class MyWindow(QtGui.QWidget):
        def __init__(self, data_list, *args):
            super(MyWindow, self).__init__(*args)

            # setGeometry(x_pos, y_pos, width, height)
            self.setGeometry(300, 200, 570, 450)
            self.setWindowTitle("Click on column title to sort")

            sp = QtGui.QSizePolicy(
                QtGui.QSizePolicy.MinimumExpanding, 
                QtGui.QSizePolicy.MinimumExpanding)
            sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
            self.setSizePolicy(sp)

            table_model = SpectralModelsTableModel(self, data_list)
            table_view = QtGui.QTableView()
            table_view.setModel(table_model)
            table_view.resizeColumnsToContents()

            table_view.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
            table_view.setSortingEnabled(True)

            selectionModel = table_view.selectionModel()
            selectionModel.selectionChanged.connect(self.row_selected)

            self.table_view = table_view

            layout = QtGui.QHBoxLayout(self)
            layout.setContentsMargins(10, 10, 10, 10)
            layout.addWidget(table_view)

            # MPL figure


            # Create a matplotlib widget.
            blank_widget = QtGui.QWidget(self)
            sp = QtGui.QSizePolicy(
                QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
            sp.setHorizontalStretch(0)
            sp.setVerticalStretch(0)
            sp.setHeightForWidth(blank_widget.sizePolicy().hasHeightForWidth())
            blank_widget.setSizePolicy(sp)
            blank_widget.setObjectName("norm_plot")


            self.norm_plot = mpl.MPLWidget(blank_widget, tight_layout=True,
                autofocus=True)

            mpl_layout = QtGui.QVBoxLayout(blank_widget)
            mpl_layout.addWidget(self.norm_plot, 1)
            layout.addWidget(blank_widget)

            self.setLayout(layout)


            self.mpl_axis = self.norm_plot.figure.add_subplot(111)
            self.mpl_axis.scatter([0, 1], [0.4, 0.6])

            self.norm_plot.draw()



        def row_selected(self, *args):
            indexes = self.table_view.selectionModel().selectedRows()
            for index in indexes:
                print("row {} selected".format(index.row()))

    import sys

    from smh import linelists, Session
    transitions = linelists.LineList.read("/Users/arc/research/ges/linelist/vilnius.ew")

    session = Session([
        "/Users/arc/codes/smh/hd44007red_multi.fits",
        "/Users/arc/codes/smh/hd44007blue_multi.fits",
    ])

    from smh import spectral_models as sm

    data_list = []
    for i in range(len(transitions)):
        if i % 2:
            data_list.append(sm.ProfileFittingModel(transitions[[i]], session))
        else:
            data_list.append(sm.SpectralSynthesisModel(transitions[[i]],
                session, transitions["elem1"][i]))

    app = QtGui.QApplication(sys.argv)
    window = MyWindow(data_list)
    window.show()
    #sys.exit(app.exec_())
    app.exec_()
