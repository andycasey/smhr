from PySide import QtCore, QtGui
import sys
import numpy as np
import pdb

class MyTableView(QtGui.QTableView):
    def __init__(self, parent=None, *args):
        super(MyTableView, self).__init__(parent, *args)
        self.parent = parent
        #self.setStyleSheet("QItem { font-size: 10px; }")
        self.verticalHeader().setResizeMode(QtGui.QHeaderView.Fixed)
        self.verticalHeader().setDefaultSectionSize(20)
        #self.resizeColumnsToContents()
        #self.resizeRowsToContents()
        #self.setStyleSheet("QTableView::item { border: 0px; }")

class MyCheckboxDelegate(QtGui.QStyledItemDelegate):
    """ http://stackoverflow.com/questions/3363190/qt-qtableview-how-to-have-a-checkbox-only-column """
    def createEditor(self, parent, option, index):
        return None
    

class MyTableModel(QtCore.QAbstractTableModel):
    def __init__(self, parent=None, *args):
        super(MyTableModel, self).__init__(parent, *args)
        self.parent = parent 
        self._data = np.arange(100).reshape(20,5)
        self._checked = np.zeros(20,dtype=bool)
        return None

    def data(self, index, role):
        if not index.isValid():
            return None
        row = index.row()
        column = index.column()
        if column == 0:
            if role == QtCore.Qt.CheckStateRole:
                return QtCore.Qt.Checked if self._checked[row] else QtCore.Qt.Unchecked
            elif role==QtCore.Qt.FontRole:
                return QtGui.QFont("Helvetica Neue", 10)
            else:
                return None
        column = column - 1
        if role==QtCore.Qt.FontRole:
            return QtGui.QFont("Helvetica Neue", 10)
        return self._data[row,column] if role == QtCore.Qt.DisplayRole else None
    def setData(self, index, value, role=QtCore.Qt.DisplayRole):
        print value
        if index.column() == 0:
            row = index.row()
            #value = (value != 0)
            self._checked[row] = value
            self.dataChanged.emit(index,index)
        else:
            try:
                self._data[index.row(), index.column()-1] = value
                value = True
            except:
                value = False
        return value
    def flags(self, index):
        if not index.isValid():
            return None
        if index.column() == 0:
            return  QtCore.Qt.ItemIsEditable|\
                    QtCore.Qt.ItemIsEnabled|\
                    QtCore.Qt.ItemIsUserCheckable
        else:
            return  QtCore.Qt.ItemIsSelectable|\
                    QtCore.Qt.ItemIsEnabled|\
                    QtCore.Qt.ItemIsEditable
    
    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal \
        and role == QtCore.Qt.DisplayRole:
            return str(col+1)
        return None
    def rowCount(self, parent):
        return self._data.shape[0]
    def columnCount(self, parent):
        return self._data.shape[1] + 1
    
if __name__=="__main__":
    app = QtGui.QApplication(sys.argv)
    if sys.platform == "darwin":
        # See http://successfulsoftware.net/2013/10/23/fixing-qt-4-for-mac-os-x-10-9-mavericks/
        substitutes = [
            (".Lucida Grande UI", "Lucida Grande"),
            (".Helvetica Neue DeskInterface", "Helvetica Neue")
        ]
        for substitute in substitutes:
            QtGui.QFont.insertSubstitution(*substitute)
    view = MyTableView()
    model = MyTableModel()
    view.setModel(model)
    view.show()
    sys.exit(app.exec_())
    
