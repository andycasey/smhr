#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A dialog to manage the atomic physics and spectral models available in the
current session.
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
import numpy as np
import sys
from PySide import QtCore, QtGui

from smh.linelists import LineList
from smh.spectral_models import (ProfileFittingModel, SpectralSynthesisModel)

logger = logging.getLogger(__name__)

if sys.platform == "darwin":
        
    # See http://successfulsoftware.net/2013/10/23/fixing-qt-4-for-mac-os-x-10-9-mavericks/
    substitutes = [
        (".Lucida Grande UI", "Lucida Grande"),
        (".Helvetica Neue DeskInterface", "Helvetica Neue")
    ]
    for substitute in substitutes:
        QtGui.QFont.insertSubstitution(*substitute)


class LineListTableModel(QtCore.QAbstractTableModel):

    headers = [u"Wavelength\n(Å)", "Species\n", "EP\n(eV)", "log(gf)\n", "C6\n",
        "D0\n", "Comments\n"]
    columns = ["wavelength", "element", "expot", "loggf", "damp_vdw", "dissoc_E",
        "comments"]

    def __init__(self, parent, session, *args):
        """
        An abstract model for line lists.
        """
        super(LineListTableModel, self).__init__(parent, *args)
        self.session = session

    def rowCount(self, parent):
        try:
            N = len(self.session.metadata["line_list"])
        except:
            N = 0
        return N

    def columnCount(self, parent):
        return len(self.headers)

    def data(self, index, role):
        if not index.isValid() or role != QtCore.Qt.DisplayRole:
            return None

        column = self.columns[index.column()]
        value = self.session.metadata["line_list"][column][index.row()]
        if column not in ("element", "comments"):
            return "{:.3f}".format(value)
        return value

    def setData(self, index, value, role):

        column = self.columns[index.column()]
        self.session.metadata["line_list"][column][index.row()] = value
        return value


    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal \
        and role == QtCore.Qt.DisplayRole:
            return self.headers[col]
        return None


    def sort(self, column, order):

        if "line_list" not in self.session.metadata:
            return None

        self.emit(QtCore.SIGNAL("layoutAboutToBeChanged()"))

        self.session.metadata["line_list"].sort(self.columns[column])
        if order == QtCore.Qt.DescendingOrder:
            self.session.metadata["line_list"].reverse()

        self.dataChanged.emit(self.createIndex(0, 0),
            self.createIndex(self.rowCount(0), self.columnCount(0)))
        self.emit(QtCore.SIGNAL("layoutChanged()"))


    def flags(self, index):
        if not index.isValid():
            return None
        return  QtCore.Qt.ItemIsEnabled|\
                QtCore.Qt.ItemIsEditable|\
                QtCore.Qt.ItemIsSelectable


class LineListTableView(QtGui.QTableView):

    def __init__(self, parent, session, *args):
        super(LineListTableView, self).__init__(parent, *args)
        self.session = session
        self._parent = parent


    def contextMenuEvent(self, event):
        """
        Provide a context (right-click) menu for the line list table.

        :param event:
            The mouse event that triggered the menu.
        """
        
        menu = QtGui.QMenu(self)
        import_lines = menu.addAction("Import lines..")
        menu.addSeparator()
        import_profiles = menu.addAction("Import lines for profile models..")
        import_syntheses = menu.addAction("Import files as synthesis models..")
        menu.addSeparator()
        add_profiles_action = menu.addAction("Model with profiles")
        add_synth_action = menu.addAction("Model by synthesis")
        menu.addSeparator()
        delete_action = menu.addAction("Delete")

        any_selected = len(self.selectionModel().selectedRows()) > 0
        if not any_selected:
            add_profiles_action.setEnabled(False)
            add_synth_action.setEnabled(False)
            delete_action.setEnabled(False)

        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == import_lines:
            self.import_from_filename()

        elif action == add_profiles_action:
            self.add_selected_rows_as_profile_models()
            
        elif action == add_synth_action:
            self.add_selected_rows_as_synthesis_model()

        elif action == import_profiles:
            self.add_imported_lines_as_profile_models()

        elif action == import_syntheses:
            self.add_imported_lines_as_synthesis_model()

        elif action == delete_action:
            self.delete_selected_rows()

        return None


    def add_imported_lines_as_profile_models(self):
        """ Import line list data from a file and create profile models. """

        spectral_models_to_add = []
        transitions = self.import_from_filename()
        N = len(transitions)
        for index in range(N):
            spectral_models_to_add.append(
                ProfileFittingModel(self.session, transitions["hash"][[index]]))

        self.session.metadata.setdefault("spectral_models", [])
        self.session.metadata["spectral_models"].extend(spectral_models_to_add)
        
        # Update the spectral models abstract table model.
        self._parent.models_view.model().rowsInserted.emit(
            QtCore.QModelIndex(), 0, N-1)
        self._parent.models_view.model().dataChanged.emit(
            QtCore.QModelIndex(), QtCore.QModelIndex())
        self._parent.models_view.resizeColumnsToContents()


    def add_imported_lines_as_synthesis_model(self):
        """
        Import line list data from a file and create a single synthesis model.
        """

        transitions = self.import_from_filename()
        spectral_model = SpectralSynthesisModel(self.session, 
            transitions["hash"], transitions.unique_elements)

        self.session.metadata.setdefault("spectral_models", [])
        self.session.metadata["spectral_models"].append(spectral_model)
        
        # Update the spectral models abstract table model.
        self._parent.models_view.model().rowsInserted.emit(
            QtCore.QModelIndex(), 0, 0)
        self._parent.models_view.model().dataChanged.emit(
            QtCore.QModelIndex(), QtCore.QModelIndex())
        self._parent.models_view.resizeColumnsToContents()


    def add_selected_rows_as_profile_models(self):
        """ Add the selected rows as profile spectral models. """

        spectral_models_to_add = []
        for row in self.selectionModel().selectedRows():
            spectral_models_to_add.append(
                ProfileFittingModel(self.session,
                    self.session.metadata["line_list"]["hash"][[row.row()]]))

        self.session.metadata.setdefault("spectral_models", [])
        self.session.metadata["spectral_models"].extend(spectral_models_to_add)
        
        # Update the spectral models abstract table model.
        self._parent.models_view.model().rowsInserted.emit(
            QtCore.QModelIndex(), 0, len(spectral_models_to_add)-1)
        self._parent.models_view.model().dataChanged.emit(
            QtCore.QModelIndex(), QtCore.QModelIndex())
        self._parent.models_view.resizeColumnsToContents()

        return None


    def add_selected_rows_as_synthesis_model(self):
        """ Add the selected rows as a single spectral synthesis model. """

        row_indices = []
        for row in self.selectionModel().selectedRows():
            row_indices.append(row.row())
        row_indices = np.array(row_indices)

        # Which elements are contributing here?
        transitions = self.session.metadata["line_list"][row_indices]
        elements = transitions.unique_elements

        if len(elements) > 1:
            raise NotImplementedError("ask the user which element(s)")

        spectral_model = SpectralSynthesisModel(self.session,
            transitions["hash"], elements)
            
        self.session.metadata.setdefault("spectral_models", [])
        self.session.metadata["spectral_models"].append(spectral_model)

        # Update the spectral models abstract table model.
        self._parent.models_view.model().rowsInserted.emit(
            QtCore.QModelIndex(), 0, 0)
        self._parent.models_view.model().dataChanged.emit(
            QtCore.QModelIndex(), QtCore.QModelIndex())
        self._parent.models_view.resizeColumnsToContents()

        return None


    def delete_selected_rows(self):
        """ Delete the rows selected in the table. """

        N_skipped = 0
        mask = np.ones(len(self.session.metadata["line_list"]), dtype=bool)
        all_hashes_used = []
        for spectral_model in self.session.metadata.get("spectral_models", []):
            all_hashes_used.extend(spectral_model._transition_hashes)
        all_hashes_used = list(set(all_hashes_used))

        for row in self.selectionModel().selectedRows():

            # Is this selected transition used in any spectral models?
            t_hash = self.session.metadata["line_list"]["hash"][row.row()]
            if t_hash in all_hashes_used:
                N_skipped += 1
            else:
                mask[row.row()] = False

        # Anything to do?
        if N_skipped > 0:
            # Display a warning dialog.
            QtGui.QMessageBox.warning(self, "Lines cannot be deleted",
                "{} line(s) selected for deletion cannot be removed because "
                "they are associated with existing spectral models.\n\n"
                "Delete the associated spectral models before deleting "
                "the lines.".format(N_skipped))

        if np.all(mask):
            return None

        self.session.metadata["line_list"] \
            = self.session.metadata["line_list"][mask]

        # TODO: There *must* be a better way to do this..
        for i, index in enumerate(np.sort(np.where(~mask)[0])):
            self.model().rowsRemoved.emit(
                QtCore.QModelIndex(), index - i, index - i)

        self.clearSelection()

        return None
    

    def import_from_filename(self):
        """ Import atomic physics data from a line list file. """

        filenames, selected_filter = QtGui.QFileDialog.getOpenFileNames(self,
            caption="Select files", dir="")
        if not filenames:
            return None

        # Load from files.
        line_list = LineList.read(filenames[0])
        for filename in filenames[1:]:
            line_list = line_list.merge(LineList.read(filename), in_place=False)

        # Merge the line list with any existing line list in the session.
        if self.session.metadata.get("line_list", None) is None:
            self.session.metadata["line_list"] = line_list
            N = len(line_list)
        else:
            N = len(self.session.metadata["line_list"]) - len(line_list)
            self.session.metadata["line_list"] \
                = self.session.metadata["line_list"].merge(
                    line_list, in_place=False)
        self.model().rowsInserted.emit(QtCore.QModelIndex(), 0, N - 1)
        self.model().dataChanged.emit(
            QtCore.QModelIndex(), QtCore.QModelIndex())

        return line_list



class SpectralModelsTableModel(QtCore.QAbstractTableModel):

    headers = [u"Wavelength\n(Å)", "Elements\n", "Model type\n",
        "Use for\nstellar parameters", "Use for\nchemical abundances"]

    def __init__(self, parent, session, *args):
        """
        An abstract table model for spectral models.
        """
        super(SpectralModelsTableModel, self).__init__(parent, *args)
        self.session = session
        self._parent = parent


    def rowCount(self, parent):
        try:
            N = len(self.session.metadata["spectral_models"])
        except:
            N = 0
        return N


    def columnCount(self, parent):
        return len(self.headers)


    def data(self, index, role):

        if role == QtCore.Qt.CheckStateRole and index.isValid() \
        and index.column() in (3, 4):

            sm = self.session.metadata["spectral_models"][index.row()]
            if index.column() == 3:
                attr = "use_for_stellar_parameter_inference"
            else:
                attr = "use_for_stellar_composition_inference"
            return QtCore.Qt.Checked if getattr(sm, attr) else QtCore.Qt.Unchecked

        if role != QtCore.Qt.DisplayRole or not index.isValid():
            return None

        spectral_model = self.session.metadata["spectral_models"][index.row()]
        column = index.column()

        if column == 0: # Wavelength (approx.)
            return spectral_model._repr_wavelength

        elif column == 1: # Element(s).
            return spectral_model._repr_element

        elif column == 2: # Model type.

            if isinstance(spectral_model, SpectralSynthesisModel):
                return "Spectral synthesis"
            elif isinstance(spectral_model, ProfileFittingModel):
                return "Profile fitting"
            else:
                return "Unknown"
        else:
            return None


    def setData(self, index, value, role):
        try:
            a = {
                3: "use_for_stellar_parameter_inference",
                4: "use_for_stellar_composition_inference"
            }[index.column()]

        except KeyError:
            return False

        else:
            self.session.metadata["spectral_models"][index.row()].metadata[a] = value
            self.dataChanged.emit(index, index)
            
            return value
    

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal \
        and role == QtCore.Qt.DisplayRole:
            return self.headers[col]
        return None


    def sort(self, column, order):

        if "spectral_models" not in self.session.metadata:
            return None

        self.emit(QtCore.SIGNAL("layoutAboutToBeChanged()"))

        sorters = {
            0: lambda sm: sm.transitions["wavelength"].mean(),
            1: lambda sm: sm._repr_element,
            2: lambda sm: isinstance(sm, SpectralSynthesisModel),
            3: lambda sm: sm.use_for_stellar_parameter_inference,
            4: lambda sm: sm.use_for_stellar_composition_inference
        }

        self.session.metadata["spectral_models"].sort(key=sorters[column])

        if order == QtCore.Qt.DescendingOrder:
            self.session.metadata["spectral_models"].reverse()

        self.dataChanged.emit(self.createIndex(0, 0),
            self.createIndex(self.rowCount(0), self.columnCount(0)))
        self.emit(QtCore.SIGNAL("layoutChanged()"))


    def flags(self, index):
        if not index.isValid():
            return None

        if index.column() in (3, 4):
            return  QtCore.Qt.ItemIsEditable|\
                    QtCore.Qt.ItemIsEnabled|\
                    QtCore.Qt.ItemIsUserCheckable

        else:
            return  QtCore.Qt.ItemIsSelectable|\
                    QtCore.Qt.ItemIsEnabled
                    





class SpectralModelsTableView(QtGui.QTableView):

    def __init__(self, parent, session, *args):
        super(SpectralModelsTableView, self).__init__(parent, *args)
        self.session = session
        self._parent = parent


    def contextMenuEvent(self, event):
        """
        Provide a context (right-click) menu for the spectral models table.

        :param event:
            The mouse event that triggered the menu.
        """
        
        menu = QtGui.QMenu(self)
        delete_action = menu.addAction("Delete")

        any_selected = len(self.selectionModel().selectedRows()) > 0
        if not any_selected:
            delete_action.setEnabled(False)

        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == delete_action:
            self.delete_selected_rows()

        return None


    def delete_selected_rows(self):
        """ Delete the selected spectral models. """

        delete_indices = [_.row() for _ in self.selectionModel().selectedRows()]

        self.session.metadata["spectral_models"] = [sm \
            for i, sm in enumerate(self.session.metadata["spectral_models"]) \
                if i not in delete_indices]

        # TODO: There *must* be a better way to do this..
        for i, index in enumerate(delete_indices):
            self.model().rowsRemoved.emit(
                QtCore.QModelIndex(), index - i, index - i)

        self.clearSelection()
        return None



class TransitionsDialog(QtGui.QDialog):

    def __init__(self, session, *args):
        """
        Initialise a dialog to manage the transitions (atomic physics and
        spectral models) for the given session.

        :param session:
            The session that will be inspected for transitions.
        """

        super(TransitionsDialog, self).__init__(*args)

        self.session = session

        self.setGeometry(900, 400, 900, 400)
        self.move(QtGui.QApplication.desktop().screen().rect().center() \
            - self.rect().center())
        self.setWindowTitle("Manage transitions")

        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, 
            QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)

        parent_vbox = QtGui.QVBoxLayout(self)
        tabs = QtGui.QTabWidget(self)

        # Line list tab.
        self.linelist_tab = QtGui.QWidget()
        self.linelist_view = LineListTableView(self, session)
        self.linelist_view.setModel(
            LineListTableModel(self, session))
        self.linelist_view.setSelectionBehavior(
            QtGui.QAbstractItemView.SelectRows)
        self.linelist_view.setSortingEnabled(True)
        self.linelist_view.resizeColumnsToContents()

        QtGui.QVBoxLayout(self.linelist_tab).addWidget(self.linelist_view)
        tabs.addTab(self.linelist_tab, "Line list")

        self.models_tab = QtGui.QWidget()
        self.models_view = SpectralModelsTableView(self, session)
        self.models_view.setModel(
            SpectralModelsTableModel(self, session))
        self.models_view.setSelectionBehavior(
            QtGui.QAbstractItemView.SelectRows)
        self.models_view.setSortingEnabled(True)
        self.models_view.resizeColumnsToContents()

        QtGui.QVBoxLayout(self.models_tab).addWidget(self.models_view)
        tabs.addTab(self.models_tab, "Spectral models")

        parent_vbox.addWidget(tabs)

        # A horizontal line.
        hr = QtGui.QFrame(self)
        hr.setFrameShape(QtGui.QFrame.HLine)
        hr.setFrameShadow(QtGui.QFrame.Sunken)
        parent_vbox.addWidget(hr)

        # Buttons.
        hbox = QtGui.QHBoxLayout()
        btn_import = QtGui.QPushButton(self)
        btn_import.setText("Import..")
        hbox.addWidget(btn_import)

        btn_export = QtGui.QPushButton(self)
        btn_export.setText("Export..")
        hbox.addWidget(btn_export)

        # Spacer with a minimum width.
        hbox.addItem(QtGui.QSpacerItem(40, 20, 
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))

        btn_save_as_default = QtGui.QPushButton(self)
        btn_save_as_default.setText("Save as default")
        hbox.addWidget(btn_save_as_default)

        btn_ok = QtGui.QPushButton(self)
        btn_ok.setText("OK")
        btn_ok.setFocus()
        hbox.addWidget(btn_ok)

        parent_vbox.addLayout(hbox)

        # Connect the buttons.
        btn_import.clicked.connect(self.import_transitions)
        btn_export.clicked.connect(self.export_transitions)
        btn_save_as_default.clicked.connect(self.save_as_default)
        btn_ok.clicked.connect(self.close)

        return None


    def import_transitions(self):
        """ Import transitions (line lists and spectral models) from a file. """
        raise NotImplementedError


    def export_transitions(self):
        """ Export transitions (line lists and spectral models) to a file. """
        raise NotImplementedError


    def save_as_default(self):
        """
        Save the current line list and spectral models as the defaults for
        future SMH sessions.
        """
        raise NotImplementedError



if __name__ == "__main__":

    # This is just for development testing.
    app = QtGui.QApplication(sys.argv)
    window = TransitionsDialog(None)
    window.exec_()

    

