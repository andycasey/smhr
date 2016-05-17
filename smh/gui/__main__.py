#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The main GUI window for Spectroscopy Made Hard. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

from PySide import QtCore, QtGui

# Import functionality related to each tab
import rv, normalization, summary, stellar_parameters

import smh


class Ui_MainWindow(QtGui.QMainWindow):
    """
    The main GUI window for Spectroscopy Made Hard.
    """

    def __init__(self, session_path=None):
        super(Ui_MainWindow, self).__init__()

        self.unsaved_session_changes = False
        self.session = None

        # Load a session already?
        if session_path is not None:
            self.open_session(session_path)

        self.setObjectName("smh")
        self.resize(1200, 600)
        self.move(QtGui.QApplication.desktop().screen().rect().center() \
            - self.rect().center())

        # Initialise the menus and associated actions.
        self.__init_menus__()

        # Set up the UI.
        self.__init_ui__()



    def __init_menus__(self):
        """
        Initialize main window menus and associated actions.
        """

        # File menu.
        new_session = QtGui.QAction("&New", self,
            shortcut=QtGui.QKeySequence.New,
            statusTip="Create a new session",
            triggered=self.new_session)

        open_session = QtGui.QAction("&Open...", self,
            shortcut=QtGui.QKeySequence.Open,
            statusTip="Open an existing session from disk",
            triggered=self.open_session)

        save_session = QtGui.QAction("&Save", self,
            shortcut=QtGui.QKeySequence.Save,
            statusTip="Save the session to disk",
            triggered=self.save_session)

        save_session_as = QtGui.QAction("Save &As", self,
            statusTip="Save the session to a new file",
            triggered=self.save_session_as)

        file_menu = self.menuBar().addMenu("&File")
        file_menu.addAction(new_session)
        file_menu.addAction(open_session)
        file_menu.addAction(save_session)
        file_menu.addAction(save_session_as)

        # Export menu.
        self._menu_export_normalized_spectrum \
            = QtGui.QAction("Normalized rest-frame spectrum", self,
                statusTip="Export a normalized, rest-frame spectrum resampled "
                          "onto a common wavelength mapping",
                triggered=self.export_normalized_spectrum)
        self._menu_export_normalized_spectrum.setEnabled(False)
        export_menu = self.menuBar().addMenu("&Export")
        export_menu.addAction(self._menu_export_normalized_spectrum)

        self.statusbar = QtGui.QStatusBar(self)
        self.statusbar.setObjectName("statusbar")
        self.statusbar.showMessage("Spectroscopy Made Harder v{0} ({1})".format(
            smh.__version__, smh.__git_status__))
        self.setStatusBar(self.statusbar)

        return True


    def new_session(self):
        """ Initialise new session. """

        # Do we already have a session open with unsaved changes?
        if self.session is not None and self.unsaved_session_changes:
            response = QtGui.QMessageBox.question(self, "Are you sure?",
                "You have unsaved changes.\n\n"
                "Are you sure you want to start a new session?", 
                QtGui.QMessageBox.StandardButton.Yes \
                | QtGui.QMessageBox.StandardButton.No)

            if not response == QtGui.QMessageBox.Yes:
                return

        # Get filenames of input spectra.
        filenames, selected_filter = QtGui.QFileDialog.getOpenFileNames(self,
            caption="Select input spectra", dir="")
        if not filenames: return

        # Create a session.
        self.session = smh.Session(filenames)

        # Disable all tabs except for Summary and RV.
        for i in range(self.tabs.count()):
            self.tabs.setTabEnabled(i, i < 2)

        # Re-populate widgets in all tabs.
        self.summary_tab._populate_widgets()
        self.rv_tab.update_from_new_session()
        self.normalization_tab._populate_widgets()

        self._update_window_title()

        return None


    def _update_window_title(self):
        """
        Update the window title.
        """

        joiner, prefix = (" - ", "Spectroscopy Made Hard")
        if self.session is None:
            title = prefix
            
        else:
            try:
                object_name = self.session.metadata["OBJECT"]

            except (AttributeError, KeyError):
                title = joiner.join([prefix, "Unnamed"])

            else:
                title = joiner.join([prefix, object_name])

        self.setWindowTitle(title)

        return None


    def open_session(self, path=None):
        """ Open existing session. """

        if path is None:
            path, _ = QtGui.QFileDialog.getOpenFileName(self,
                caption="Select session", dir="", filter="*.smh")
            if not path: return

        raise NotImplementedError

        print("Open session")
        return None



    def save_session(self):
        """ Save session. """
        print("Save session")
        return None


    def save_session_as(self):
        """ Save session as new filename. """
        print("Save session as")
        return None


    def export_normalized_spectrum(self):
        """ Export a normalized, rest-frame spectrum. """

        self.session.normalized_spectrum.write("test.txt")
        print("wrote to test.txt")



    def __init_ui__(self):
        """
        Set up the primary user interface (not the stuff in tabs).
        """
        
        # Create the central widget with a vertical layout.
        cw = QtGui.QWidget(self)
        cw_vbox = QtGui.QVBoxLayout(cw)

        # Create an empty frame for padding at the top of the application.
        top_frame_pad = QtGui.QFrame(cw)
        top_frame_pad.setMinimumSize(QtCore.QSize(0, 10))
        top_frame_pad.setFrameShape(QtGui.QFrame.NoFrame)
        top_frame_pad.setFrameShadow(QtGui.QFrame.Plain)
        top_frame_pad.setLineWidth(0)

        cw_vbox.addWidget(top_frame_pad)


        # Create the primary widget for all the main tabs.
        self.tabs = QtGui.QTabWidget(cw)
        # TODO: review whether this is necessary.
        #self.tabs.setMinimumSize(QtCore.QSize(300, 0))
        self.tabs.setTabPosition(QtGui.QTabWidget.North)
        self.tabs.setUsesScrollButtons(False)

        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.MinimumExpanding)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(self.tabs.sizePolicy().hasHeightForWidth())
        self.tabs.setSizePolicy(sp)

        # Create summary tab.
        self.summary_tab = summary.SummaryTab(self)
        self.tabs.addTab(self.summary_tab, "Summary")

        # Create radial velocity tab
        self.rv_tab = rv.RVTab(self)
        self.tabs.addTab(self.rv_tab, "Radial velocity")
        
        # Create normalization tab.
        self.normalization_tab = normalization.NormalizationTab(self)
        self.tabs.addTab(self.normalization_tab, "Normalization")

        # Create stellar parameters tab.
        self.stellar_parameters_tab \
            = stellar_parameters.StellarParametersTab(self)
        self.tabs.addTab(self.stellar_parameters_tab, "Stellar parameters")

        # Add remaining empty tabs.
        extra_tab_names = ("Chemical abundances", )

        for tab_name in extra_tab_names:
            tab = QtGui.QWidget()
            self.tabs.addTab(tab, tab_name)
        
        # Disable all tabs except the first one.
        for i in range(self.tabs.count()):
            self.tabs.setTabEnabled(i, i == 0)

        cw_vbox.addWidget(self.tabs)
        self.setCentralWidget(cw)

        self.tabs.setCurrentIndex(0)
        #QtCore.QMetaObject.connectSlotsByName(self)

        self._update_window_title()


if __name__ == '__main__':

    import sys

    if sys.platform == "darwin":
            
        # See http://successfulsoftware.net/2013/10/23/fixing-qt-4-for-mac-os-x-10-9-mavericks/
        substitutes = [
            (".Lucida Grande UI", "Lucida Grande"),
            (".Helvetica Neue DeskInterface", "Helvetica Neue")
        ]
        for substitute in substitutes:
            QtGui.QFont.insertSubstitution(*substitute)

    app = QtGui.QApplication(sys.argv)
    window = Ui_MainWindow()
    window.show()
    sys.exit(app.exec_())
