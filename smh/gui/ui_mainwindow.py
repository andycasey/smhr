#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The main GUI window for Spectroscopy Made Hard. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import sys
import logging
import os
from PySide import QtCore, QtGui
import yaml

# Import functionality related to each tab
import rv, normalization, summary, stellar_parameters, chemical_abundances

# Functions related to warnings and exceptions.
import exception

import smh
from linelist_manager import TransitionsDialog
from isotope_manager import IsotopeDialog

logger = logging.getLogger(__name__)

class Ui_MainWindow(QtGui.QMainWindow):
    """
    The main GUI window for Spectroscopy Made Hard.
    """

    def __init__(self, session_path=None, spectrum_filenames=None):
        super(Ui_MainWindow, self).__init__()

        self.unsaved_session_changes = False
        self.session = None
        self.session_path = None

        if session_path is not None and spectrum_filenames is not None:
            raise WTFError

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

        if spectrum_filenames is not None:
            print("DEBUGGING ONLY")
            self.new_session(spectrum_filenames)
            self.rv_tab.cross_correlate_and_correct()
            self.normalization_tab.normalize_and_stitch()


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
        file_menu.addSeparator()

        file_menu.addAction(open_session)
        self.open_recent_menu = QtGui.QMenu("Open &recent", self)


        # Read recently opened from the default settings path.
        with open(smh.Session._default_settings_path, "rb") as fp:
            recently_opened = yaml.load(fp).get("_gui_recently_opened", [])

        self.update_recently_opened(recently_opened)
        file_menu.addMenu(self.open_recent_menu)
        
        file_menu.addSeparator()
        file_menu.addAction(save_session)
        file_menu.addAction(save_session_as)

        self.action_transitions_manager = QtGui.QAction("&Transitions..", self,
            statusTip="Manage line lists and spectral models",
            triggered=self.transitions_manager)
        self.action_isotopes_manager = QtGui.QAction("&Isotopes..", self,
            statusTip="Manage isotopes for elements and molecules",
            triggered=self.isotopes_manager)
        self.action_transitions_manager.setEnabled(False)
        self.action_isotopes_manager.setEnabled(False)
        edit_menu = self.menuBar().addMenu("&Edit")
        edit_menu.addAction(self.action_transitions_manager)
        edit_menu.addAction(self.action_isotopes_manager)

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
        self._default_statusbar_message = "Spectroscopy Made Harder v{} ({})"\
            .format(smh.__version__, smh.__git_status__)
        self.statusbar.showMessage(self._default_statusbar_message)
        self.setStatusBar(self.statusbar)

        return True


    def add_to_recently_opened(self, path):
        """
        Add the specified path to the recently opened list.

        :param path:
            The path of the recently opened file to add to the list.
        """

        with open(smh.Session._default_settings_path, "rb") as fp:
            default_settings = yaml.load(fp)

        default_settings.setdefault("_gui_recently_opened", [])

        # Only show unique entries.
        try:
            default_settings["_gui_recently_opened"].remove(path)

        except ValueError:
            None

        default_settings["_gui_recently_opened"].insert(0, path)
        default_settings["_gui_recently_opened"] \
            = default_settings["_gui_recently_opened"][-5:]

        with open(smh.Session._default_settings_path, "wb") as fp:
            fp.write(yaml.dump(default_settings))

        self.update_recently_opened(default_settings["_gui_recently_opened"])
        return None


    def clear_recently_opened(self):
        """
        Clear the recently opened list.
        """

        with open(smh.Session._default_settings_path, "rb") as fp:
            default_settings = yaml.load(fp)

        default_settings["_gui_recently_opened"] = []
        with open(smh.Session._default_settings_path, "wb") as fp:
            fp.write(yaml.dump(default_settings))

        self.update_recently_opened([])
        return None


    def update_recently_opened(self, paths):
        """
        Update the recently opened menu with new entries.

        :param paths:
            The disk paths of the recently opened sessions.
        """

        self.open_recent_menu.clear()

        for i, path in enumerate(paths):
            action = QtGui.QAction(os.path.basename(path), self)
            self.connect(action, QtCore.SIGNAL("triggered()"),
                lambda *_: self.open_session(path))
            self.open_recent_menu.addAction(action)

        if len(paths) == 0:
            self._no_recent_sessions = QtGui.QAction(
                "(No recent sessions)", self, triggered=lambda *_: None)
            self._no_recent_sessions.setEnabled(False)
            self.open_recent_menu.addAction(self._no_recent_sessions)

        self.open_recent_menu.addSeparator()
        self.open_recent_menu.addAction(QtGui.QAction("Clear recently opened",
            self, statusTip="Clear the recently opened sessions",
            triggered=self.clear_recently_opened))

        return None


    def new_session(self, filenames=None):
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
        if filenames is None:
            filenames, selected_filter = QtGui.QFileDialog.getOpenFileNames(
                self, caption="Select input spectra", dir="")
            if not filenames:
                return None

        # Create a session.
        self.session = smh.Session(filenames)
        self.session_path = None

        # Import default session settings
        with open(smh.Session._default_settings_path, "rb") as fp:
            defaults = yaml.load(fp)
        self.session.metadata.update(defaults)

        # Disable all tabs except for Summary and RV.
        for i in range(self.tabs.count()):
            self.tabs.setTabEnabled(i, i < 2)

        # Enable relevant menu actions.
        self.action_transitions_manager.setEnabled(True)
        self.action_isotopes_manager.setEnabled(True)

        # Re-populate widgets in all tabs.
        self.summary_tab._populate_widgets()
        self.rv_tab.update_from_new_session()
        self.normalization_tab._populate_widgets()
        self.stellar_parameters_tab.populate_widgets()

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

        print("opening {}".format(path))


        if path is None:
            path, _ = QtGui.QFileDialog.getOpenFileName(self,
                caption="Select session", dir="", filter="*.smh")
            if not path: return


        self.add_to_recently_opened(path)
        self.session_path = path

        print("Opening session from {}".format(path))
        self.session = smh.Session.load(path)

        # Enable relevant menu actions.
        self.action_transitions_manager.setEnabled(True)
        self.action_isotopes_manager.setEnabled(True)

        # Refresh GUI of all tabs 
        # TODO does this all work?
        # Summary tab
        for i in range(self.tabs.count()):
            self.tabs.setTabEnabled(i, i <= 1)
        self.summary_tab._populate_widgets()

        if "rv" not in self.session.metadata: return None
        #self.rv_tab.new_session_loaded()
        # TODO put all these in RV tab
        self.rv_tab._populate_widgets()
        self.rv_tab.draw_template(refresh=True)
        self.rv_tab.redraw_ccf(refresh=True)
        if "rv_measured" not in self.session.metadata["rv"]: return None
        
        if "normalization" not in self.session.metadata: return None
        self.tabs.setTabEnabled(2, True)
        #self.normalization_tab.new_session_loaded()
        # TODO put all these in normalization tab
        self.normalization_tab._populate_widgets()
        #self.normalization_tab.draw_order()
        #self.normalization_tab.current_order_index = 0
        # TODO seems to not save the stitched normalized spectrum?
        self.normalization_tab.draw_continuum(refresh=True)
        if "continuum" not in self.session.metadata["normalization"]: return None

        self.tabs.setTabEnabled(3, True)
        self.tabs.setTabEnabled(4, True)
        #self.stellar_parameters_tab.new_session_loaded()
        # TODO there are likely more things needed here!
        self.stellar_parameters_tab.populate_widgets()

        self.chemical_abundances_tab.new_session_loaded()
        
        return None



    def save_session(self):
        """ Save session. """
        if self.session_path is None:
            self.save_session_as()
            return
        print("Saving to {}".format(self.session_path))
        self.session.save(self.session_path, overwrite=True)
        return None


    def save_session_as(self, path=None):
        """ Save session as new filename. """
        print("Save session as")
        if path is None:
            path, _ = QtGui.QFileDialog.getSaveFileName(self,
                caption="Enter filename", dir="", filter="*.smh")
            if not path: return
        self.session_path = path
        print("Saving to {}".format(self.session_path))
        self.session.save(path, overwrite=True)
        return None

    def export_normalized_spectrum(self):
        """ Export a normalized, rest-frame spectrum. """

        self.session.normalized_spectrum.write("test.txt")
        print("wrote to test.txt")


    def transitions_manager(self):
        """
        Open the transitions manager dialog to edit line lists and spectral
        models.
        """

        # Ensure to update the proxy data models when the transitions dialog has
        # been closed.
        window = TransitionsDialog(self.session, callbacks=[
            self.session.index_spectral_models,
            self.stellar_parameters_tab.proxy_spectral_models.reset,
            ])
        window.exec_()

        return None


    def isotopes_manager(self):
        """
        Open the isotopes manager dialog.
        """
        window = IsotopeDialog(self.session)
        window.exec_()
        return None


    def __init_ui__(self):
        """
        Set up the primary user interface (not the stuff in tabs).
        """
        
        # Create the central widget with a vertical layout.
        cw = QtGui.QWidget(self)
        cw_vbox = QtGui.QVBoxLayout(cw)
        cw_vbox.setContentsMargins(10, 20, 10, 10)



        # Create the primary widget for all the main tabs.
        self.tabs = QtGui.QTabWidget(cw)
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

        # Create chemical abundances tab
        # BUT IT'S XBOX HUGE
        self.chemical_abundances_tab \
            = chemical_abundances.ChemicalAbundancesTab(self)
        self.tabs.addTab(self.chemical_abundances_tab, "Chemical abundances")

        # Disable all tabs except the first one.
        for i in range(self.tabs.count()):
            self.tabs.setTabEnabled(i, i == 0)

        cw_vbox.addWidget(self.tabs)
        self.setCentralWidget(cw)

        self.tabs.setCurrentIndex(0)
        self._update_window_title()

