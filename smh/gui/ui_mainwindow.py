#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The main GUI window for Spectroscopy Made Hard. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import sys
import logging
import os
from PySide2 import (QtCore, QtGui as QtGui2, QtWidgets as QtGui)
import yaml
import numpy as np

# Import functionality related to each tab
import rv, normalization, summary, stellar_parameters, chemical_abundances, review

import smh
#from balmer import BalmerLineFittingDialog
from balmer import *
from linelist_manager import TransitionsDialog
from isotope_manager import IsotopeDialog
from plotting import SummaryPlotDialog, SNRPlotDialog

logger = logging.getLogger(__name__)

class Ui_MainWindow(QtGui.QMainWindow):
    """
    The main GUI window for Spectroscopy Made Hard.
    """

    def __init__(self, parent=None, session_path=None, spectrum_filenames=None):
        super(Ui_MainWindow, self).__init__()

        self.parent = parent
        
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
        desktop = QtGui.QApplication.desktop()
        self.move(desktop.screen().rect().center() - self.rect().center())

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
            shortcut=QtGui2.QKeySequence.New,
            statusTip="Create a new session",
            triggered=self.new_session)

        open_session = QtGui.QAction("&Open...", self,
            shortcut=QtGui2.QKeySequence.Open,
            statusTip="Open an existing session from disk",
            triggered=self.open_session)

        save_session = QtGui.QAction("&Save", self,
            shortcut=QtGui2.QKeySequence.Save,
            statusTip="Save the session to disk",
            triggered=self.save_session)

        save_session_as = QtGui.QAction("Save &As", self,
            statusTip="Save the session to a new file",
            triggered=self.save_session_as)

        file_menu = self.menuBar().addMenu("&File")
        file_menu.addAction(new_session)
        file_menu.addSeparator()

        file_menu.addAction(open_session)
        
        file_menu.addSeparator()
        file_menu.addAction(save_session)
        file_menu.addAction(save_session_as)

        self.action_transitions_manager = QtGui.QAction("&Transitions..", self,
            statusTip="Manage line lists and spectral models",
            triggered=self.transitions_manager)
        self.action_isotopes_manager = QtGui.QAction("&Isotopes..", self,
            statusTip="Manage isotopes for elements and molecules",
            triggered=self.isotopes_manager)
        self.action_comparison_spectrum = QtGui.QAction("&Comparison Spectrum..", self,
            statusTip="Plot a comparison spectrum (in cyan)",
            triggered=self.comparison_spectrum_dialog)
        self.action_clear_comparison_spectrum = QtGui.QAction("&Clear Comparison Spectrum", self,
            statusTip="Remove comparison spectrum",
            triggered=self.clear_comparison_spectrum)
        self.action_transitions_manager.setEnabled(False)
        self.action_isotopes_manager.setEnabled(False)
        self.action_comparison_spectrum.setEnabled(True)
        self.action_clear_comparison_spectrum.setEnabled(True)
        edit_menu = self.menuBar().addMenu("&Edit")
        edit_menu.addAction(self.action_transitions_manager)
        edit_menu.addAction(self.action_isotopes_manager)
        edit_menu.addAction(self.action_comparison_spectrum)
        edit_menu.addAction(self.action_clear_comparison_spectrum)

        # Advanced menu
        advanced_menu = self.menuBar().addMenu("&Advanced")
        self._action_fit_balmer_lines = QtGui.QAction(
            "Fit Balmer line(s)", self, enabled=False,
            statusTip="Interactively fit Balmer line wing(s)",
            triggered=self.show_balmer_line_dialog)
        advanced_menu.addAction(self._action_fit_balmer_lines)


        # Plot menu
        plot_menu = self.menuBar().addMenu("&Plot")
        ncap_summary_plot = QtGui.QAction("&SNR Plot", self,
            statusTip="Make SNR plot",
            triggered=self.snr_plot)
        plot_menu.addAction(ncap_summary_plot)
        summary_plot = QtGui.QAction("&Summary Plot", self,
            statusTip="Make summary plot",
            triggered=self.summary_plot)
        plot_menu.addAction(summary_plot)
        ncap_summary_plot = QtGui.QAction("&Ncap Summary Plot", self,
            statusTip="Make ncap summary plot",
            triggered=self.ncap_summary_plot)
        plot_menu.addAction(ncap_summary_plot)

        # Export menu.
        self._menu_export_normalized_spectrum \
            = QtGui.QAction("Normalized rest-frame spectrum", self,
                statusTip="Export a normalized, rest-frame spectrum resampled "
                          "onto a common wavelength mapping",
                triggered=self.export_normalized_spectrum)
        self._menu_export_unnormalized_spectrum \
            = QtGui.QAction("Unnormalized rest-frame spectrum", self,
                statusTip="Export a coadded rest-frame spectrum resampled "
                          "onto a common wavelength mapping",
                triggered=self.export_unnormalized_spectrum)
        self._menu_export_stitched_continuum \
            = QtGui.QAction("Stitched continuum", self,
                statusTip="Export a coadded continuum",
                triggered=self.export_stitched_continuum)
        #self._menu_export_normalized_spectrum.setEnabled(False)
        self._menu_print_abundance_table \
            = QtGui.QAction("Print abundance table", self,
                statusTip="",
                triggered=self.print_abundance_table)
        self._menu_export_abundance_table \
            = QtGui.QAction("Export abundance table", self,
                statusTip="",
                triggered=self.export_abundance_table)
        self._menu_export_spectral_model_measurements \
            = QtGui.QAction("Export spectral model measurements", self,
                statusTip="",
                triggered=self.export_spectral_model_measurements)
        export_menu = self.menuBar().addMenu("&Export")
        export_menu.addAction(self._menu_export_normalized_spectrum)
        export_menu.addAction(self._menu_export_unnormalized_spectrum)
        export_menu.addAction(self._menu_export_stitched_continuum)
        export_menu.addAction(self._menu_print_abundance_table)
        export_menu.addAction(self._menu_export_abundance_table)
        export_menu.addAction(self._menu_export_spectral_model_measurements)

        self.statusbar = QtGui.QStatusBar(self)
        self._default_statusbar_message = "Spectroscopy Made Harder v{} ({})"\
            .format(smh.__version__, smh.__git_status__)
        self.statusbar.showMessage(self._default_statusbar_message)
        self.setStatusBar(self.statusbar)

        return True



    def closeEvent(self, event):
        """
        Check to see if the user wants to save their session before closing.

        :param event:
            The close event.
        """
        ## TODO in PyQt5 something causes this to display twice sometimes???

        if self.session is None:
            # There is no session.
            event.accept()
            return None

        reply = QtGui.QMessageBox.question(self, "Message", 
            "Any unsaved changes will be lost. Are you sure you want to quit?", 
            QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

        return None




    def new_session(self, filenames=None):
        """
        Initialise a new session.

        :param filenames: [optional]
            A list of input spectra.
        """

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
        print("filenames:",filenames)
        if filenames is None or filenames is False:
            filenames, selected_filter = QtGui.QFileDialog.getOpenFileNames(
                self, caption="Select input spectra", directory="", filter="*")
            print("filenames:",filenames)
            if not filenames:
                return None

        # Create a session.
        self.session = smh.Session(filenames)
        self.session_path = None

        # Import default session settings
        with open(smh.Session._default_settings_path, "rb") as fp:
            try:
                defaults = yaml.load(fp, yaml.FullLoader)
            except AttributeError:
                defaults = yaml.load(fp)
        self.session.metadata.update(defaults)

        # TODO: WE SHOULD REMOVE THIS: THE GUI SHOULD READ FROM .SETTINGS()
        

        # Disable all tabs except for Summary and RV.
        for i in range(self.tabs.count()):
            self.tabs.setTabEnabled(i, i < 2)

        # Enable/disable relevant menu actions.
        self.action_transitions_manager.setEnabled(True)
        self.action_isotopes_manager.setEnabled(True)
        self.action_comparison_spectrum.setEnabled(True)
        self.action_clear_comparison_spectrum.setEnabled(True)
        self._action_fit_balmer_lines.setEnabled(False)

        # Re-populate widgets in all tabs.
        self.summary_tab._populate_widgets()
        self.rv_tab.update_from_new_session()
        self.normalization_tab._populate_widgets()
        self.stellar_parameters_tab.new_session_loaded()
        self.chemical_abundances_tab.new_session_loaded()
        self.review_tab.new_session_loaded()

        self._update_window_title(os.path.basename(filenames[0]))

        
        return None


    def _update_window_title(self, default_title="Unnamed"):
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
                title = joiner.join([prefix, default_title])

            else:
                title = joiner.join([prefix, object_name])

        self.setWindowTitle(title)

        return None


    def open_session(self, path=None):
        """ Open existing session. """

        print("Testing path:",path)
        if path is None or path is False:
            path, _ = QtGui.QFileDialog.getOpenFileName(self,
                caption="Select session", directory="", filter="*.smh")
            if not path: return
        print("We got:",path)


        self.session_path = path

        self.session = smh.Session.load(path)

        # Enable relevant menu actions.
        self.action_transitions_manager.setEnabled(True)
        self.action_isotopes_manager.setEnabled(True)
        self.action_comparison_spectrum.setEnabled(True)
        self.action_clear_comparison_spectrum.setEnabled(True)
        self._action_fit_balmer_lines.setEnabled(False)


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
        if "rv_measured" not in self.session.metadata["rv"] \
        and "rv_applied" not in self.session.metadata["rv"]: return None
        
        self.normalization_tab._populate_widgets()
        print("populated normalization widgets")
        if "normalization" not in self.session.metadata: return None
        self.tabs.setTabEnabled(2, True)
        #self.normalization_tab.new_session_loaded()
        # TODO put all these in normalization tab
        #self.normalization_tab.draw_order()
        #self.normalization_tab.current_order_index = 0
        # TODO seems to not save the stitched normalized spectrum?
        self.normalization_tab.draw_continuum(refresh=True)
        if "continuum" not in self.session.metadata["normalization"]: return None

        self.tabs.setTabEnabled(3, True)
        self.tabs.setTabEnabled(4, True)
        self.tabs.setTabEnabled(5, True)
        
        self.stellar_parameters_tab.new_session_loaded()
        self.chemical_abundances_tab.new_session_loaded()
        self.review_tab.new_session_loaded()
        
        self._update_window_title(os.path.basename(self.session_path))

        # Bring to top after loading
        self.raise_()
        self.activateWindow()
        self.showNormal()
        return None



    def save_session(self):
        """ Save session. """

        if self.session_path is None:
            self.save_session_as()
            return

        logger.info("Saving to {}".format(self.session_path))
        self.session.save(self.session_path, overwrite=True)
        return None


    def save_session_as(self, path=None):
        """
        Save session as new filename.

        :param path: [optional]
            The filename where to save the session.
        """

        if path is None or path is False:
            path, _ = QtGui.QFileDialog.getSaveFileName(self,
                caption="Enter filename", directory="", filter="*.smh")
            if not path: return

        self.session_path = path
        self.save_session()
        return None


    def show_balmer_line_dialog(self):
        """ Show an interactive dialog for fitting Balmer line profiles. """

        try:
            v = self.session.metadata["rv"]["rv_applied"]
        except (KeyError, TypeError):
            v = 0

        # Create a copy of the input spectra, because we will shift them.
        spectra = []
        for spectrum in self.session.input_spectra:
            s = spectrum.copy()
            s.redshift(v=v)
            spectra.append(s)

        labels = ["Order {}".format(i) for i in range(1, 1 + len(spectra))]
        
        # Normalized spectrum?
        if hasattr(self.session, "normalized_spectrum"):
            spectra += [self.session.normalized_spectrum]
            labels += ["Normalized rest-frame spectrum"]

        dialog = BalmerLineFittingDialog(spectra,
            observed_spectra_labels=labels, session=self.session)
        dialog.exec_()

        return None


    def export_normalized_spectrum(self):
        """ Export a normalized, rest-frame spectrum. """
        if self.session is None: return
        path, _ = QtGui.QFileDialog.getSaveFileName(self,
            caption="Enter normalized rest frame spectrum filename", directory="", filter="")
        if not path: return
        self.session.export_normalized_spectrum(path)

    def export_unnormalized_spectrum(self):
        """ Export a stitched, unnormalized, rest-frame spectrum. """
        if self.session is None: return
        path, _ = QtGui.QFileDialog.getSaveFileName(self,
            caption="Enter unnormalized rest frame spectrum filename", directory="", filter="")
        if not path: return
        self.session.export_unnormalized_spectrum(path)

    def export_stitched_continuum(self):
        """ Export a stitched continuum. """
        if self.session is None: return
        path, _ = QtGui.QFileDialog.getSaveFileName(self,
            caption="Enter continuum filename", directory="", filter="")
        if not path: return
        self.session.export_stitched_continuum(path)

    def print_abundance_table(self):
        """ Print abundance table to console (HACK) """
        summary_dict = self.session.summarize_spectral_models()
        keys = summary_dict.keys()
        keys = np.sort(keys)
        print("species   N logeps  err [X/H] [X/Fe]")
        for key in keys:
            num_models, logeps, stdev, stderr, XH, XFe = summary_dict[key]
            print("{:7.1f} {:3} {:6.2f} {:4.2f} {:5.2f} {:6.2f}".format(key,num_models,logeps,stdev,XH,XFe))

    
    def export_abundance_table(self):
        if self.session is None: return
        path, _ = QtGui.QFileDialog.getSaveFileName(self,
            caption="Enter abundance table filename", directory="", filter="")
        if not path: return
        self.session.export_abundance_table(path)
    
    def export_spectral_model_measurements(self):
        if self.session is None: return
        path, _ = QtGui.QFileDialog.getSaveFileName(self,
            caption="Enter spectral model measurements filename", directory="", filter="")
        if not path: return
        self.session.export_spectral_model_measurements(path)
    
    def transitions_manager(self):
        """
        Open the transitions manager dialog to edit line lists and spectral
        models.
        """

        # Ensure to update the proxy data models when the transitions dialog has
        # been closed.
        window = TransitionsDialog(self.session, callbacks=[
            self.stellar_parameters_tab.new_session_loaded,
            self.chemical_abundances_tab.new_session_loaded,
            self.review_tab.new_session_loaded
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


    def comparison_spectrum_dialog(self):
        """
        Open a dialog to pick comparison spectrum
        """
        path, _ = QtGui.QFileDialog.getOpenFileName(self,
            caption="Pick comparison spectrum", directory="", filter="")
        if not path: return
        spectrum = smh.specutils.Spectrum1D.read(path)
        self.stellar_parameters_tab.specfig.update_comparison_spectrum(spectrum)
        self.chemical_abundances_tab.figure.update_comparison_spectrum(spectrum)
        return None


    def clear_comparison_spectrum(self):
        """
        Remove the comparison spectrum
        """
        self.stellar_parameters_tab.specfig.update_comparison_spectrum(None)
        self.chemical_abundances_tab.figure.update_comparison_spectrum(None)
        return None


    def summary_plot(self):
        """
        Make a summary plot in a popup dialog
        """
        window = SummaryPlotDialog(self.session, self)
        window.exec_()
        return None

    def ncap_summary_plot(self):
        """
        Make a ncap summary plot in a popup dialog
        """
        window = SummaryPlotDialog(self.session, self, ncap=True)
        window.exec_()
        return None

    def snr_plot(self):
        """
        Make a snr plot (=1/ivar) in a popup dialog
        """
        window = SNRPlotDialog(self.session, self)
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
        self.chemical_abundances_tab \
            = chemical_abundances.ChemicalAbundancesTab(self)
        self.tabs.addTab(self.chemical_abundances_tab, "Line Measurements")

        # Create review tab
        self.review_tab \
            = review.ReviewTab(self)
        self.tabs.addTab(self.review_tab, "Review")

        # Disable all tabs except the first one.
        for i in range(self.tabs.count()):
            self.tabs.setTabEnabled(i, i == 0)

        cw_vbox.addWidget(self.tabs)
        self.setCentralWidget(cw)

        self.tabs.setCurrentIndex(0)
        self._update_window_title()

    def transition_dialog_callback(self):
        self.stellar_parameters_tab.new_session_loaded()
        self.chemical_abundances_tab.new_session_loaded()
        self.review_tab.new_session_loaded()

    def refresh_all_guis(self):
        raise NotImplementedError
