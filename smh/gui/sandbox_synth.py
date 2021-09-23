#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The main GUI window for Spectroscopy Made Hard. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
from PySide2 import (QtCore, QtWidgets as QtGui)
import yaml
import time, os, sys

import numpy as np

# Import functionality related to each tab
import rv, normalization, summary, stellar_parameters, chemical_abundances

# Functions related to warnings and exceptions.
import exception

import smh
from smh.linelists import LineList
import smh.spectral_models as sm

logger = logging.getLogger(__name__)

from ui_mainwindow import *

datadir = os.path.dirname(os.path.abspath(__file__))+'/../tests/test_data'
if __name__ == '__main__':

    import sys

    # Create the app and clean up any style bugs.
    try:
        app = QtGui.QApplication(sys.argv)

    except RuntimeError:
        # For development.
        None

    if sys.platform == "darwin":
            
        # See http://successfulsoftware.net/2013/10/23/fixing-qt-4-for-mac-os-x-10-9-mavericks/
        substitutes = [
            (".Lucida Grande UI", "Lucida Grande"),
            (".Helvetica Neue DeskInterface", "Helvetica Neue")
        ]
        for substitute in substitutes:
            QtGui.QFont.insertSubstitution(*substitute)

    # Create a global exception hook.
    sys._excepthook = sys.excepthook

    # Allow certain exceptions to be ignored, and these can be added to through
    # the GUI.
    ignore_exception_messages = []
    def exception_hook(exception_type, message, traceback):
        """
        An exception hook that will display a GUI and optionally allow the user
        to submit a GitHub issue.

        :param exception_type:
            The type of exception that was raised.

        :param message:
            The exception message.

        :param traceback:
            The traceback of the exception.
        """

        # Show the exception in the terminal.
        sys._excepthook(exception_type, message, traceback)

        # Should this exception be ignored?
        if message.__repr__() in ignore_exception_messages:
            return None

        # Load a GUI that shows the exception.
        exception_gui = exception.ExceptionWidget(
            exception_type, message, traceback)
        exception_gui.exec_()

        # Ignore future exceptions of this kind?
        if exception_gui.ignore_in_future:
            ignore_exception_messages.append(message.__repr__())

        return None

    sys.excepthook = exception_hook

    # Run the main application window.
    app.window = Ui_MainWindow(spectrum_filenames=[
        datadir+"/spectra/hd122563_1blue_multi_090205_oldbutgood.fits",
        datadir+"/spectra/hd122563_1red_multi_090205_oldbutgood.fits"
    ])
    # Enable all tabs
    for i in range(app.window.tabs.count()):
        app.window.tabs.setTabEnabled(i, True)
    session = app.window.session
    
    # TODO load a synthesis model here
    transitions = LineList.read(datadir+'/linelists/lin4554new')
    session.metadata["line_list"] = transitions
    synth = sm.SpectralSynthesisModel(session, transitions['hash'], 'Ba')
    session.metadata["spectral_models"] = [synth]

    with open(smh.Session._default_settings_path, "rb") as fp:
        defaults = yaml.load(fp)
        
    app.window.chemical_abundances_tab.new_session_loaded()
    app.window.tabs.setCurrentIndex(4)
    app.window.chemical_abundances_tab.fit_all()

    app.window.show()
    sys.exit(app.exec_())
