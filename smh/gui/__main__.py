#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The main GUI window for Spectroscopy Made Hard. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import sys

import logging
from PySide2 import (QtCore, QtGui as QtGui2, QtWidgets as QtGui)
import yaml
import os

# Fixing pyside2 bug: https://www.loekvandenouweland.com/content/pyside2-big-sur-does-not-show-window.html
os.environ['QT_MAC_WANTS_LAYER'] = '1'

# Hacky path solution to override the original path to smhr
here = os.path.abspath(os.path.dirname(__file__))
smhr_rpa_path = '/'.join(here.split('/')[:-2])
sys.path.insert(0,smhr_rpa_path)

# Functions related to warnings and exceptions.
import exception

import smh

logger = logging.getLogger(__name__)
logger.addHandler(smh.handler)

from ui_mainwindow import *

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
            QtGui2.QFont.insertSubstitution(*substitute)

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
    #app.window = Ui_MainWindow(spectrum_filenames=[
    #    "/Users/arc/Downloads/hd122563_1blue_multi_090205_oldbutgood.fits",
    #    "/Users/arc/Downloads/hd122563_1red_multi_090205_oldbutgood.fits"
    #])
    app.window = Ui_MainWindow(app)
    app.window.show()
    app.window.raise_()

    # DEBUG
    #testpath = os.path.dirname(os.path.abspath(__file__))+"/../tests/test_data/test_G64-12_v02.smh"
    #logger.debug("Loading {}".format(testpath))
    #app.window.open_session(path=testpath)
    
    # DEBUG
    #testpath = os.path.dirname(os.path.abspath(__file__))+"/../tests/test_data/spectra/hd122563.fits"
    #logger.debug("Loading {}".format(testpath))
    #app.window.new_session(testpath)
    #app.window.rv_tab.cross_correlate_and_correct()
    #app.window.session.metadata["normalization"]["continuum"] = [1]
    #app.window.session.metadata["normalization"]["normalization_kwargs"] = [{}]
    #app.window.normalization_tab.normalize_and_stitch()
    #app.window.tabs.setCurrentIndex(3)
    #app.window.session.import_master_list(os.path.dirname(os.path.abspath(__file__))+"/../tests/test_data/my_master_list.txt")
    #
    #session = app.window.session
    
    sys.exit(app.exec_())
