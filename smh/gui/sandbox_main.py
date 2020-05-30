#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The main GUI window for Spectroscopy Made Hard. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
logging.basicConfig()
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
    
    with open(smh.Session._default_settings_path, "rb") as fp:
        defaults = yaml.load(fp)
    datadir = os.path.dirname(os.path.abspath(__file__))+'/../tests/test_data'

    # Load in line_list
    ll = LineList.read(os.path.dirname(os.path.abspath(__file__))+'/../tests/test_data/linelists/complete.list')
    session.metadata['line_list'] = ll[190:200] #strong lines
    #session.metadata['line_list'] = ll[90:100]

    # Load in spectral_models from linelist
    #from linelist_manager import TransitionsDialog
    print("Loading spectral models..."); start = time.time()
    sm = []
    for hash in session.metadata["line_list"]["hash"]:
        sm.append(smh.spectral_models.ProfileFittingModel(session, [hash]))
    session.metadata['spectral_models'] = sm
    print("Done! {:.1f}s".format(time.time()-start))

    print("Fitting lines..."); start = time.time()
    app.window.chemical_abundances_tab.fit_all()
    print("Done! {:.1f}s".format(time.time()-start))
    app.window.tabs.setCurrentIndex(4)

    print("Measuring lines..."); start = time.time()
    app.window.chemical_abundances_tab.measure_all()
    _current_abundances = []
    _current_EW = []
    for m in session.metadata['spectral_models']:
        try:
            _current_abundances.append(m.metadata['fitted_result'][-1]['abundances'][0])
            _current_EW.append(m.metadata['fitted_result'][-1]['equivalent_width'][0])
        except KeyError:
            _current_abundances.append(np.nan)
            _current_EW.append(np.nan)
    print("Done! {:.1f}s".format(time.time()-start))
    
    #print("Measuring uncertainty from session..."); start = time.time()
    #abundances, uncertainties = session.measure_abundances()
    #print("Done! {:.1f}s".format(time.time()-start))
    
    #total_diff = np.nansum(np.abs(abundances-np.array(_current_abundances)))
    #assert total_diff == 0, total_diff

    app.window.show()
    sys.exit(app.exec_())
