#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The radial velocity tab view for the Spectroscopy Made Hard GUI. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import os
import sys
from PySide2 import (QtCore, QtWidgets as QtGui)

import mpl
from matplotlib import (gridspec, pyplot as plt)
# E. Holmbeck added this line to account for a newer version of matplotlib
plt.rcParams['errorbar.capsize'] = 3
from smh import (Session, isoutils)
from smh.linelists import LineList

import logging
logger = logging.getLogger(__name__)

__all__ = ["SummaryPlotDialog"]

class SummaryPlotDialog(QtGui.QDialog):
    def __init__(self, session, parent, ncap=False):
        super(SummaryPlotDialog, self).__init__()
        # Center it on the parent location
        rect = parent.geometry()
        x, y = rect.center().x(), rect.center().y()
        w = 1000
        h = 640
        self.setGeometry(x-w/2, y-h/2, w, h)
        vbox = QtGui.QVBoxLayout(self)
        figure_widget = mpl.MPLWidget(None, tight_layout=True)
        vbox.addWidget(figure_widget)
        if ncap:
            session.make_ncap_summary_plot(figure_widget.figure)
        else:
            session.make_summary_plot(figure_widget.figure)
        
class SNRPlotDialog(QtGui.QDialog):
    def __init__(self, session, parent, ncap=False):
        super(SNRPlotDialog, self).__init__()
        # Center it on the parent location
        rect = parent.geometry()
        x, y = rect.center().x(), rect.center().y()
        w = 1000
        h = 640
        self.setGeometry(x-w/2, y-h/2, w, h)
        vbox = QtGui.QVBoxLayout(self)
        figure_widget = mpl.MPLWidget(None, tight_layout=True)
        figure_widget.enable_interactive_zoom()
        vbox.addWidget(figure_widget)
        session.make_snr_plot(figure_widget.figure)

