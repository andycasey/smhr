#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Functionality to use matplotlib figures in PySide GUIs. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import matplotlib

matplotlib.use("Qt4Agg")
# HACK TODO: Use a .matplotlibrc file
#matplotlib.rcParams["backend.qt4"] = "PySide"
matplotlib.rcParams["lines.antialiased"] = True
matplotlib.rcParams["text.antialiased"] = True

import matplotlib.pyplot as plt

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from PySide import QtGui

class MPLWidget(FigureCanvas):

    def __init__(self, parent=None ,xlabel='x',ylabel='y',title='Title'):
        super(MPLWidget, self).__init__(Figure())

        self.setParent(parent)

        self.figure = Figure(dpi=72, tight_layout=True)
        self.canvas = FigureCanvas(self.figure)

        # Get background of parent widget.

        # Because all these figures will be in a tab, we need to get the color
        # right. It seems impossible to get the *actual* color of the parent
        # background when the widget is in a tab, but it seems it is just 10
        # points darker.
        if parent is not None:
            bg_color = [(_ - 10)/255. for _ in \
                parent.palette().color(QtGui.QPalette.Window).getRgb()[:3]]
            self.figure.patch.set_facecolor(bg_color)

        self.axes = self.figure.add_subplot(111, axisbg="#FFFFFF")
        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)
        self.axes.set_title(title)
