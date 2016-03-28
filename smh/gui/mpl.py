#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Functionality to use matplotlib figures in PySide GUIs. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os
import matplotlib

# Load our matplotlibrc file.
matplotlib.rc_file(os.path.join(os.path.dirname(__file__), "matplotlibrc"))

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.backends.backend_qt4agg \
#    import NavigationToolbar2QTAgg as NavigationToolbar

from matplotlib.figure import Figure

from PySide import QtGui


class MPLWidget(FigureCanvas):
    """
    A widget to contain a matplotlib figure.
    """

    def __init__(self, parent=None, toolbar=False, tight_layout=True):
        super(MPLWidget, self).__init__(Figure())

        self.setParent(parent)

        self.figure = Figure(tight_layout=tight_layout)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = None #if not toolbar else NavigationToolbar(self, parent)

        # Get background of parent widget.

        # Because all these figures will be in a tab, we need to get the color
        # right. It seems impossible to get the *actual* color of the parent
        # background when the widget is in a tab, but it seems it is just 10
        # points darker.
        if parent is not None:
            bg_color = [(_ - 10)/255. for _ in \
                parent.palette().color(QtGui.QPalette.Window).getRgb()[:3]]
            self.figure.patch.set_facecolor(bg_color)

        return None
