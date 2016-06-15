#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Functionality to use matplotlib figures in PySide GUIs. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os
import matplotlib
from warnings import simplefilter
import time
import numpy as np

# Ignore warnings from matplotlib about fonts not being found.
simplefilter("ignore", UserWarning)

# Load our matplotlibrc file.
matplotlib.rc_file(os.path.join(os.path.dirname(__file__), "matplotlibrc"))

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.backends.backend_qt4agg \
#    import NavigationToolbar2QTAgg as NavigationToolbar

from matplotlib.figure import Figure

from PySide import QtCore, QtGui

DOUBLE_CLICK_INTERVAL = 0.1 # MAGIC HACK

class MPLWidget(FigureCanvas):
    """
    A widget to contain a matplotlib figure.
    """

    def __init__(self, parent=None, toolbar=False, tight_layout=True,
        autofocus=False, background_hack=True, **kwargs):
        """
        A widget to contain a matplotlib figure.

        :param autofocus: [optional]
            If set to `True`, the figure will be in focus when the mouse hovers
            over it so that keyboard shortcuts/matplotlib events can be used.

        Methods for zooming:
        axis_right_mouse_press: call for MPLWidget.mpl_connect("button_press_event")
        update_zoom_box: connect to MPLWidget.mpl_connect("motion_notify_event")
                         MPLWidget.mpl_disconnect as needed
        axis_right_mouse_release: call for MPLWidget.mpl_connect("button_release_event")
        unzoom_on_z_press: MPLWidget.mpl_connect("key_press_event")
        reset_zoom_limits(ax): call whenever you reset limits on an axis and want it to be zoomable
        """
        super(MPLWidget, self).__init__(Figure())
        
        self.figure = Figure(tight_layout=tight_layout)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setParent(parent)

        # Focus the canvas initially.
        self.canvas.setFocusPolicy(QtCore.Qt.WheelFocus)
        self.canvas.setFocus()

        self.toolbar = None #if not toolbar else NavigationToolbar(self, parent)

        # Get background of parent widget.

        # Because all these figures will be in a tab, we need to get the color
        # right. It seems impossible to get the *actual* color of the parent
        # background when the widget is in a tab, but it seems it is just 10
        # points darker.
        #if background_hack is not None:
        #    bg_color = [(_ - 10)/255. for _ in \
        #        parent.palette().color(QtGui.QPalette.Window).getRgb()[:3]]
        #    background_hack.figure.patch.set_facecolor(bg_color)

        if background_hack:
            self.figure.patch.set_facecolor([(_ - 10)/255. for _ in \
                self.palette().color(QtGui.QPalette.Window).getRgb()[:3]])


        if autofocus:
            self._autofocus_cid = self.canvas.mpl_connect(
                "figure_enter_event", self._focus)

        # State for zoom box
        self.old_xlim = {}
        self.old_ylim = {}
        self.right_clicked_axis = None
        self.zoom_box_lines = {}
        
        return None

    def _focus(self, event):
        """ Set the focus of the canvas. """
        self.canvas.setFocus()

    ######################
    # Methods for zooming
    ######################

    def get_current_axis(self,event):
        for i,ax in enumerate(self.figure.axes):
            if event.inaxes in [ax]:
                return i,ax
        raise RuntimeError("Cannot find correct axis")
    
    def axis_right_mouse_press(self, event):
        """
        Right mouse button pressed in axis
        """
        
        if event.button != 3: return None

        self.x1 = event.xdata
        self.y1 = event.ydata
        try:
            i,ax = self.get_current_axis(event)
        except RuntimeError:
            return None

        # Save original limits if not previously saved
        if i not in self.old_xlim or \
        self.old_xlim[i] is None:
            self.old_xlim[i] = ax.get_xlim()
            self.old_ylim[i] = ax.get_ylim()

        # Mark this as currently rightclicked axis
        self.right_clicked_axis = ax

        # Create lines if needed
        if i not in self.zoom_box_lines:
            self.zoom_box_lines[i] = [ax.plot([np.nan,np.nan],[np.nan,np.nan],'k--')[0],
                                      ax.plot([np.nan,np.nan],[np.nan,np.nan],'k--')[0],
                                      ax.plot([np.nan,np.nan],[np.nan,np.nan],'k--')[0],
                                      ax.plot([np.nan,np.nan],[np.nan,np.nan],'k--')[0]]
        
        # Create zoom box signal
        self._interactive_zoom_box_signal = (
            time.time(),
            self.mpl_connect("motion_notify_event",
                             self.update_zoom_box)
        )
        return None
    
    def update_zoom_box(self, event):
        """
        Updated the zoom box
        """
        self.x2 = event.xdata
        self.y2 = event.ydata
        try:
            i,ax = self.get_current_axis(event)
        except RuntimeError:
            return None

        # Don't update if not in current axis
        if ax != self.right_clicked_axis:
            return None

        # Skip if doubleclick
        signal_time, signal_cid = self._interactive_zoom_box_signal
        if time.time() - signal_time < DOUBLE_CLICK_INTERVAL:
            return None

        #  Plot dashed box
        self.zoom_box_lines[i][0].set_xdata([self.x1,self.x2])
        self.zoom_box_lines[i][0].set_ydata([self.y1,self.y1])
        self.zoom_box_lines[i][1].set_xdata([self.x1,self.x2])
        self.zoom_box_lines[i][1].set_ydata([self.y2,self.y2])
        self.zoom_box_lines[i][2].set_xdata([self.x1,self.x1])
        self.zoom_box_lines[i][2].set_ydata([self.y1,self.y2])
        self.zoom_box_lines[i][3].set_xdata([self.x2,self.x2])
        self.zoom_box_lines[i][3].set_ydata([self.y1,self.y2])
        self.draw()
        return None

    def axis_right_mouse_release(self, event):
        """
        Right mouse button released in axis
        """

        if event.button != 3: return None

        # Skip if not mouse pressed
        try:
            signal_time, signal_cid = self._interactive_zoom_box_signal
        except AttributeError:
            return None

        # Skip if doubleclick
        if time.time() - signal_time < DOUBLE_CLICK_INTERVAL:
            return None

        self.x2 = event.xdata
        self.y2 = event.ydata
        try:
            i,ax = self.get_current_axis(event)
        except RuntimeError:
            # Pass because still have to clean up
            pass
        else:
            # Apply the zoom
            if ax == self.right_clicked_axis:
                x1 = min(self.x1, self.x2)
                x2 = max(self.x1, self.x2)
                y1 = min(self.y1, self.y2)
                y2 = max(self.y1, self.y2)
                ax.set_xlim([x1,x2])
                ax.set_ylim([y1,y2])
        
        # Clean up
        for i,_ax in enumerate(self.figure.axes):
            if _ax == self.right_clicked_axis:
                break # set i
        else:
            raise RuntimeError("Cannot find right clicked axis")
        self.zoom_box_lines[i][0].set_xdata([np.nan,np.nan])
        self.zoom_box_lines[i][0].set_xdata([np.nan,np.nan])
        self.zoom_box_lines[i][1].set_xdata([np.nan,np.nan])
        self.zoom_box_lines[i][1].set_xdata([np.nan,np.nan])
        self.zoom_box_lines[i][2].set_xdata([np.nan,np.nan])
        self.zoom_box_lines[i][2].set_xdata([np.nan,np.nan])
        self.zoom_box_lines[i][3].set_xdata([np.nan,np.nan])
        self.zoom_box_lines[i][3].set_xdata([np.nan,np.nan])
        self.right_clicked_axis = None

        self.mpl_disconnect(signal_cid)
        self.draw()
        del self._interactive_zoom_box_signal

        return None
        
    def unzoom_on_z_press(self, event):
        """
        Unzoom event when keyboard "z" is pressed
        """
        
        if event.key not in ["z","Z"]: return None
        # right_clicked_axis is not None if you are 
        # currently holding down the right mouse button
        if self.right_clicked_axis is not None: return None
        
        #print("Key press",event,event.key)
        
        try:
            i,ax = self.get_current_axis(event)
        except RuntimeError:
            return None
        
        # Skip if no zoom applied yet
        if i not in self.old_xlim: return None
        
        # Get original limits
        old_xlim = self.old_xlim[i]
        old_ylim = self.old_ylim[i]
        ax.set_xlim(old_xlim)
        ax.set_ylim(old_ylim)
        self.draw()
        
        # Reset zoom
        self.old_xlim[i] = None
        self.old_ylim[i] = None

    def reset_zoom_limits(self, ax=None):
        """
        Need to call this after setting x/y limits not through MPLWidget
        Otherwise MPLWidget.unzoom will reset to old limits
        """
        for i,_ax in enumerate(self.figure.axes):
            if ax is None:
                self.old_xlim[i] = None
                self.old_ylim[i] = None
            elif ax==_ax:
                self.old_xlim[i] = None
                self.old_ylim[i] = None
                return None
        if ax is None: return None
        raise ValueError("Could not identify axis to reset zoom limits for")
