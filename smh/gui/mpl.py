#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Functionality to use matplotlib figures in PySide2 GUIs. """

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

from PySide2 import (QtCore, QtGui as QtGui2, QtWidgets as QtGui)

DOUBLE_CLICK_INTERVAL = 0.1 # MAGIC HACK

class MPLWidget(FigureCanvas):
    """
    A widget to contain a matplotlib figure.
    """

    __double_click_interval = 0.1
    __right_double_click_interval = 0.2

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
                self.palette().color(QtGui2.QPalette.Window).getRgb()[:3]])


        if autofocus:
            self._autofocus_cid = self.canvas.mpl_connect(
                "figure_enter_event", self._focus)

        # State for zoom box
        self._right_click_zoom_box = {}

        # State for shift. (Antimasks)
        self.shift_key_pressed = False
        
        # State for space. (Pan)
        self.space_key_pressed = False 
       

    def _focus(self, event):
        """ Set the focus of the canvas. """
        print("_focus", event)
        self.canvas.setFocus()

    ######################
    # Methods for zooming
    ######################

    def get_current_axis(self,event):
        for i,ax in enumerate(self.figure.axes):
            if event.inaxes in [ax]:
                return i,ax
        raise RuntimeError("Cannot find correct axis")
    

    def enable_interactive_zoom(self):

        self._interactive_zoom_history = {}
        self._interactive_zoom_signals = [
            self.mpl_connect("button_press_event", self._interactive_zoom_press),
            self.mpl_connect("button_release_event", self._interactive_zoom_release),
        ]


    def disable_interactive_zoom(self):
        try:
            for cid in self._interactive_zoom_signals:
                self.mpl_disconnect(cid)

        except AttributeError:
            return None

        else:
            del self._interactive_zoom_signals

        return None



    def _interactive_zoom_press(self, event):
        """
        Right-mouse button pressed in axis.

        :param event:
            A matplotlib event.
        """
        
        if event.button != 3 or event.inaxes is None:
            return None

        self._interactive_zoom_initial_bounds = [event.xdata, event.ydata]
        self._interactive_zoom_axis_index = self.figure.axes.index(event.inaxes)

        # Create lines if needed
        if self._interactive_zoom_axis_index not in self._right_click_zoom_box:
            self._right_click_zoom_box[self._interactive_zoom_axis_index] \
                = event.inaxes.plot([np.nan], [np.nan], "k:")[0]

        # Create zoom box signal
        self._interactive_zoom_box_signal = (
            time.time(),
            self.mpl_connect("motion_notify_event", self._update_interactive_zoom)
        )
        return None
    

    def _update_interactive_zoom(self, event):
        """
        Updated the zoom box.
        """

        if event.inaxes is None \
        or self.figure.axes[self._interactive_zoom_axis_index] != event.inaxes:
            return None

        # Skip if doubleclick
        # TODO: return to this.
        signal_time, signal_cid = self._interactive_zoom_box_signal
        if time.time() - signal_time < self.__right_double_click_interval:
            return None

        xs, ys = self._interactive_zoom_initial_bounds
        xe, ye = event.xdata, event.ydata

        self._right_click_zoom_box[self._interactive_zoom_axis_index].set_data(
            np.array([
                [xs, xe, xe, xs, xs],
                [ys, ys, ye, ye, ys]
            ])
        )
        self.draw()
        return None


    def _interactive_zoom_release(self, event):
        """
        Right mouse button released in axis.
        """

        if event.button != 3 or event.inaxes is None \
        or self.figure.axes[self._interactive_zoom_axis_index] != event.inaxes:

            # The right-click mouse button was released outside of the axis that
            # we want. Disconnect the signal and let them try again.
            try:
                signal_time, signal_cid = self._interactive_zoom_box_signal

            except AttributeError:
                return None

            else:
                self._right_click_zoom_box[self._interactive_zoom_axis_index].set_data(
                    np.array([[np.nan], [np.nan]]))

                self.mpl_disconnect(signal_cid)
                self.draw()
                del self._interactive_zoom_box_signal

            return None

        self._interactive_zoom_history.setdefault(
            self._interactive_zoom_axis_index, 
            [(event.inaxes.get_xlim(), event.inaxes.get_ylim())])

        # If this is a double-click, go to the previous zoom in history. 
        signal_time, signal_cid = self._interactive_zoom_box_signal
        if time.time() - signal_time < self.__right_double_click_interval:
            
            if len(self._interactive_zoom_history[self._interactive_zoom_axis_index]) > 1:
                xlim, ylim = self._interactive_zoom_history[self._interactive_zoom_axis_index].pop(-1)
            else:
                xlim, ylim = self._interactive_zoom_history[self._interactive_zoom_axis_index][-1]

        else:        
            # Set the zoom.
            xlim = np.sort([self._interactive_zoom_initial_bounds[0], event.xdata])
            ylim = np.sort([self._interactive_zoom_initial_bounds[1], event.ydata])

            self._interactive_zoom_history[self._interactive_zoom_axis_index].append(
                [xlim, ylim])

        # Hide the zoom box.
        self._right_click_zoom_box[self._interactive_zoom_axis_index].set_data(
            np.array([[np.nan], [np.nan]]))

        event.inaxes.set_xlim(xlim)
        event.inaxes.set_ylim(ylim)

        self.mpl_disconnect(signal_cid)
        self.draw()

        del self._interactive_zoom_box_signal

        return None


    def reset_zoom_limits(self):
        # Get the first entry from the zoom history.
        for index, limits in self._interactive_zoom_history.items():
            xlim, ylim = limits.pop(0)
            self.figure.axes[index].set_xlim(xlim)
            self.figure.axes[index].set_ylim(ylim)

        self._interactive_zoom_history = {}
        self.draw()
        return None


    def _clear_zoom_history(self, axis_index=None):

        if axis_index is None:
            self._interactive_zoom_history = {}

        elif axis_index in self._interactive_zoom_history:
            del self._interactive_zoom_history[axis_index]
            
        return None


    def enable_drag_to_mask(self, axes):
        """
        Enable drag-to-mask on axes in this figure.

        :param axes:
            The axis or list of axes in this figure where drag-to-mask should
            be enabled. If a list of axes are given, then masking in one will
            display in all others.
        """

        if isinstance(axes, matplotlib.axes.Axes):
            axes = [axes]

        for i, ax in enumerate(axes):
            if ax not in self.figure.axes:
                raise ValueError(
                    "axes (zero-index {}) is not in this figure".format(i))

        self._drag_to_mask_axes = axes

        # Connect events.
        self.mpl_connect("button_press_event", self._drag_to_mask_press)
        self.mpl_connect("button_release_event", self._drag_to_mask_release)

        self.dragged_masks = []

        self._mask_interactive_region = dict(zip(axes, [None] * len(axes)))

        self._masked_regions = {}
        [self._masked_regions.setdefault(ax, []) for ax in axes]

        return None



    def _drag_to_mask_press(self, event):
        """
        Event for triggering when the left mouse button has been clicked, just
        before a region is dragged to be masked.

        :param event:
            The matplotlib event.
        """

        if event.inaxes not in self._drag_to_mask_axes or event.button != 1:
            return None

        if not event.dblclick:

            # Single left-hand mouse button click.
            xmin, xmax, ymin, ymax = (event.xdata, np.nan, -1e8, +1e8)

            for ax in self._mask_interactive_region.keys():
                interactive_region = self._mask_interactive_region[ax]
                if interactive_region is None:
                    self._mask_interactive_region[ax] = ax.axvspan(**{
                            "xmin": xmin,
                            "xmax": xmax,
                            "ymin": ymin,
                            "ymax": ymax,
                            "facecolor": "r",
                            "edgecolor": "none",
                            "alpha": 0.25,
                            "zorder": -1
                        })

                else:
                    interactive_region.set_xy([
                        [xmin, ymin],
                        [xmin, ymax],
                        [xmax, ymax],
                        [xmax, ymin],
                        [xmin, ymin]
                    ])

            # Set the signal and the time.
            self._mask_interactive_region_signal = (
                time.time(),
                self.mpl_connect(
                    "motion_notify_event", self._drag_to_mask_motion)
            )
        
        else:
            # Double click.

            # Matches any existing masked regions?
            delete_mask_index = None
            for i, (xmin, xmax) in enumerate(self.dragged_masks):
                if xmax >= event.xdata >= xmin:
                    # Remove this mask.
                    delete_mask_index = i
                    break

            if delete_mask_index is not None:
                self.dragged_masks.pop(delete_mask_index)
                self._draw_dragged_masks()

            return None

        return None


    def _drag_to_mask_motion(self, event):
        """
        Event for triggering when the mouse has been moved while the left-hand
        button is held, indicating a change in the region to be masked.

        :param event:
            The matplotlib event.
        """

        if event.xdata is None or event.inaxes not in self._drag_to_mask_axes \
        or event.button != 1:
            return None 

        signal_time, signal_cid = self._mask_interactive_region_signal
        if time.time() - signal_time > self.__double_click_interval:

            # Update all interactive masks.
            for ax in self._mask_interactive_region.keys():
                data = self._mask_interactive_region[ax].get_xy()

                # Update xmax.
                data[2:4, 0] = event.xdata
                self._mask_interactive_region[ax].set_xy(data)

            self.draw()

        return None


    def _drag_to_mask_release(self, event):
        """
        Event for triggering when the left-hand mouse button has been released
        and a region has been marked to be masked.

        :param event:
            The matplotlib event.
        """

        try:
            signal_time, signal_cid = self._mask_interactive_region_signal

        except AttributeError:
            return None
        
        xy = self._mask_interactive_region.values()[0].get_xy()
        
        xdata = event.xdata if event.xdata is not None else xy[2, 0]

        # If the two mouse events were within some time interval,
        # then we should not add a mask because those signals were probably
        # part of a double-click event.
        if  time.time() - signal_time > self.__double_click_interval \
        and np.abs(xy[0,0] - xdata) > 0:
            
            # Update the cache with the new mask.
            limits = [xy[0, 0], xy[2, 0]]
            self.dragged_masks.append([min(limits), max(limits)])

            # Update the view.
            self._draw_dragged_masks()

        # Clear the interactive masks.
        xy[:, 0] = np.nan
        for ax in self._mask_interactive_region.keys():
            self._mask_interactive_region[ax].set_xy(xy)

        self.mpl_disconnect(signal_cid)
        del self._mask_interactive_region_signal

        self.draw()

        return None


    def _draw_dragged_masks(self):
        """ Draw the dragged masks in the relevant axes. """

        ymin, ymax = (-1e8, +1e8)
        for i, (xmin, xmax) in enumerate(self.dragged_masks):

            for ax in self._masked_regions.keys():
                try:
                    masked_region = self._masked_regions[ax][i]

                except IndexError:
                    self._masked_regions[ax].append(ax.axvspan(**{
                        "xmin": xmin,
                        "xmax": xmax,
                        "ymin": ymin,
                        "ymax": ymax,
                        "facecolor": "r",
                        "edgecolor": "none",
                        "alpha": 0.25,
                        "zorder": -1
                    }))

                else:
                    masked_region.set_xy([
                        [xmin, ymin],
                        [xmin, ymax],
                        [xmax, ymax],
                        [xmax, ymin],
                        [xmin, ymin]
                    ])

        # Clear off any additional ones.
        N = len(self.dragged_masks)
        for ax in self._masked_regions.keys():
            for masked_region in self._masked_regions[ax][N:]:
                masked_region.set_xy([
                    [np.nan, np.nan],
                    [np.nan, np.nan],
                    [np.nan, np.nan],
                    [np.nan, np.nan],
                    [np.nan, np.nan],
                ])

        self.draw()

        return None


    def key_press_flags(self, event):
        if event.key == "shift":
            self.shift_key_pressed = True
        elif event.key == "space":
            self.space_key_pressed = True
        return None
        
    def key_release_flags(self, event):
        if event.key == "shift":
            self.shift_key_pressed = False
        elif event.key == "space":
            self.space_key_pressed = False
        return None        
