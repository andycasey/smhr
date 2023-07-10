#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The normalization tab in Spectroscopy Made Hard """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["NormalizationTab"]

import logging
import numpy as np
import sys
from PySide2 import (QtCore, QtGui as QtGui2, QtWidgets as QtGui)
from time import time

from matplotlib import gridspec

import mpl

# This is a bad idea, but it's April 1st.
from style_utils import wavelength_to_hex


logger = logging.getLogger(__name__)

c = 299792458e-3 # km/s

# The minimum time (in seconds) between a mouse click/release to differentiate
# a single click from a double click
#DOUBLE_CLICK_INTERVAL = 0.1 # MAGIC HACK
# E. Holmbeck changed this since we don't have double-click anymore
DOUBLE_CLICK_INTERVAL = 0.0 # MAGIC HACK

# The pixel tolerance to select and remove an additional point.
PIXEL_PICKER_TOLERANCE = 30 # MAGIC HACK


def dict_updated(default, new, exclude=None):
    updated = {}
    exclude = exclude or ()
    for key, value in default.items():
        if key in new and key not in exclude and new[key] != value:
            updated[key] = (value, new[key])
    return updated


class NormalizationTab(QtGui.QWidget):

    def __init__(self, parent):
        """
        Create a tab for the normalization of spectra.
        
        :param parent:
            The parent widget.
        """

        super(NormalizationTab, self).__init__(parent)
        self.parent = parent

        # Establish the GUI for this tab.
        sp = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, 
            QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)

        # Create a top-level horizontal layout to contain a MPL figure and
        # a vertical layout of settings..
        tab_layout = QtGui.QHBoxLayout(self)
        tab_layout.setContentsMargins(20, 20, 20, 0)
        
        settings_widget = QtGui.QWidget()
        settings_layout = QtGui.QVBoxLayout(settings_widget)
        settings_widget.setFixedWidth(300)
        
        # Start the grid layout for the normalization tab.
        settings_grid_layout = QtGui.QGridLayout()

        # Normalization function.
        self.function_label = QtGui.QLabel(self)
        self.function_label.setText("Function")
        settings_grid_layout.addWidget(self.function_label, 0, 0, 1, 1)
        
        # Put the normalization function combo box in a horizontal layout with 
        # a spacer.
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.function = QtGui.QComboBox(self)
        self.function.setObjectName("norm_function")
        hbox.addWidget(self.function)
        settings_grid_layout.addLayout(hbox, 0, 1, 1, 1)

        for each in ("polynomial", "spline"):
            self.function.addItem(each.title())

        # Normalization function order.
        self.order_label = QtGui.QLabel(self)
        self.order_label.setText("Order")
        settings_grid_layout.addWidget(self.order_label, 1, 0, 1, 1)
        
        # Put the normalization order combo box in a horizontal layout with a
        # spacer
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.order = QtGui.QComboBox(self)
        self.order.setMaximumSize(QtCore.QSize(50, 16777215))
        self.order.setObjectName("norm_order")
        hbox.addWidget(self.order)
        settings_grid_layout.addLayout(hbox, 1, 1, 1, 1)

        orders = range(1, 10)
        for order in orders:
            self.order.addItem("{0:.0f}".format(order))

        # Maximum number of iterations.
        self.max_iter_label = QtGui.QLabel(self)
        self.max_iter_label.setText("Maximum iterations")
        settings_grid_layout.addWidget(self.max_iter_label, 2, 0, 1, 1)

        # Put the maxium number of iterations in a horizontal layout with a 
        # spacer.
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_max_iter = QtGui.QComboBox(self)
        self.norm_max_iter.setMaximumSize(QtCore.QSize(50, 16777215))
        self.norm_max_iter.setObjectName("norm_norm_max_iter")
        hbox.addWidget(self.norm_max_iter)
        settings_grid_layout.addLayout(hbox, 2, 1, 1, 1)

        norm_max_iters = range(1, 10)
        for iteration in norm_max_iters:
            self.norm_max_iter.addItem("{0:.0f}".format(iteration))


        # Low sigma clipping.
        self.low_sigma_clip_label = QtGui.QLabel(self)
        self.low_sigma_clip_label.setText("Low sigma clip")
        settings_grid_layout.addWidget(self.low_sigma_clip_label, 3, 0, 1, 1)

        # Put the low sigma line edit box in a horizontal layout with a spacer.
        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(-1, -1, 5, -1)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.low_sigma_clip = QtGui.QLineEdit(self)
        self.low_sigma_clip.setMaximumSize(QtCore.QSize(40, 16777215))
        self.low_sigma_clip.setAlignment(QtCore.Qt.AlignCenter)
        self.low_sigma_clip.setObjectName("norm_low_sigma_clip")
        self.low_sigma_clip.setValidator(
            QtGui2.QDoubleValidator(0, 1000, 2, self.low_sigma_clip))

        hbox.addWidget(self.low_sigma_clip)
        settings_grid_layout.addLayout(hbox, 3, 1, 1, 1)


        # High sigma clipping.
        self.high_sigma_clip_label = QtGui.QLabel(self)
        self.high_sigma_clip_label.setText("High sigma clip")
        settings_grid_layout.addWidget(self.high_sigma_clip_label, 4, 0, 1, 1)

        # Put the high sigma line edit box in a horizontal layout with a spacer.
        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(-1, -1, 5, -1)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.high_sigma_clip = QtGui.QLineEdit(self)
        self.high_sigma_clip.setMaximumSize(QtCore.QSize(40, 16777215))
        self.high_sigma_clip.setAlignment(QtCore.Qt.AlignCenter)
        self.high_sigma_clip.setObjectName("norm_high_sigma_clip")
        self.high_sigma_clip.setValidator(
            QtGui2.QDoubleValidator(0, 1000, 2, self.high_sigma_clip))
        hbox.addWidget(self.high_sigma_clip)
        settings_grid_layout.addLayout(hbox, 4, 1, 1, 1)
        

        # Knot spacing.
        self.knot_spacing_label = QtGui.QLabel(self)
        settings_grid_layout.addWidget(self.knot_spacing_label, 5, 0, 1, 1)
        self.knot_spacing_label.setText(u"Knot spacing (Å)")

        # Put the knot spacing lint edit box in a horizontal layout with a spacer
        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(-1, -1, 5, -1)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.knot_spacing = QtGui.QLineEdit(self)
        self.knot_spacing.setMaximumSize(QtCore.QSize(40, 16777215))
        self.knot_spacing.setAlignment(QtCore.Qt.AlignCenter)
        self.knot_spacing.setValidator(
            QtGui2.QDoubleValidator(0, 10000, 0, self.knot_spacing))
        self.knot_spacing.setObjectName("norm_knot_spacing")
        hbox.addWidget(self.knot_spacing)
        settings_grid_layout.addLayout(hbox, 5, 1, 1, 1)

        # -----------------------------------------------------------------
        # E. Holmbeck added these lines
        # Blue trimming.
        self.blue_trim_label = QtGui.QLabel(self)
        settings_grid_layout.addWidget(self.blue_trim_label, 6, 0, 1, 1)
        self.blue_trim_label.setText(u"Blue trim (pixels)")

        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(-1, -1, 5, -1)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.blue_trim = QtGui.QLineEdit(self)
        self.blue_trim.setMaximumSize(QtCore.QSize(40, 16777215))
        self.blue_trim.setAlignment(QtCore.Qt.AlignCenter)
        self.blue_trim.setValidator(
            QtGui2.QDoubleValidator(0, 10000, 0, self.blue_trim))
        self.blue_trim.setObjectName("blue_trim")
        hbox.addWidget(self.blue_trim)
        settings_grid_layout.addLayout(hbox, 6, 1, 1, 1)

        # Red trimming.
        self.red_trim_label = QtGui.QLabel(self)
        settings_grid_layout.addWidget(self.red_trim_label, 7, 0, 1, 1)
        self.red_trim_label.setText(u"Red trim (pixels)")

        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(-1, -1, 5, -1)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.red_trim = QtGui.QLineEdit(self)
        self.red_trim.setMaximumSize(QtCore.QSize(40, 16777215))
        self.red_trim.setAlignment(QtCore.Qt.AlignCenter)
        self.red_trim.setValidator(
            QtGui2.QDoubleValidator(0, 10000, 0, self.red_trim))
        self.red_trim.setObjectName("red_trim")
        hbox.addWidget(self.red_trim)
        settings_grid_layout.addLayout(hbox, 7, 1, 1, 1)
        # -----------------------------------------------------------------
        
        # End of the grid in the normalization tab.
        settings_layout.addLayout(settings_grid_layout)

        # Add a label.
        label = QtGui.QLabel(self)
        label.setText("Global continuum mask:")
        settings_layout.addWidget(label)

        # Add options for continuum mask.
        hbox = QtGui.QHBoxLayout()
        self.continuum_mask = QtGui.QComboBox(self)
        self.continuum_mask.setObjectName("contiuum_mask")
        hbox.addWidget(self.continuum_mask)
        hbox.addItem(
            QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
                QtGui.QSizePolicy.Minimum))
        self.edit_continuum_masks = QtGui.QPushButton(self)
        self.edit_continuum_masks.setObjectName("edit_continuum_masks")
        self.edit_continuum_masks.setText("Edit masks..")
        hbox.addWidget(self.edit_continuum_masks)

        settings_layout.addLayout(hbox)

        # Add a 'normalize and stitch button'
        self.stitch_btn = QtGui.QPushButton(self)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(self.stitch_btn.sizePolicy().hasHeightForWidth())
        self.stitch_btn.setSizePolicy(sp)
        self.stitch_btn.setMinimumSize(QtCore.QSize(250, 0))
        font = QtGui2.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.stitch_btn.setFont(font)
        self.stitch_btn.setCursor(QtGui2.QCursor(QtCore.Qt.PointingHandCursor))
        self.stitch_btn.setDefault(True)
        self.stitch_btn.setObjectName("stitch_btn")
        self.stitch_btn.setText("Normalize and stitch orders")

        settings_layout.addWidget(self.stitch_btn)

        # Add a spacer.
        settings_layout.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding))

        tab_layout.addWidget(settings_widget)

        # Create a matplotlib widget.
        self.norm_plot = mpl.MPLWidget(None, tight_layout=True, matchbg=self)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(self.norm_plot.sizePolicy().hasHeightForWidth())
        self.norm_plot.setSizePolicy(sp)
        self.norm_plot.setFocusPolicy(QtCore.Qt.StrongFocus)

        self.order_slide = QtGui.QSlider(self)
        self.order_slide.setGeometry(QtCore.QRect(230, 200, 160, 22))
        self.order_slide.setOrientation(QtCore.Qt.Horizontal)
        self.order_slide.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.order_slide.setMaximum(15)
        self.order_slide.setOrientation(QtCore.Qt.Horizontal)
        self.order_slide.setTickInterval(1)
        self.current_order_label = QtGui.QLabel(self)
        self.current_order_label.setText("")
        self.order_slide.valueChanged.connect(self.update_order_figure)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.norm_plot)
        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(self.order_slide)
        hbox.addWidget(self.current_order_label)
        vbox.addLayout(hbox)

        tab_layout.addLayout(vbox)

        # Set up the plot.
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        self.ax_order = self.norm_plot.figure.add_subplot(gs[0])
        # Line for the data.
        self.ax_order.plot([np.nan], [np.nan], c='k', zorder=3)#, drawstyle='steps-mid')
        # Line for the continuum.
        self.ax_order.plot([], [], linestyle="--", linewidth=2, c='r', zorder=4)
        # Points for the continuum knots.
        self.ax_order.plot([], [], 'o', mfc='none', mec='r', zorder=5, mew=1, ms=10)

        # Line for the neighbouring order(s) (joined by a NaN).
        self.ax_order.plot([np.nan], [np.nan], c='#666666', zorder=1, drawstyle='steps-mid')
        # Line for the neighbouring order(s) continuum (joined by a NaN)
        self.ax_order.plot([np.nan], [np.nan], c='b', zorder=2)

        # Additional point markers.
        self.ax_order.scatter([], [], facecolor="k", zorder=5, picker=5)

        self.ax_order.set_xticklabels([])
        self.ax_order.set_ylabel("Flux")

        self.ax_order_norm = self.norm_plot.figure.add_subplot(gs[1])
        self.ax_order_norm.axhline(1, linestyle=":", c="#666666", zorder=1)
        self.ax_order_norm.plot([np.nan], [np.nan], c='k', zorder=2)

        # TODO: Make (0, 1.2) a default view setting.
        self.ax_order_norm.set_ylim(0, 1.2)
        self.ax_order_norm.set_yticks([0, 0.5, 1.0])
        self.ax_order_norm.set_xlabel(u"Wavelength (Å)")

        self.norm_plot.draw()

        # Create signals.
        self.stitch_btn.clicked.connect(self.normalize_and_stitch)

        self.norm_plot.mpl_connect(
            "key_press_event", self.figure_key_press)
        self.norm_plot.mpl_connect(
            "button_press_event", self.figure_mouse_press)
        self.norm_plot.mpl_connect(
            "button_release_event", self.figure_mouse_release)

        # Zoom box
        #self.norm_plot.mpl_connect(
        #    "button_press_event", self.norm_plot.axis_right_mouse_press)
        #self.norm_plot.mpl_connect(
        #    "button_release_event", self.norm_plot.axis_right_mouse_release)
        #self.norm_plot.mpl_connect(
        #    "key_press_event", self.norm_plot.unzoom_on_z_press)
        self.norm_plot.enable_interactive_zoom()
        
        self.function.currentIndexChanged.connect(
            self.update_normalization_function)
        self.order.currentIndexChanged.connect(
            self.update_normalization_order)
        self.norm_max_iter.currentIndexChanged.connect(
            self.update_normalization_max_iterations)
        self.low_sigma_clip.textChanged.connect(
            self.update_low_sigma_clip)
        self.high_sigma_clip.textChanged.connect(
            self.update_high_sigma_clip)
        self.knot_spacing.textChanged.connect(self.update_knot_spacing)
        # -----------------------------------------------------------------
        # E. Holmbeck added connect to trim
        self.blue_trim.textChanged.connect(self.update_blue_trim)
        self.red_trim.textChanged.connect(self.update_red_trim)
        # -----------------------------------------------------------------

        self.low_sigma_clip.textChanged.connect(self.check_state)
        self.high_sigma_clip.textChanged.connect(self.check_state)
        self.knot_spacing.textChanged.connect(self.check_state)
        # -----------------------------------------------------------------
        # E. Holmbeck added connect to trim
        self.blue_trim.textChanged.connect(self.check_state)
        self.red_trim.textChanged.connect(self.check_state)
        # -----------------------------------------------------------------

        return None


    def normalize_and_stitch(self):
        """
        Normalize any remaining orders, and stitch them together.
        """

        # Normalize any remaining orders.
        index = self.current_order_index 
        for i in range(len(self.parent.session.input_spectra)):
            self.update_order_index(i) # TODO: This is clumsy.
            self.fit_continuum(clobber=False)

        # Go back to original order.
        self.update_order_index(index)

        # Stitch and stack all orders.
        self.parent.session.stitch_and_stack()

        # Enable the menu-bar and the next three tabs.
        self.parent._menu_export_normalized_spectrum.setEnabled(True)
        self.parent.tabs.setTabEnabled(self.parent.tabs.indexOf(self) + 1, True)
        self.parent.tabs.setTabEnabled(self.parent.tabs.indexOf(self) + 2, True)
        self.parent.tabs.setTabEnabled(self.parent.tabs.indexOf(self) + 3, True)

        self.parent.stellar_parameters_tab.new_session_loaded()
        self.parent.chemical_abundances_tab.new_session_loaded()

        return None


    def check_state(self, *args, **kwargs):
        """
        Update the background color of a QLineEdit object based on whether the
        input is valid.
        """

        # TODO: Implement from
        # http://stackoverflow.com/questions/27159575/pyside-modifying-widget-colour-at-runtime-without-overwriting-stylesheet

        sender = self.sender()
        validator = sender.validator()
        state = validator.validate(sender.text(), 0)[0]
        if state == QtGui2.QValidator.Acceptable:
            color = 'none' # normal background color
        elif state == QtGui2.QValidator.Intermediate:
            color = '#fff79a' # yellow
        else:
            color = '#f6989d' # red
        sender.setStyleSheet('QLineEdit { background-color: %s }' % color)


    def figure_key_press(self, event):
        """
        Key press event in the normalization figure.
        """

        # Check if we are waiting for an exclusion region first
        # (This means the mouse has been clicked, but not released in-axis yet)
        try:
            self._exclude_selected_region_signal
        except AttributeError:
            None
        else:
            return None

        # Show a new order.
        if event.key in ("left", "right"):
            offset = 1 if event.key == "right" else -1
            self.order_slide.setValue(self.order_slide.value() + offset)
            self.update_order_figure()

            return None

        # Scale the continuum up/down.
        if event.key in ("up", "down"):
            """
            clip = self._cache["input"]["high_sigma_clip"]
            if event.key == "up":
                clip = max(clip-0.01, 0)
            if event.key == "down":
                clip += 0.01
            self._cache["input"]["high_sigma_clip"] = clip
            self.high_sigma_clip.setText(
                str(self._cache["input"]["high_sigma_clip"]))
            
            """
            scale = self._cache["input"].get("scale", 1.0)
            sign = +1 if event.key == "up" else -1

            self._cache["input"]["scale"] = scale + sign * 0.01

            self.fit_continuum(True)
            self.draw_continuum(True)

            return None


        # 'd': No normalization for this order.
        if event.key in ("d", "D"):
            try:
                idx, session = self.current_order_index, self.parent.session

            except AttributeError:
                return None

            session.metadata["normalization"]["continuum"][idx] = 1
            session.metadata["normalization"]["normalization_kwargs"][idx] = {}

            self.draw_continuum(True)

            return None


        # 'c': Clear the scale, excluded regions and points for this order.
        if event.key in "cC":
            for key in ("scale", "exclude", "additional_points"):
                if key in self._cache["input"]:
                    del self._cache["input"][key]

            # Force refit.
            self.fit_continuum(clobber=True)
            self.draw_continuum(refresh=False)
            self.update_continuum_mask(refresh=True)
            self.norm_plot.reset_zoom_limits()

            return True


        # undo/remove the last mask
        if event.key in ("u", "U"):
            if "exclude" in self._cache["input"]:
                exclude_regions = self._cache["input"]["exclude"]
                if len(exclude_regions) > 0:
                    exclude_regions = exclude_regions[1:]
                    self._cache["input"]["exclude"] = exclude_regions
                    
                    self.fit_continuum(clobber=True)
                    self.draw_continuum(refresh=False)
                    self.update_continuum_mask(refresh=True)
        
        # 'r': Reset the zoom limits without refitting/clearing masks
        if event.key in "rR":
            self.norm_plot.reset_zoom_limits()
            self.draw_continuum(refresh=False)
            self.update_continuum_mask(refresh=True)

            return True


        # 'f': Refit without resetting the zoom limits.
        # Also can be used to recenter the bottom plot
        if event.key in "fF":
            # Force refit.
            self.fit_continuum(clobber=True)
            self.draw_continuum(refresh=False)
            self.update_continuum_mask(refresh=True)

            return True
        
        # 'a': add point
        if event.key in "aA":
            points = np.vstack([
                self.ax_order.collections[0].get_offsets(),
                [event.xdata, event.ydata]
            ])
            # TODO: set size by their weight?
            self.ax_order.collections[0].set_offsets(points)
            
            idx = self.current_order_index
            N = points.shape[0]
            # TODO: adhere to the knot weights
            self._cache["input"]["additional_points"] \
                = np.hstack((points, 100 * np.ones(N).reshape((N, 1))))

            self.fit_continuum(clobber=True)
            self.draw_continuum(refresh=False)
            self.update_continuum_mask(refresh=True)
            return True
            
        # 'x': clear all added points
        if event.key in "xX":
            for key in ["additional_points"]:
                if key in self._cache["input"]:
                    del self._cache["input"][key]
            
            self.fit_continuum(clobber=True)
            self.draw_continuum(refresh=False)
            self.update_continuum_mask(refresh=True)
            return True
            
    def figure_mouse_press(self, event):
        """
        Function to handle event left clicks (single or double click).
        
        :param event:
            The matplotlib event signal.
        """
        
        # Add/remove an additional point?
        if event.dblclick:

            logger.info("Removed double-click to add point. Use 'a' key instead")
            if event.button == 1:
                # Add a point.
                # Removed adding points because APJ strongly recommends not using this
                """
                points = np.vstack([
                    self.ax_order.collections[0].get_offsets(),
                    [event.xdata, event.ydata]
                ])
                # TODO: set size by their weight?
                self.ax_order.collections[0].set_offsets(points)
                """
                pass
            
            else:
                # Are we within <tolerance of a point?
                points = self.ax_order.collections[0].get_offsets()

                # Need to scale x-distance to convert to pixels.
                idx = self.current_order.dispersion.searchsorted(event.xdata)
                xscale = np.nanmean(
                    np.diff(self.current_order.dispersion[idx-5:idx+5]))

                """
                bbox = self.ax_order.get_window_extent().transformed(
                    self.norm_plot.dpi_scale_trans.inverted())
                width =  bbox.width * self.norm_plot.dpi
                height = bbox.height * self.norm_plot.dpi
                print(width, height)
                """
                # TODO: Fix this distance thing.

                distance = np.sqrt(
                      ((points[:, 0] - event.xdata)/xscale)**2 \
                    + (points[:, 1] - event.ydata)**2)
                
                if distance.size > 0:

                    index = np.argmin(distance)
                    if distance[index] < PIXEL_PICKER_TOLERANCE:
                        # Remove that point.
                        keep = np.ones(points.shape[0], dtype=bool)
                        keep[index] = False
                        self.ax_order.collections[0].set_offsets(points[keep])

                    else:
                        print("Closest point {} px away".format(distance[index]))

            # Update the cache.
            """
            # Removed here because don't want to add points this way
            idx = self.current_order_index
            N = points.shape[0]
            # TODO: adhere to the knot weights
            self._cache["input"]["additional_points"] \
                = np.hstack((points, 100 * np.ones(N).reshape((N, 1))))
            """
            self.fit_continuum(clobber=True)
            self.draw_continuum(refresh=True)

            return None
        
        if event.button != 1: return None
        # Single click.
        # Set up/update the excluded region.
        xmin, xmax, ymin, ymax = (event.xdata, np.nan, -1e8, +1e8)
        try:
            self._exclude_selected_region
        except AttributeError:
            self._exclude_selected_region = self.ax_order.axvspan(**{
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
            self._exclude_selected_region.set_xy([
                [xmin, ymin],
                [xmin, ymax],
                [xmax, ymax],
                [xmax, ymin],
                [xmin, ymin]
            ])

        # Set the signal and the time.
        self._exclude_selected_region_signal = (
            time(),
            self.norm_plot.mpl_connect(
                "motion_notify_event", self.update_exclude_selected_region)
            )
        return None



    def figure_mouse_release(self, event):
        """
        A callback function that is executed when the left mouse button is released.
        This signal typically triggers the extent of a region to mask in the
        current order.

        :param event:
            The matplotlib event.
        """

        if event.button != 1: return None
        try:
            signal_time, signal_cid = self._exclude_selected_region_signal
        except AttributeError:
            return None
        
        xy = self._exclude_selected_region.get_xy()
        
        if xy[0,0] > xy[2,0]:
            x1, x2 = xy[0,0], xy[2,0]
            xy[0,0] = xy[1,0] = xy[4,0] = x2
            xy[2,0] = xy[3,0] = x1
        
        # If the two mouse events were within some time interval,
        # then we should not add a mask because those signals were probably
        # part of a double-click event.
        if  time() - signal_time > DOUBLE_CLICK_INTERVAL \
        and np.abs(xy[0,0] - xy[2,0]) > 0:
            
            # Update the cache with the new mask.
            _ =  np.array(self._cache["input"].get("exclude", np.array([])))
            _.shape = (-1, 2)
            
            self._cache["input"]["exclude"] = np.vstack((
                np.array([xy[0,0], xy[2, 0]]).reshape(-1, 2), _))
            
            
            #print("exclude",self._cache["input"]["exclude"])

            # Fit and re-draw the continuum, and its mask.
            self.fit_continuum(clobber=True)
            self.update_continuum_mask(refresh=False)
            self.draw_continuum(refresh=False)

        xy[:, 0] = np.nan

        self._exclude_selected_region.set_xy(xy)
        self.norm_plot.mpl_disconnect(signal_cid)
        self.norm_plot.draw()
        del self._exclude_selected_region_signal
        return None


    def update_exclude_selected_region(self, event):
        """
        Update the visible selected exclusion region for this order. This
        function is linked to a callback for when the mouse position moves.

        :param event:
            The matplotlib motion event to show the current mouse position.
        """
        if event.xdata is None:
            return
        
        signal_time, signal_cid = self._exclude_selected_region_signal
        if time() - signal_time > DOUBLE_CLICK_INTERVAL: 
            
            data = self._exclude_selected_region.get_xy()

            # Update xmax.
            data[2:4, 0] = event.xdata
            self._exclude_selected_region.set_xy(data)

            self.norm_plot.draw()

        return None


    def _populate_widgets(self):
        """
        Populate the widgets with default settings, and a default figure view.
        """

        if self.parent.session is None:
            # No point populating the widgets with the default values from the
            # SMH file because these will be updated when a session is loaded.
            return

        keys = ("function", "order", "low_sigma_clip", "high_sigma_clip",
            "knot_spacing", "blue_trim", "red_trim", "max_iterations")
        self._cache = {
            "input": {}
        }
        for key in keys:
            self._cache["input"][key] \
                = self.parent.session.setting(("normalization", key))

        # Continuum masks.
        self._cache["masks"] \
            = self.parent.session.setting(("normalization", "masks"))
        self._cache["default_mask"] \
            = self.parent.session.setting(("normalization", "default_mask")) \
                or self._cache["masks"].keys()[0]


        # Put these values into the widgets.
        self.low_sigma_clip.setText(
            str(self._cache["input"]["low_sigma_clip"]))
        self.high_sigma_clip.setText(
            str(self._cache["input"]["high_sigma_clip"]))
        self.knot_spacing.setText(str(
            self._cache["input"]["knot_spacing"]))
        # ----------------------------------------------------------------
        # E. Holmbeck added these
        self.blue_trim.setText(str(
            self._cache["input"]["blue_trim"]))
        self.red_trim.setText(str(
            self._cache["input"]["red_trim"]))
        # ----------------------------------------------------------------

        functions = [self.function.itemText(i).lower() \
            for i in range(self.function.count())]
        self.function.setCurrentIndex(functions.index(
            self._cache["input"]["function"]))

        # Normalization order.
        orders = [int(self.order.itemText(i)) \
            for i in range(self.order.count())]
        self.order.setCurrentIndex(orders.index(
            self._cache["input"]["order"]))

        # Normalization maximum iterations.
        norm_max_iters = [int(self.norm_max_iter.itemText(i)) \
            for i in range(self.norm_max_iter.count())]
        self.norm_max_iter.setCurrentIndex(norm_max_iters.index(
            self._cache["input"]["max_iterations"]))

        # Mask names.
        for name in self._cache["masks"].keys():
            self.continuum_mask.addItem(name)

        self.continuum_mask.setCurrentIndex(
            list(self._cache["masks"].keys()).index(
                self._cache["default_mask"]))

        self.order_slide.setMaximum(len(self.parent.session.input_spectra) - 1)
        self.current_order_label.setText("Order 1 of {}".format(
            len(self.parent.session.input_spectra)))

        # Draw the widgets.
        try:
            self.order_slide.setValue(0)
            self.update_order_index(0)
            self.update_continuum_mask(refresh=False)
            self.fit_continuum(clobber=False)
            self.draw_order(refresh=False)
            self.draw_continuum(refresh=True)

        except (AttributeError, KeyError):
            # HACK
            # when loading a fresh session, it will skip all those blocks
            # I think this is okay?
            pass
        return None


    def update_rv_applied(self):
        """
        Make updates to the view when the radial velocity applied has been
        updated.
        """

        # Clear out the normalization in any other orders.
        N = len(self.parent.session.input_spectra)
        self.parent.session.metadata["normalization"] = {
            "continuum": [None] * N,
            "normalization_kwargs": [{}] * N
        }

        # Update the current order fit, and the view.
        self.update_order_index()
        self.update_continuum_mask(refresh=False)
        self.fit_continuum(clobber=True)
        self.draw_order(refresh=False)
        self.draw_continuum(refresh=True)

        return None



    def update_continuum_mask(self, refresh=False):
        """
        Draw the continuum mask (relevant for all orders).
        """

        ymin, ymax = (-1e8, 1e8)
        kwds = {
            "xmin": np.nan,
            "xmax": np.nan,
            "ymin": ymin,
            "ymax": ymax,
            "facecolor": "r",
            "edgecolor": "none",
            "alpha": 0.25,
            "zorder": -1
        }

        transform = lambda start, end, v=0: np.array([
                [start * (1 - v/c), ymin],
                [start * (1 - v/c), ymax],
                [end   * (1 - v/c), ymax],
                [end   * (1 - v/c), ymin],
                [start * (1 - v/c), ymin]
            ])

        mask = self._cache["masks"][self.continuum_mask.currentText()]

        # Any added regions to mask out? v-stack these
        try:
            self._masked_wavelengths
        except AttributeError:
            self._masked_wavelengths = []
            self._masked_wavelengths_norm = []

        # Different kind of masks: rest_wavelength, obs_wavelength, pixels
        # rest_wavelength
        # The obsered spectrum is shifted to be at rest, so the continuum masks
        # will also be in the rest frame. So we don't need to shift the
        # 'rest_wavelength' mask, but we do need to shift the 'obs_wavelength'
        # mask

        # Get the applied velocity to shift some masks.
        try:
            rv_applied = self.parent.session.metadata["rv"]["rv_applied"]
        except (AttributeError, KeyError):
            rv_applied = 0
        
        # -----------------------------------------------------------------
        # E. Holmbeck added read-in BCV from header
        try:
            vhelio = self.parent.session.metadata["rv"]["heliocentric_correction"]
            bcv_shift = self.parent.session.metadata["rv"]["barycentric_correction"]
            dop_shift = vhelio + bcv_shift
        except (AttributeError, KeyError):
            dop_shift = 0.0
        # -----------------------------------------------------------------
        _ =self.parent.session.metadata["normalization"]["normalization_kwargs"]
        masked_regions = [
            #np.array(mask.get("rest_wavelength", [])),
            np.array(mask.get("rest_wavelength", [])) * (1.0 - dop_shift/c),
            np.array(mask.get("obs_wavelength", [])) * (1.0 - rv_applied/c),
            np.array(_[self.current_order_index].get("exclude", []))
        ]
        if "pixel" in mask:
            masked_regions.append(
                # MAGIC HACK
                self.current_order.dispersion[np.array(mask["pixel"])] + 1e-3
            )

        for each in masked_regions:
            each.shape = (-1, 2)

        masked_regions = np.vstack(masked_regions)

        # Remove duplicate masked regions.
        _ = np.ascontiguousarray(masked_regions).view(
            np.dtype((
                np.void, 
                masked_regions.dtype.itemsize * masked_regions.shape[1])))
        __, idx = np.unique(_, return_index=True)
        masked_regions = masked_regions[idx]

        i = 0
        for start, end in masked_regions:
            if i >= len(self._masked_wavelengths):
                # Create a polygon in the main axis.
                self._masked_wavelengths.append(
                    self.ax_order.axvspan(**kwds))

                # And for the normalization axis.
                self._masked_wavelengths_norm.append(
                    self.ax_order_norm.axvspan(**kwds))

            polygons = (
                self._masked_wavelengths[i],
                self._masked_wavelengths_norm[i]
            )
            for polygon in polygons:
                polygon.set_xy(transform(start, end))

            i += 1

        # Any leftover polygons?
        for polygon in self._masked_wavelengths[i:]:
            polygon.set_xy(transform(np.nan, np.nan))

        for polygon in self._masked_wavelengths_norm[i:]:
            polygon.set_xy(transform(np.nan, np.nan))


        if refresh:
            self.norm_plot.draw()
        return True



    def update_knot_spacing(self):
        """ Update the knot spacing. """
        knot_spacing = self.knot_spacing.text()
        if knot_spacing:
            self._cache["input"]["knot_spacing"] = float(knot_spacing)
            self.reset_input_style_defaults()
            self.fit_continuum(True, sender=self.knot_spacing)
            self.draw_continuum(True)
            
        return None

	# -----------------------------------------------------------------
	# E. Holmbeck added these update functions
    def update_blue_trim(self):
        try:
        	trim_region = int(self.blue_trim.text())
        except ValueError:
        	return None
        
        if trim_region == 0:
            return None
            
        try:
            x, y = (self.current_order.dispersion, self.current_order.flux)
        except AttributeError:
            return None
        
        # If "exclude" doesn't exist, add it.
        try:
            exclude = self._cache["input"]["exclude"]
        except:
            self._cache["input"]["exclude"] = np.array( 
                [[x[0], x[trim_region]]])
            exclude = self._cache["input"]["exclude"]
        
        # Replace the mask that goes to the end anyway     
        if len(exclude) == 0:
            self._cache["input"]["exclude"] = np.array( 
                [[x[0], x[trim_region]]])
        else:
            mask_to_delete = []
            for i,e in enumerate(exclude):
                if e[0] == x[0] and e[1] != x[trim_region]:
                    mask_to_delete.append(i)
            
            self._cache["input"]["exclude"] = np.delete(exclude, mask_to_delete, axis=0)
            self._cache["input"]["exclude"] = np.insert(self._cache["input"]["exclude"], 
                0, [[x[0], x[trim_region]]], axis=0)

        if trim_region:        
            self._cache["input"]["blue_trim"] = trim_region
            self.reset_input_style_defaults()
            self.fit_continuum(True)
            self.draw_continuum(True)
            self.update_continuum_mask(refresh=True)

        return None

    def update_red_trim(self):
        try:
        	trim_region = int(self.red_trim.text())
        except ValueError:
        	return None

        if trim_region == 0:
            return None
        
        try:
            x, y = (self.current_order.dispersion, self.current_order.flux)
        except AttributeError:
            return None
        
        try:
            exclude = self._cache["input"]["exclude"]
        except:
            self._cache["input"]["exclude"] = np.array( 
                [[x[-trim_region], x[-1]+1e-3]])
            exclude = self._cache["input"]["exclude"]
        
        if len(exclude) == 0:
            self._cache["input"]["exclude"] = np.array( 
                [[x[-trim_region], x[-1]+1e-3]])
        else:
            mask_to_delete = []
            for i,e in enumerate(exclude):
                if e[0] != x[-trim_region] and e[1] == x[-1]+1e-3:
                    mask_to_delete.append(i)
            
            self._cache["input"]["exclude"] = np.delete(exclude, mask_to_delete, axis=0)
            self._cache["input"]["exclude"] = np.append(self._cache["input"]["exclude"], 
                [[x[-trim_region], x[-1]+1e-3]], axis=0)

        if trim_region:        
            self._cache["input"]["red_trim"] = trim_region
            self.reset_input_style_defaults()
            self.fit_continuum(True)
            self.draw_continuum(True)
            self.update_continuum_mask(refresh=True)

        return None
	# -----------------------------------------------------------------

    def update_high_sigma_clip(self):
        """ Update the high sigma clip value. """
        high_sigma = self.high_sigma_clip.text()
        if high_sigma:
            try:
                self._cache["input"]["high_sigma_clip"] = float(high_sigma)
            except ValueError:
                pass
            self.reset_input_style_defaults()
            self.fit_continuum(True)
            self.draw_continuum(True)
        return None


    def update_low_sigma_clip(self):
        """ Update the low sigma clip value. """
        low_sigma = self.low_sigma_clip.text()
        if low_sigma:
            try:
                self._cache["input"]["low_sigma_clip"] = float(low_sigma)
            except ValueError:
                pass
            self.reset_input_style_defaults()
            self.fit_continuum(True)
            self.draw_continuum(True)
        return None


    def update_normalization_function(self, refresh=True):
        """ Update the normalization function. """
        self._cache["input"]["function"] = self.function.currentText()

        indices = range(5, 10)

        # If the function is a spline, then it is limited to order 5.
        if self._cache["input"]["function"] == "Spline":
            if int(self.order.currentText()) > 5:
                # Limit it to 5.
                self.order.setCurrentIndex(4) # Index 4 = Order '5'

            # Disable the other entries.
            for i in indices:
                item = self.order.model().item(i)
                if item is not None:
                    item.setEnabled(False)

        else:
            # Enable order entries greater than 5.
            for i in indices:
                item = self.order.model().item(i)
                if item is not None:
                    item.setEnabled(True)


        self.reset_input_style_defaults()
        if refresh:
            self.fit_continuum(True)
            self.draw_continuum(True)
        return None


    def update_normalization_order(self):
        """ Update the normalization order. """
        self._cache["input"]["order"] = int(self.order.currentText())
        self.reset_input_style_defaults()
        self.fit_continuum(True)
        self.draw_continuum(True)
        return None


    def update_normalization_max_iterations(self):
        """ Update the maximum number of iterations. """
        self._cache["input"]["max_iterations"] \
            = int(self.norm_max_iter.currentText())
        self.fit_continuum(True)
        self.draw_continuum(True)
        self.reset_input_style_defaults()
        return None


    def update_order_index(self, index=None):
        """
        Update the currently selected order.
        """
        if index is None:
            index = getattr(self, "current_order_index", 0)

        session = self.parent.session
        self.current_order_index = index
        self.current_order \
            = session.input_spectra[self.current_order_index].copy()

        # Apply any RV correction.
        try:
            v = session.metadata["rv"]["rv_applied"]
        except (AttributeError, KeyError):
            v = 0

        self.current_order._dispersion *= (1 - v/c)

        # Update the view if the input settings don't match the settings used
        # to normalize the current order.
        self.check_for_different_input_settings()
        
        return None



    def update_order_figure(self, index=None):

        #offset = 1 if event.key == "right" else -1

        self.update_order_index(index or self.order_slide.value())

        self.draw_order(refresh=False)
        self.update_continuum_mask(refresh=False)
        self.fit_continuum(clobber=False)
        self.draw_continuum(refresh=True)

        self.current_order_label.setText("Order {} of {}".format(
            1 + self.order_slide.value(), len(self.parent.session.input_spectra)))
        
        # Do all of the orders have continuum? If so, update the button.
        if not any([(_ is None) for _ in \
        self.parent.session.metadata["normalization"]["continuum"]]):
            self.stitch_btn.setText("Stitch orders")

        return None 











    def check_for_different_input_settings(self):
        """
        Check whether the current input settings reflect the settings used to
        normalize the currently displayed order.

        2023-07-10 APJ made it so it just updates the display, instead of highlighting
        """

        session, index = self.parent.session, self.current_order_index

        # Is there continuum already for this new order?
        try:
            continuum = session.metadata["normalization"]["continuum"][index]
        except:
            logger.debug("check_for_different_input_settings: no continuum kwd")
        normalization_kwargs \
            = session.metadata["normalization"]["normalization_kwargs"][index]

        # These keys don't have widgets, but need to be updated.
        extra_keys = ("additional_points", "exclude")
        for key in extra_keys:
            if key in normalization_kwargs:
                self._cache["input"][key] = normalization_kwargs[key]
            elif key in self._cache["input"]:
                del self._cache["input"][key]

        if continuum is None:
            # Holmbeck: a little hacky, but it works?
            self.update_blue_trim()
            self.update_red_trim()
            return

        # If so, are the current normalization keywords different to the ones
        # used for this one?
        input_items = {
            "function": [self.function_label, self.function],
            "order": [self.order_label, self.order],
            "max_iterations": [self.max_iter_label, self.norm_max_iter],
            "low_sigma_clip": [self.low_sigma_clip_label, self.low_sigma_clip],
            "high_sigma_clip": \
                [self.high_sigma_clip_label, self.high_sigma_clip],
            "knot_spacing": [self.knot_spacing, self.knot_spacing_label],
            "blue_trim": [self.blue_trim, self.blue_trim_label],
            "red_trim": [self.red_trim, self.red_trim_label],
        }

        diff = dict_updated(self._cache["input"], normalization_kwargs,
            exclude=("additional_points", "exclude"))

        # update the cache values
        self._cache["input"].update(normalization_kwargs)
        # update the view for all textboxes
        self.low_sigma_clip.setText(str(self._cache["input"]["low_sigma_clip"]))
        self.high_sigma_clip.setText(str(self._cache["input"]["high_sigma_clip"]))
        self.knot_spacing.setText(str(self._cache["input"]["knot_spacing"]))
        self.blue_trim.setText(str(self._cache["input"]["blue_trim"]))
        self.red_trim.setText(str(self._cache["input"]["red_trim"]))
        # update the view for dropdown boxes
        functions = [self.function.itemText(i).lower() \
            for i in range(self.function.count())]
        orders = [int(self.order.itemText(i)) \
            for i in range(self.order.count())]
        norm_max_iters = [int(self.norm_max_iter.itemText(i)) \
            for i in range(self.norm_max_iter.count())]
        #print(functions, orders, norm_max_iters)
        self.function.setCurrentIndex(functions.index(
            self._cache["input"]["function"].lower()))
        self.order.setCurrentIndex(orders.index(
            self._cache["input"]["order"]))
        self.norm_max_iter.setCurrentIndex(norm_max_iters.index(
            self._cache["input"]["max_iterations"]))
        ## Update the valid order numbers in the GUI
        self.update_normalization_function(False)
        
        # By default, everything should be styled normally.
        self.reset_input_style_defaults(sum(input_items.values(), []))

        """
        for key, (current, used) in diff.items():
            if key in input_items:
                # Update the font-weight of those objects.
                items = input_items[key]
                for item in items:
                    item.setStyleSheet("{0} {{ font-weight: bold }}".format(
                        item.__class__.__name__))
                    item.setStatusTip("Order {0} was normalized using {1} ="
                        " {2} (not {3})"\
                        .format(1 + index, key, used, current))
        """

            
        return None


    def reset_input_style_defaults(self, items=None):
        """
        Reset the styling inputs.
        """
        items = items or (
            self.function_label, self.function,
            self.order_label, self.order,
            self.max_iter_label, self.norm_max_iter,
            self.low_sigma_clip_label, self.low_sigma_clip,
            self.high_sigma_clip_label, self.high_sigma_clip,
            self.knot_spacing_label, self.knot_spacing,
            self.blue_trim_label, self.blue_trim,
            self.red_trim_label, self.red_trim,
        )
        # Ensure all the things are styled normally.
        for item in items:
            item.setStyleSheet('{0} {{ font-weight: normal }}'\
                .format(item.__class__.__name__))
            item.setStatusTip("")
        self.parent.statusbar.showMessage(self.parent._default_statusbar_message)
        return None


    def draw_order(self, refresh=False):
        """
        Draw the current order.

        Note: The order in `self.current_order` is already in the rest frame. 
        """

        x, y = (self.current_order.dispersion, self.current_order.flux)
        
        self.ax_order.lines[0].set_data([x, y])
        
        # Show a few percent either side.
        percent = 2
        trimming = (x[-1] - x[0]) * percent/100.
        self.ax_order.set_xlim(x[0] - trimming, x[-1] + trimming)

        trimming = (np.nanmax(y) - np.nanmin(y)) * percent/100.
        self.ax_order.set_ylim(np.nanmin(y) - trimming, np.nanmax(y) + trimming)

        self.norm_plot.reset_zoom_limits()

        #self.ax_order.set_title("Order {0} of {1}".format(
        #    1 + self.current_order_index, 
        #    len(self.parent.session.input_spectra)))

        if refresh:
            self.norm_plot.draw()
        return None


    def fit_continuum(self, clobber, sender=None):
        """
        Update continuum for the current order.

        :param clobber:
            Clobber any existing continuum determination for this order.
        """

        # Any existing continuum determination?
        try:
            index, session = (self.current_order_index, self.parent.session)
        except AttributeError:
            return None

        try:
            continuum = session.metadata["normalization"]["continuum"][index]
        except KeyError:
            # Nothing to do
            return
        if continuum is not None and not clobber:
            # Nothing to do.
            return
        #print("gui.normalization.fit_continuum: fitting {}".format(index))

        kwds = self._cache["input"].copy()
        kwds["full_output"] = True

        # Add in the global continuum masks.
        global_mask = self._cache["masks"][self.continuum_mask.currentText()]
        try:
            rv_applied = self.parent.session.metadata["rv"]["rv_applied"]
        except (AttributeError, KeyError):
            rv_applied = 0
        
        # -----------------------------------------------------------------
        # E. Holmbeck added read-in BCV from header
        try:
            vhelio = self.parent.session.metadata["rv"]["heliocentric_correction"]
            bcv_shift = self.parent.session.metadata["rv"]["barycentric_correction"]
            dop_shift = vhelio + bcv_shift
        except (AttributeError, KeyError):
            dop_shift = 0.0
        # -----------------------------------------------------------------

        if np.isnan(dop_shift):
            mask_kinds = [
                (0,  global_mask.get("rest_wavelength", [])),
                (rv_applied, global_mask.get("obs_wavelength", []))
            ]
        else:
            mask_kinds = [
                (dop_shift,  global_mask.get("rest_wavelength", [])),
                (rv_applied, global_mask.get("obs_wavelength", []))
            ]

        regions = []
        for v, masked_regions in mask_kinds:
            for region in np.array(masked_regions):
                start, end = np.array(region) * (1 - v/c)

                if end >= self.current_order.dispersion[0] \
                and self.current_order.dispersion[-1] >= start:
                    regions.append((start, end))

        if "pixel" in global_mask:
            # MAGIC HACK
            regions.extend(1e-3 + \
                self.current_order.dispersion[np.array(global_mask["pixel"])])

        if kwds.get("exclude", None) is None:
            kwds["exclude"] = np.array(regions)

        else:
            _ = np.array(kwds["exclude"])
            _.shape = (-1, 2)
            kwds["exclude"] = np.vstack((_, np.array(regions).reshape(-1, 2)))

        try:
            normalized_spectrum, continuum, left, right \
                = self.current_order.fit_continuum(**kwds)
            #if kwds["knot_spacing"] == 201:
            #    raise KeyError("what") #HACK #TESTING #TODO
        except:
            logger.exception("No continuum could be fit.")
            self.parent.statusbar.showMessage(
                "Exception occurred while trying to fit the continuum.")
            raise
            continuum = np.nan

            # Did a user input something bad? Let them know..
            if sender is not None:
                print("Sender is ", sender) #TODO: This is broken..
                sender.setStyleSheet("{0} {{ background-color: #f6989d; }}"\
                    .format(sender.__class__.__name__))
        
        session.metadata["normalization"]["continuum"][index] = continuum
        session.metadata["normalization"]["normalization_kwargs"][index] = kwds

        return None


    def draw_continuum(self, refresh=False):
        """
        Draw the continuum for the current order.

        Note: The order in `self.current_order` is already in the rest-frame.
        """

        try:
            index = self.current_order_index
        except AttributeError:
            return None

        meta = self.parent.session.metadata["normalization"]
        try:
            continuum = meta["continuum"][index]
        except:
            logger.debug("draw_continuum: no continuum kw")
            return None
        kwds = meta["normalization_kwargs"][index]

        self.ax_order.lines[1].set_data([
            self.current_order.dispersion, continuum])

        # Set the color of the line.
        #self.ax_order.lines[1].set_color(wavelength_to_hex(
        #    np.nanmean(self.current_order.dispersion)))
        self.ax_order.lines[1].set_color("r")
        
        # Add points showing knot locations
        try:
            knots = self.current_order.get_knots(kwds["knot_spacing"], kwds["exclude"])
            knoty = continuum[np.searchsorted(self.current_order.dispersion, knots)]
            self.ax_order.lines[2].set_data([knots, knoty])
        except KeyError:
            self.ax_order.lines[2].set_data([[np.nan], [np.nan]])
        
        # Update the normalization preview in the lower axis.
        self.ax_order_norm.lines[1].set_data([
            self.current_order.dispersion, 
            self.current_order.flux/continuum])
        self.ax_order_norm.set_xlim(self.ax_order.get_xlim())

        # Draw the additional points.
        ap = meta["normalization_kwargs"][index].get("additional_points", None)
        if ap is None:
            ap = np.array([[],[]]).T
        else:
            ap = ap[:, :2]
        self.ax_order.collections[0].set_offsets(ap)

        if refresh:
            self.norm_plot.draw()
        return None

