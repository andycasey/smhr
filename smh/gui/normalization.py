#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The normalization tab in Spectroscopy Made Hard """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["NormalizationTab"]

import logging
import numpy as np
import sys
from PySide import QtCore, QtGui
from time import time

from matplotlib import gridspec

import mpl

# This is a bad idea, but it's April 1st.
from style_utils import wavelength_to_hex


logger = logging.getLogger(__name__)

c = 299792458e-3 # km/s

# The minimum time (in seconds) between a mouse click/release to differentiate
# a single click from a double click
DOUBLE_CLICK_INTERVAL = 0.1 # MAGIC HACK

# The pixel tolerance to select and remove an additional point.
PIXEL_PICKER_TOLERANCE = 30 # MAGIC HACK


def dict_updated(default, new):
    updated = {}
    for key, value in default.items():
        if key in new and new[key] != value:
            updated[key] = (value, new[key])
    return updated


class NormalizationTab(QtGui.QWidget):

    def __init__(self, parent=None):
        super(NormalizationTab, self).__init__(parent)
        self.parent = parent

        # Establish the GUI for this tab.
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)

        # Create a top-level horizontal layout to contain a MPL figure and
        # a vertical layout of settings..
        tab_layout = QtGui.QHBoxLayout(self)
        tab_layout.setContentsMargins(10, 10, 10, 10)
        
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
            QtGui.QDoubleValidator(0, 1000, 2, self.low_sigma_clip))

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
            QtGui.QDoubleValidator(0, 1000, 2, self.high_sigma_clip))
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
            QtGui.QIntValidator(0, 10000, self.knot_spacing))
        self.knot_spacing.setObjectName("norm_knot_spacing")
        hbox.addWidget(self.knot_spacing)
        settings_grid_layout.addLayout(hbox, 5, 1, 1, 1)


        # End of the grid in the normalization tab.
        settings_layout.addLayout(settings_grid_layout)

        # Add a label.
        label = QtGui.QLabel(self)
        label.setText("Continuum mask:")
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


        # Add a spacer.
        settings_layout.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding))

        # Add a 'normalize and stitch button'
        self.stitch_btn = QtGui.QPushButton(self)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(self.stitch_btn.sizePolicy().hasHeightForWidth())
        self.stitch_btn.setSizePolicy(sp)
        self.stitch_btn.setMinimumSize(QtCore.QSize(300, 0))
        self.stitch_btn.setMaximumSize(QtCore.QSize(300, 16777215))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.stitch_btn.setFont(font)
        self.stitch_btn.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.stitch_btn.setDefault(True)
        self.stitch_btn.setObjectName("stitch_btn")
        self.stitch_btn.setText("Normalize and stitch orders")
        if sys.platform == "darwin":
            self.stitch_btn.setStyleSheet('QPushButton {color: white}')

        settings_layout.addWidget(self.stitch_btn)


        tab_layout.addWidget(settings_widget)

        # Create a matplotlib widget.
        blank_widget = QtGui.QWidget(self)
        blank_widget.setStatusTip(
            "Use left/right arrows to move between orders, "
            "double-click for all keyboard shortcuts")
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(blank_widget.sizePolicy().hasHeightForWidth())
        blank_widget.setSizePolicy(sp)
        blank_widget.setObjectName("norm_plot")


        self.norm_plot = mpl.MPLWidget(blank_widget, tight_layout=True,
            autofocus=True)

        layout = QtGui.QVBoxLayout(blank_widget)
        layout.addWidget(self.norm_plot, 1)
        tab_layout.addWidget(blank_widget)

        # Set up the plot.
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        self.ax_order = self.norm_plot.figure.add_subplot(gs[0])
        # Line for the data.
        self.ax_order.plot([], [], c='k', zorder=3)#, drawstyle='steps-mid')
        # Line for the continuum.
        self.ax_order.plot([], [], linestyle="--", linewidth=2, c='r', zorder=4)

        # Line for the neighbouring order(s) (joined by a NaN).
        self.ax_order.plot([], [], c='#666666', zorder=1, drawstyle='steps-mid')
        # Line for the neighbouring order(s) continuum (joined by a NaN)
        self.ax_order.plot([], [], c='b', zorder=2)

        # Additional point markers.
        self.ax_order.scatter([], [], facecolor="k", zorder=5, picker=5)

        # Regions

        self.ax_order.set_xticklabels([])
        self.ax_order.set_yticklabels([])
        self.ax_order.set_ylabel("Flux")

        self.ax_order_norm = self.norm_plot.figure.add_subplot(gs[1])
        self.ax_order_norm.axhline(1, linestyle=":", c="#666666", zorder=1)
        self.ax_order_norm.plot([], [], c='k', zorder=2)
        self.ax_order_norm.set_ylim(0, 1.2)
        self.ax_order_norm.set_yticks([0, 0.5, 1.0])
        self.ax_order_norm.set_xlabel(u"Wavelength (Å)")

        #self.norm_plot.figure.tight_layout(w_pad=0, h_pad=0, pad=0.4)
        self.norm_plot.draw()

        # Create signals.
        # Note that key_press_event is linked to norm_plot.canvas, while the
        # mouse events are linked to norm_plot.
        # I don't know why, but that's how it works.
        self.norm_plot.canvas.mpl_connect(
            "key_press_event", self.figure_key_press)
        self.norm_plot.mpl_connect(
            "button_press_event", self.figure_mouse_press)
        self.norm_plot.mpl_connect(
            "button_release_event", self.figure_mouse_release)
        
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

        self.low_sigma_clip.textChanged.connect(self.check_state)
        self.high_sigma_clip.textChanged.connect(self.check_state)
        self.knot_spacing.textChanged.connect(self.check_state)
        





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
        if state == QtGui.QValidator.Acceptable:
            color = 'none' # normal background color
        elif state == QtGui.QValidator.Intermediate:
            color = '#fff79a' # yellow
        else:
            color = '#f6989d' # red
        sender.setStyleSheet('QLineEdit { background-color: %s }' % color)



    def figure_key_press(self, event):
        """
        Key press event in the normalization figure.
        """

        # Don't allow the user to flip between orders if they are selecting a
        # region to exclude.
        if hasattr(self, "_waiting_on_e_key") and event.key not in "eE":
            return None

        # Show a new order.
        if event.key in ("left", "right"):
            offset = 1 if event.key == "right" else -1

            # TODO: deal with discarded order indices, etc.

            self.update_order_index(np.clip(self.current_order_index + offset,
                0, len(self.parent.session.input_spectra) - 1))
            self.draw_order()
            self.fit_continuum(False)
            self.draw_continuum(True)

            # Do all of the orders have continuum? If so, update the button.
            if not any([(_ is None) for _ in \
            self.parent.session.metadata["normalization"]["continuum"]]):
                self.stitch_btn.setText("Stitch orders")

            return None

        # Scale the continuum up/down.
        if event.key in ("up", "down"):
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


    def figure_mouse_press(self, event):
        """
        Function to handle event clicks (single or double click).
        
        :param event:
            The matplotlib event signal.
        """
        
        # Add/remove an additional point?
        if event.dblclick:

            if event.button == 1:
                # Add a point.
                points = np.vstack([
                    self.ax_order.collections[0].get_offsets(),
                    [event.xdata, event.ydata]
                ])
                # TODO: set size by their weight?
                self.ax_order.collections[0].set_offsets(points)

            else:
                # Are we within <tolerance of a point?
                points = self.ax_order.collections[0].get_offsets()

                # Need to scale x-distance to convert to pixels.
                idx = self.current_order.dispersion.searchsorted(event.xdata)
                xscale = np.nanmean(
                    np.diff(self.current_order.dispersion[idx-5:idx+5]))

                bbox = self.ax_order.get_window_extent().transformed(
                    self.norm_plot.dpi_scale_trans.inverted())
                width =  bbox.width * self.norm_plot.dpi
                height = bbox.height * self.norm_plot.dpi

                print(width, height)
                # TODO: Fix this distance thing.

                distance = np.sqrt(
                      ((points[:, 0] - event.xdata)/xscale)**2 \
                    + (points[:, 1] - event.ydata)**2)
                
                index = np.argmin(distance)
                if distance[index] < PIXEL_PICKER_TOLERANCE:
                    # Remove that point.
                    keep = np.ones(points.shape[0], dtype=bool)
                    keep[index] = False
                    self.ax_order.collections[0].set_offsets(points[keep])

                else:
                    print("Closest point {} px away".format(distance[index]))

            # Update the cache.
            nm = self.parent.session.metadata["normalization"]
            idx = self.current_order_index
            N = points.shape[0]
            nm["normalization_kwargs"][idx]["additional_points"] \
                = np.hstack((points, 200 * np.ones(N).reshape((N, 1))))
            self.fit_continuum(clobber=True)
            self.draw_continuum(refresh=True)

            return None

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
        
        # If the two events were within <1 second, then we should not add a
        # mask.
        signal_time, signal_cid = self._exclude_selected_region_signal
        data = self._exclude_selected_region.get_xy()
        
        if time() - signal_time < DOUBLE_CLICK_INTERVAL: 
            # Don't do anything because it was probably a double click.
            data[:, 0] = np.nan
            
        elif np.abs(data[0,0] - event.xdata) > 0:
            print("ok do something!")
            data[2:4, 0] = event.xdata
            print("ADDED MASK of size: ", np.abs(data[0, 0] - event.xdata))



        self._exclude_selected_region.set_xy(data)
        self.norm_plot.mpl_disconnect(signal_cid)
        self.norm_plot.draw()
        


    def update_exclude_selected_region(self, event):
        """
        Update the selected exclusion region for this order.

        :param event:
            The matplotlib motion event to show the current mouse position.
        """
        
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
            "knot_spacing", "max_iterations")
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
            self._cache["masks"].keys().index(
                self._cache["default_mask"]))



        # Draw the widgets.
        self.update_order_index(0)
        self.update_continuum_mask()
        self.fit_continuum(False)
        self.draw_order()
        self.draw_continuum(True)
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
        self.update_continuum_mask()
        self.fit_continuum(clobber=True)
        self.draw_order()
        self.draw_continuum(refresh=True)

        return None


    def update_continuum_mask(self):
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

        # Transformation function.
        transform = lambda start, end, v=0: np.array([
                [start * (1 - v/c), ymin],
                [start * (1 - v/c), ymax],
                [end   * (1 - v/c), ymax],
                [end   * (1 - v/c), ymin],
                [start * (1 - v/c), ymin]
            ])

        mask = self._cache["masks"][self.continuum_mask.currentText()]

        # Get the applied velocity to shift some masks.
        try:
            v = self.parent.session.metadata["rv"]["rv_applied"]
        except (AttributeError, KeyError):
            v = 0

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
        for i, (start, end) in enumerate(mask.get("rest_wavelength", [])):
            if i >= len(self._masked_wavelengths):
                # Create a polygon in the main figure.
                self._masked_wavelengths.append(self.ax_order.axvspan(**kwds))

                # And for the normalization preview.
                self._masked_wavelengths_norm.append(
                    self.ax_order_norm.axvspan(**kwds))

            polygons = \
                (self._masked_wavelengths[i], self._masked_wavelengths_norm[i])
            for polygon in polygons:
                polygon.set_xy(transform(start, end))

        # obs_wavelength
        for j, (start, end) in enumerate(mask.get("obs_wavelength", [])):
            if i + j + 1 >= len(self._masked_wavelengths):
                # Create a polygon in the main figure.
                self._masked_wavelengths.append(self.ax_order.axvspan(**kwds))

                # And for the normalization preview.
                self._masked_wavelengths_norm.append(
                    self.ax_order_norm.axvspan(**kwds))

            for polygon in (
                    self._masked_wavelengths[i + j + 1], 
                    self._masked_wavelengths_norm[i + j + 1]):
                polygon.set_xy(transform(start, end, v))

        # Any leftover polygons?
        for polygon in self._masked_wavelengths[i + j + 2:]:
            polygon.set_xy(transform(np.nan, np.nan))

        for polygon in self._masked_wavelengths_norm[i + j + 2:]:
            polygon.set_xy(transform(np.nan, np.nan))

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
        

    def update_high_sigma_clip(self):
        """ Update the high sigma clip value. """
        high_sigma = self.high_sigma_clip.text()
        if high_sigma:
            self._cache["input"]["high_sigma_clip"] = float(high_sigma)
            self.reset_input_style_defaults()
            self.fit_continuum(True)
            self.draw_continuum(True)
        return None


    def update_low_sigma_clip(self):
        """ Update the low sigma clip value. """
        low_sigma = self.low_sigma_clip.text()
        if low_sigma:
            self._cache["input"]["low_sigma_clip"] = float(low_sigma)
            self.reset_input_style_defaults()
            self.fit_continuum(True)
            self.draw_continuum(True)
        return None


    def update_normalization_function(self):
        """ Update the normalization function. """
        self._cache["input"]["function"] = self.function.currentText()
        self.reset_input_style_defaults()
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
            index = self.current_order_index

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


    def check_for_different_input_settings(self):
        """
        Check whether the current input settings reflect the settings used to
        normalize the currently displayed order.
        """

        session, index = self.parent.session, self.current_order_index

        # Is there continuum already for this new order?
        continuum = session.metadata["normalization"]["continuum"][index]
        normalization_kwargs \
            = session.metadata["normalization"]["normalization_kwargs"][index]

        if continuum is None: return

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
        }

        diff = dict_updated(self._cache["input"], normalization_kwargs)
        if continuum is not None and diff:
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
        else:
            # Ensure all the things are styled normally.
            self.reset_input_style_defaults(sum(input_items.values(), []))

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
        )
        # Ensure all the things are styled normally.
        for item in items:
            item.setStyleSheet('{0} {{ font-weight: normal }}'\
                .format(item.__class__.__name__))
            item.setStatusTip("")
        return None


    def draw_order(self, refresh=False):
        """
        Draw the current order.

        Note: The order in `self.current_order` is already in the rest frame. 
        """

        x, y = (self.current_order.dispersion, self.current_order.flux)
        
        self.ax_order.lines[0].set_data([x, y])
        self.ax_order.set_xlim(x[0], x[-1])
        self.ax_order.set_ylim(np.nanmin(y), np.nanmax(y))

        self.ax_order.set_title("Order {0} of {1}".format(
            1 + self.current_order_index, 
            len(self.parent.session.input_spectra)))

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

        continuum = session.metadata["normalization"]["continuum"][index]
        if continuum is not None and not clobber:
            # Nothing to do.
            return

        kwds = self._cache["input"].copy()
        kwds["full_output"] = True

        # Add in any additonal points/masked region to kwds.
        existing = \
            session.metadata["normalization"]["normalization_kwargs"][index]
        for key in ("additional_points", ):
            kwds[key] = existing.get(key, None)

        try:
            normalized_spectrum, continuum, left, right \
                = self.current_order.fit_continuum(**kwds)
            if kwds["knot_spacing"] == 201:
                raise KeyError("what") #HACK #TESTING #TODO
        except:
            logger.exception("No continuum could be fit.")
            self.parent.statusbar.showMessage(
                "Exception occurred while trying to fit the continuum.")

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
        continuum = meta["continuum"][index]

        self.ax_order.lines[1].set_data([
            self.current_order.dispersion, continuum])

        # Set the color of the line.
        self.ax_order.lines[1].set_color(wavelength_to_hex(
            np.nanmean(self.current_order.dispersion)))

        # Update the normalization preview in the lower axis.
        self.ax_order_norm.lines[1].set_data([
            self.current_order.dispersion, 
            self.current_order.flux/continuum])
        self.ax_order_norm.set_xlim(self.ax_order.get_xlim())

        # Draw the additional points.
        ap = meta["normalization_kwargs"][index].get("additional_points", None)
        if ap is None:
            ap = np.array([])
        else:
            ap = ap[:, :2]
        self.ax_order.collections[0].set_offsets(ap)

        if refresh:
            self.norm_plot.draw()
        return None

