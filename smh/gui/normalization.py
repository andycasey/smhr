#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The normalization tab in Spectroscopy Made Hard """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import sys
from PySide import QtCore, QtGui

from matplotlib import gridspec

import mpl

__all__ = ["NormalizationTab"]


c = 299792458e-3 # km/s


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
        self.norm_function_label = QtGui.QLabel(self)
        self.norm_function_label.setText("Function")
        settings_grid_layout.addWidget(self.norm_function_label, 0, 0, 1, 1)
        
        # Put the normalization function combo box in a horizontal layout with 
        # a spacer.
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_function = QtGui.QComboBox(self)
        self.norm_function.setObjectName("rv_norm_function")
        hbox.addWidget(self.norm_function)
        settings_grid_layout.addLayout(hbox, 0, 1, 1, 1)

        norm_functions = ("polynomial", "spline")
        for each in norm_functions:
            self.norm_function.addItem(each.title())

        # Normalization function order.
        self.norm_order_label = QtGui.QLabel(self)
        self.norm_order_label.setText("Order")
        settings_grid_layout.addWidget(self.norm_order_label, 1, 0, 1, 1)
        
        # Put the normalization order combo box in a horizontal layout with a
        # spacer
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_order = QtGui.QComboBox(self)
        self.norm_order.setMaximumSize(QtCore.QSize(50, 16777215))
        self.norm_order.setObjectName("rv_norm_order")
        hbox.addWidget(self.norm_order)
        settings_grid_layout.addLayout(hbox, 1, 1, 1, 1)

        norm_orders = range(1, 10)
        for order in norm_orders:
            self.norm_order.addItem("{0:.0f}".format(order))

        # Maximum number of iterations.
        self.norm_max_iter_label = QtGui.QLabel(self)
        self.norm_max_iter_label.setText("Maximum iterations")
        settings_grid_layout.addWidget(self.norm_max_iter_label, 2, 0, 1, 1)

        # Put the maxium number of iterations in a horizontal layout with a 
        # spacer.
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_max_iter = QtGui.QComboBox(self)
        self.norm_max_iter.setMaximumSize(QtCore.QSize(50, 16777215))
        self.norm_max_iter.setObjectName("rv_norm_max_iter")
        hbox.addWidget(self.norm_max_iter)
        settings_grid_layout.addLayout(hbox, 2, 1, 1, 1)

        norm_max_iters = range(1, 10)
        for iteration in norm_max_iters:
            self.norm_max_iter.addItem("{0:.0f}".format(iteration))


        # Low sigma clipping.
        label = QtGui.QLabel(self)
        label.setText("Low sigma clip")
        settings_grid_layout.addWidget(label, 3, 0, 1, 1)

        # Put the low sigma line edit box in a horizontal layout with a spacer.
        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(-1, -1, 5, -1)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_low_sigma = QtGui.QLineEdit(self)
        self.norm_low_sigma.setMaximumSize(QtCore.QSize(40, 16777215))
        self.norm_low_sigma.setAlignment(QtCore.Qt.AlignCenter)
        self.norm_low_sigma.setObjectName("rv_norm_low_sigma")
        self.norm_low_sigma.setValidator(
            QtGui.QDoubleValidator(0, 1000, 2, self.norm_low_sigma))

        hbox.addWidget(self.norm_low_sigma)
        settings_grid_layout.addLayout(hbox, 3, 1, 1, 1)


        # High sigma clipping.
        label = QtGui.QLabel(self)
        label.setText("High sigma clip")
        settings_grid_layout.addWidget(label, 4, 0, 1, 1)

        # Put the high sigma line edit box in a horizontal layout with a spacer.
        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(-1, -1, 5, -1)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_high_sigma = QtGui.QLineEdit(self)
        self.norm_high_sigma.setMaximumSize(QtCore.QSize(40, 16777215))
        self.norm_high_sigma.setAlignment(QtCore.Qt.AlignCenter)
        self.norm_high_sigma.setObjectName("rv_norm_high_sigma")
        self.norm_high_sigma.setValidator(
            QtGui.QDoubleValidator(0, 1000, 2, self.norm_high_sigma))
        hbox.addWidget(self.norm_high_sigma)
        settings_grid_layout.addLayout(hbox, 4, 1, 1, 1)
        

        # Knot spacing.
        self.norm_knot_spacing_label = QtGui.QLabel(self)
        settings_grid_layout.addWidget(self.norm_knot_spacing_label, 5, 0, 1, 1)
        self.norm_knot_spacing_label.setText(u"Knot spacing (Å)")

        # Put the knot spacing lint edit box in a horizontal layout with a spacer
        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(-1, -1, 5, -1)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_knot_spacing = QtGui.QLineEdit(self)
        self.norm_knot_spacing.setMaximumSize(QtCore.QSize(40, 16777215))
        self.norm_knot_spacing.setAlignment(QtCore.Qt.AlignCenter)
        self.norm_knot_spacing.setValidator(
            QtGui.QIntValidator(0, 10000, self.norm_knot_spacing))
        self.norm_knot_spacing.setObjectName("rv_norm_knot_spacing")
        hbox.addWidget(self.norm_knot_spacing)
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
        self.ax_order.plot([], [], c='r', zorder=4)

        # Line for the neighbouring order(s) (joined by a NaN).
        self.ax_order.plot([], [], c='#666666', zorder=1, drawstyle='steps-mid')
        # Line for the neighbouring order(s) continuum (joined by a NaN)
        self.ax_order.plot([], [], c='b', zorder=2)

        # Additional point markers.
        # TODO
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
        self.norm_plot.canvas.mpl_connect(
            "key_press_event", self.figure_key_press)


        

        self.norm_function.currentIndexChanged.connect(
            self.update_normalization_function)
        self.norm_order.currentIndexChanged.connect(
            self.update_normalization_order)
        self.norm_max_iter.currentIndexChanged.connect(
            self.update_normalization_max_iterations)
        self.norm_low_sigma.textChanged.connect(
            self.update_low_sigma_clip)
        self.norm_high_sigma.textChanged.connect(
            self.update_high_sigma_clip)
        self.norm_knot_spacing.textChanged.connect(self.update_knot_spacing)

        self.norm_low_sigma.textChanged.connect(self.check_state)
        self.norm_high_sigma.textChanged.connect(self.check_state)
        self.norm_knot_spacing.textChanged.connect(self.check_state)
        



        return None


    def check_state(self, *args, **kwargs):
        """
        Update the background color of a QLineEdit object based on whether the
        input is valid.
        """

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

        print(event, event.key)

        if event.key in ("left", "right"):
            offset = 1 if event.key == "right" else -1

            self.update_order_index(np.clip(self.current_order_index + offset,
                0, len(self.parent.session.input_spectra) - 1))
            self.draw_order()
            self.fit_continuum(False)
            self.draw_continuum(True)

            return None

        # TODO: deal with discarded order indices, etc.


    def _populate_widgets(self):
        """
        Populate the widgets with default settings, and a default figure view.
        """

        if self.parent.session is None:
            # No point populating the widgets with the default values from the
            # SMH file because these will be updated when a session is loaded.
            return

        keys = ("function", "order", "sigma_clip", "knot_spacing",
            "max_iterations")
        self._cache = {
            "input": {}
        }
        for key in keys:
            self._cache["input"][key] = self.parent.session._default(
                None, ("normalization", key))

        # Put these values into the widgets.
        self.norm_low_sigma.setText(str(self._cache["input"]["sigma_clip"][0]))
        self.norm_high_sigma.setText(str(self._cache["input"]["sigma_clip"][1]))
        self.norm_knot_spacing.setText(str(
            self._cache["input"]["knot_spacing"]))

        norm_functions = [self.norm_function.itemText(i).lower() \
            for i in range(self.norm_function.count())]
        self.norm_function.setCurrentIndex(norm_functions.index(
            self._cache["input"]["function"]))

        # Normalization order.
        norm_orders = [int(self.norm_order.itemText(i)) \
            for i in range(self.norm_order.count())]
        self.norm_order.setCurrentIndex(norm_orders.index(
            self._cache["input"]["order"]))

        # Normalization maximum iterations.
        norm_max_iters = [int(self.norm_max_iter.itemText(i)) \
            for i in range(self.norm_max_iter.count())]
        self.norm_max_iter.setCurrentIndex(norm_max_iters.index(
            self._cache["input"]["max_iterations"]))

        # Draw the widgets.
        self.update_order_index(0)
        self.fit_continuum(False)
        self.draw_order()
        self.draw_continuum(True)
        return None


    def update_knot_spacing(self):
        """ Update the knot spacing. """
        knot_spacing = self.norm_knot_spacing.text()
        if knot_spacing:
            self._cache["input"]["knot_spacing"] = float(knot_spacing)
            self.fit_continuum(True)
            self.draw_continuum(True)
            self.reset_font_weights()

        return None
        

    def update_high_sigma_clip(self):
        """ Update the high sigma clip value. """
        high_sigma = self.norm_high_sigma.text()
        if high_sigma:
            self._cache["input"]["sigma_clip"][1] = float(high_sigma)
            self.fit_continuum(True)
            self.draw_continuum(True)
            self.reset_font_weights()
        return None


    def update_low_sigma_clip(self):
        """ Update the low sigma clip value. """
        low_sigma = self.norm_low_sigma.text()
        if low_sigma:
            self._cache["input"]["sigma_clip"][0] = float(low_sigma)
            self.fit_continuum(True)
            self.draw_continuum(True)
            self.reset_font_weights()
        return None


    def update_normalization_function(self):
        """ Update the normalization function. """
        self._cache["input"]["function"] = self.norm_function.currentText()
        self.fit_continuum(True)
        self.draw_continuum(True)
        self.reset_font_weights()
        return None


    def update_normalization_order(self):
        """ Update the normalization order. """
        self._cache["input"]["order"] = int(self.norm_order.currentText())
        self.fit_continuum(True)
        self.draw_continuum(True)
        self.reset_font_weights()
        return None


    def update_normalization_max_iterations(self):
        """ Update the maximum number of iterations. """
        self._cache["input"]["max_iterations"] \
            = int(self.norm_max_iter.currentText())
        self.fit_continuum(True)
        self.draw_continuum(True)
        self.reset_font_weights()
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
            rv_applied = session.metadata["rv"]["rv_applied"]
        except (AttributeError, KeyError):
            rv_applied = 0

        self.current_order._dispersion *= (1 - rv_applied/c)

        # Update the view if the input settings don't match the settings used
        # to normalize the current order.
        self.check_for_different_input_settings()

        # Do all of the orders have continuum? If so, update the button.
        if not None in session.metadata["normalization"]["continuum"]:
            self.stitch_btn.setText("Stitch orders")

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
            "function": [self.norm_function_label, self.norm_function],
            "order": [self.norm_order_label, self.norm_order],
            "knot_spacing": \
                [self.norm_knot_spacing, self.norm_knot_spacing_label],
            "max_iterations": [self.norm_max_iter_label, self.norm_max_iter],

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
            self.reset_font_weights(sum(input_items.values(), []))

        return None


    def reset_font_weights(self, items=None):

        items = items or (
            self.norm_function_label, self.norm_function,
            self.norm_order_label, self.norm_order,
            self.norm_knot_spacing_label, self.norm_knot_spacing,
            self.norm_max_iter_label, self.norm_max_iter,
        )
        # Ensure all the things are styled normally.
        for item in items:
            item.setStyleSheet('{0} {{ font-weight: normal }}'\
                .format(item.__class__.__name__))
            item.setStatusTip("")


    def draw_order(self, refresh=False):
        """
        Draw the current order.
        """

        x, y = self.current_order.dispersion, self.current_order.flux
        self.ax_order.lines[0].set_data([x, y])
        self.ax_order.set_xlim(x[0], x[-1])
        self.ax_order.set_ylim(np.nanmin(y), np.nanmax(y))

        self.ax_order.set_title("Order {0} of {1}".format(
            1 + self.current_order_index, 
            len(self.parent.session.input_spectra)))

        if refresh:
            self.norm_plot.draw()
        return None


    def fit_continuum(self, clobber):
        """
        Update continuum for the current order.
        """

        # Any existing continuum determination?
        index, session = (self.current_order_index, self.parent.session)

        continuum = session.metadata["normalization"]["continuum"][index]
        if continuum is not None and not clobber:
            # Nothing to do.
            return

        kwds = self._cache["input"].copy()
        kwds["full_output"] = True

        # Add in any additonal points/masked region to kwds.

        normalized_spectrum, continuum, left, right \
            = self.current_order.fit_continuum(**kwds)

        session.metadata["normalization"]["continuum"][index] = continuum
        session.metadata["normalization"]["normalization_kwargs"][index] = kwds

        return None


    def draw_continuum(self, refresh=False):
        """
        Draw the continuum for the current order.
        """

        meta = self.parent.session.metadata["normalization"]
        continuum = meta["continuum"][self.current_order_index]

        self.ax_order.lines[1].set_data([
            self.current_order.dispersion, continuum])

        # Update the normalization preview in the lower axis.
        self.ax_order_norm.lines[1].set_data([
            self.current_order.dispersion, self.current_order.flux/continuum])
        self.ax_order_norm.set_xlim(self.ax_order.get_xlim())

        if refresh:
            self.norm_plot.draw()
        return None

