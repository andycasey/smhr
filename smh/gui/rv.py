#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The radial velocity tab view for the Spectroscopy Made Hard GUI. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import sys
import yaml
from PySide import QtCore, QtGui

import mpl
from matplotlib import (gridspec, pyplot as plt)

from smh import (Session, specutils)

__all__ = ["initialise_tab"]


class RVTab(QtGui.QWidget):


    def __init__(self, parent=None):
        super(RVTab, self).__init__(parent)

        self.parent = parent

        # Create the overall RV tab.
        rv_tab = self #QtGui.QWidget()
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)

        # Create a top-level horizontal layout to contain a matplotlib figure and
        # a vertical layout of settings..
        rv_tab_layout = QtGui.QHBoxLayout(self)

        # This vertical layout will be for input settings.
        rv_settings_vbox = QtGui.QVBoxLayout()
        rv_settings_vbox.setSizeConstraint(QtGui.QLayout.SetMinimumSize)

        # A two-tab setting for 'template' and 'normalization' settings in the
        # radial velocity determination.
        rv_settings_tabs = QtGui.QTabWidget(self)

        sp = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(rv_settings_tabs.sizePolicy().hasHeightForWidth())

        rv_settings_tabs.setSizePolicy(sp)
        rv_settings_tabs.setMinimumSize(QtCore.QSize(300, 240))
        rv_settings_tabs.setMaximumSize(QtCore.QSize(300, 240))

        rv_settings_tabs.setMovable(True)
        rv_settings_tabs.setObjectName("rv_settings_tabs")


        # A tab containing template settings.
        template_tab = QtGui.QWidget()
        template_tab.setObjectName("tv_template_tab")

        # TODO: Should template_tab and template_tab_widget actually be one entity?
        template_tab_widget = QtGui.QWidget(template_tab)
        template_tab_widget.setGeometry(QtCore.QRect(0, 0, 300, 210))

        template_tab_layout = QtGui.QVBoxLayout(template_tab_widget)
        template_tab_layout.setSizeConstraint(QtGui.QLayout.SetDefaultConstraint)
        template_tab_layout.setContentsMargins(10, 10, 10, 10)

        # Add a label at the top of the template settings tab.
        label = QtGui.QLabel(template_tab_widget)
        label.setMaximumSize(QtCore.QSize(240, 16777215))
        label.setText("Template spectrum filename:")
        template_tab_layout.addWidget(label)

        # Inner horizontal layout for the template path and select.
        template_path_layout = QtGui.QHBoxLayout()

        # Template path line edit (read-only).
        self.template_path = QtGui.QLineEdit(template_tab_widget)
        self.template_path.setObjectName("rv_template_path")
        self.template_path.setReadOnly(True)
        template_path_layout.addWidget(self.template_path)

        # Template path select button.
        rv_select_template_btn = QtGui.QPushButton(template_tab_widget)
        rv_select_template_btn.setObjectName("rv_select_template_btn")
        rv_select_template_btn.setText("...")
        template_path_layout.addWidget(rv_select_template_btn)

        # Add this horizontal layout to the template tab.
        template_tab_layout.addLayout(template_path_layout)

        # Add a label for the wavelength regions
        label = QtGui.QLabel(template_tab_widget)
        label.setMaximumSize(QtCore.QSize(240, 16777215))
        label.setText("Wavelength region:")
        template_tab_layout.addWidget(label)

        # Wavelength region for CCF.
        wl_region_layout = QtGui.QHBoxLayout()
        self.wl_region = QtGui.QComboBox(template_tab_widget)
        self.wl_region.setSizeAdjustPolicy(QtGui.QComboBox.AdjustToContents)
        self.wl_region.setObjectName("rv_wl_region")
        wl_region_layout.addWidget(self.wl_region)

        # Edit button for the wavelength regions.
        rv_wl_region_edit = QtGui.QPushButton(template_tab_widget)
        rv_wl_region_edit.setMaximumSize(QtCore.QSize(80, 16777215))
        rv_wl_region_edit.setObjectName("rv_wl_region_edit")
        wl_region_layout.addWidget(rv_wl_region_edit)
        template_tab_layout.addLayout(wl_region_layout)
        rv_wl_region_edit.setText("Edit list")
        
        # Add a horizontal line.
        hr = QtGui.QFrame(template_tab_widget)
        hr.setFrameShape(QtGui.QFrame.HLine)
        hr.setFrameShadow(QtGui.QFrame.Sunken)
        template_tab_layout.addWidget(hr)

        # Add a flexible spacer.
        template_tab_layout.addItem(QtGui.QSpacerItem(
            20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding))

        # Add a cross-correlate button.
        rv_cross_correlate_btn = QtGui.QPushButton(template_tab_widget)
        rv_cross_correlate_btn.setObjectName("rv_cross_correlate_btn")
        rv_cross_correlate_btn.setText("Cross-correlate")
        template_tab_layout.addWidget(rv_cross_correlate_btn)

        # End of the template tab.
        rv_settings_tabs.addTab(template_tab, "Template")

        # Start the normalization tab.
        norm_tab = QtGui.QWidget()
        norm_tab.setObjectName("rv_normalization_tab")

        norm_tab_widget = QtGui.QWidget(norm_tab)
        norm_tab_widget.setGeometry(QtCore.QRect(0, 0, 300, 210))

        norm_tab_layout = QtGui.QVBoxLayout(norm_tab_widget)
        norm_tab_layout.setContentsMargins(10, 10, 10, 10)

        # Start the grid layout for the normalization tab.
        norm_tab_grid_layout = QtGui.QGridLayout()


        # Normalization function.
        label = QtGui.QLabel(norm_tab_widget)
        label.setText("Function")
        norm_tab_grid_layout.addWidget(label, 0, 0, 1, 1)
        
        # Put the normalization function combo box in a horizontal layout with a
        # spacer.
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_function = QtGui.QComboBox(norm_tab_widget)
        self.norm_function.setObjectName("rv_norm_function")
        hbox.addWidget(self.norm_function)
        norm_tab_grid_layout.addLayout(hbox, 0, 1, 1, 1)

        norm_functions = ("polynomial", "spline")
        for each in norm_functions:
            self.norm_function.addItem(each.title())


        # Normalization function order.
        label = QtGui.QLabel(norm_tab_widget)
        label.setText("Order")
        norm_tab_grid_layout.addWidget(label, 1, 0, 1, 1)
        
        # Put the normalization order combo box in a horizontal layout with a spacer
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_order = QtGui.QComboBox(norm_tab_widget)
        self.norm_order.setMaximumSize(QtCore.QSize(50, 16777215))
        self.norm_order.setObjectName("rv_norm_order")
        hbox.addWidget(self.norm_order)
        norm_tab_grid_layout.addLayout(hbox, 1, 1, 1, 1)

        norm_orders = range(1, 10)
        for order in norm_orders:
            self.norm_order.addItem("{0:.0f}".format(order))


        # Maximum number of iterations.
        label = QtGui.QLabel(norm_tab_widget)
        label.setText("Maximum iterations")
        norm_tab_grid_layout.addWidget(label, 2, 0, 1, 1)

        # Put the maxium number of iterations in a horizontal layout with a spacer.
        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_max_iter = QtGui.QComboBox(norm_tab_widget)
        self.norm_max_iter.setMaximumSize(QtCore.QSize(50, 16777215))
        self.norm_max_iter.setObjectName("rv_norm_max_iter")
        hbox.addWidget(self.norm_max_iter)
        norm_tab_grid_layout.addLayout(hbox, 2, 1, 1, 1)

        norm_max_iters = range(1, 10)
        for iteration in norm_max_iters:
            self.norm_max_iter.addItem("{0:.0f}".format(iteration))


        # Low sigma clipping.
        label = QtGui.QLabel(norm_tab_widget)
        label.setText("Low sigma clip")
        norm_tab_grid_layout.addWidget(label, 3, 0, 1, 1)

        # Put the low sigma line edit box in a horizontal layout with a spacer.
        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(-1, -1, 5, -1)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_low_sigma = QtGui.QLineEdit(norm_tab_widget)
        self.norm_low_sigma.setMaximumSize(QtCore.QSize(40, 16777215))
        self.norm_low_sigma.setAlignment(QtCore.Qt.AlignCenter)
        self.norm_low_sigma.setObjectName("rv_norm_low_sigma")
        hbox.addWidget(self.norm_low_sigma)
        norm_tab_grid_layout.addLayout(hbox, 3, 1, 1, 1)


        # High sigma clipping.
        label = QtGui.QLabel(norm_tab_widget)
        label.setText("High sigma clip")
        norm_tab_grid_layout.addWidget(label, 4, 0, 1, 1)

        # Put the high sigma line edit box in a horizontal layout with a spacer.
        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(-1, -1, 5, -1)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_high_sigma = QtGui.QLineEdit(norm_tab_widget)
        self.norm_high_sigma.setMaximumSize(QtCore.QSize(40, 16777215))
        self.norm_high_sigma.setAlignment(QtCore.Qt.AlignCenter)
        self.norm_high_sigma.setObjectName("rv_norm_high_sigma")
        hbox.addWidget(self.norm_high_sigma)
        norm_tab_grid_layout.addLayout(hbox, 4, 1, 1, 1)
        

        # Knot spacing.
        label = QtGui.QLabel(norm_tab_widget)
        norm_tab_grid_layout.addWidget(label, 5, 0, 1, 1)
        label.setText("Knot spacing (A)")
        # TODO: Put unicode Angstroms character in

        # Put the knot spacing lint edit box in a horizontal layout with a spacer
        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(-1, -1, 5, -1)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self.norm_knot_spacing = QtGui.QLineEdit(norm_tab_widget)
        self.norm_knot_spacing.setMaximumSize(QtCore.QSize(40, 16777215))
        self.norm_knot_spacing.setAlignment(QtCore.Qt.AlignCenter)
        self.norm_knot_spacing.setObjectName("rv_norm_knot_spacing")
        hbox.addWidget(self.norm_knot_spacing)
        norm_tab_grid_layout.addLayout(hbox, 5, 1, 1, 1)

        # End of the grid in the normalization tab.
        norm_tab_layout.addLayout(norm_tab_grid_layout)
        norm_tab_layout.addItem(QtGui.QSpacerItem(
            20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding))

        rv_settings_tabs.addTab(norm_tab, "Order normalization")
        rv_settings_vbox.addWidget(rv_settings_tabs)


        # Horizontal layout for the radial velocity measured/corrected.
        hbox = QtGui.QHBoxLayout()
        hbox.setSizeConstraint(QtGui.QLayout.SetMaximumSize)
        hbox.setContentsMargins(10, 0, 10, -1)

        label = QtGui.QLabel(self)
        label.setText("Radial velocity:")
        hbox.addWidget(label)

        # Radial velocity measured.
        self.rv_applied = QtGui.QLineEdit(self)
        sp = QtGui.QSizePolicy(QtGui.QSizePolicy.Ignored, QtGui.QSizePolicy.Fixed)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(self.rv_applied.sizePolicy().hasHeightForWidth())
        self.rv_applied.setSizePolicy(sp)
        self.rv_applied.setMinimumSize(QtCore.QSize(50, 16777215))
        self.rv_applied.setAlignment(QtCore.Qt.AlignCenter)
        self.rv_applied.setObjectName("rv_applied")
        hbox.addWidget(self.rv_applied)

        # Units/uncertainty label.
        label = QtGui.QLabel(self)
        label.setObjectName("rv_measured_units_label")
        label.setText("km/s")
        hbox.addWidget(label)

        # Correct for radial velocity button.
        rv_correct_btn = QtGui.QPushButton(self)
        rv_correct_btn.setObjectName("rv_correct_btn")
        rv_correct_btn.setText("Correct")
        hbox.addWidget(rv_correct_btn)
        rv_settings_vbox.addLayout(hbox)


        # Add a spacer until the big button.
        rv_settings_vbox.addItem(QtGui.QSpacerItem(
            20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding))


        # The cross-correlate and correct button.
        rv_ccc_btn = QtGui.QPushButton(self)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(rv_ccc_btn.sizePolicy().hasHeightForWidth())
        rv_ccc_btn.setSizePolicy(sp)
        rv_ccc_btn.setMinimumSize(QtCore.QSize(300, 0))
        rv_ccc_btn.setMaximumSize(QtCore.QSize(300, 16777215))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        rv_ccc_btn.setFont(font)
        rv_ccc_btn.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        rv_ccc_btn.setDefault(True)
        rv_ccc_btn.setObjectName("rv_ccc_btn")
        rv_ccc_btn.setText("Cross-correlate and correct")
        if sys.platform == "darwin":
            rv_ccc_btn.setStyleSheet('QPushButton {color: white}')

        rv_settings_vbox.addWidget(rv_ccc_btn)

        rv_tab_layout.addLayout(rv_settings_vbox)

        # Create a matplotlib widget.
        blank_widget = QtGui.QWidget(self)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(blank_widget.sizePolicy().hasHeightForWidth())
        blank_widget.setSizePolicy(sp)
        blank_widget.setObjectName("blank_widget")

        self.rv_plot = mpl.MPLWidget(blank_widget)
        layout = QtGui.QVBoxLayout(blank_widget)        
        layout.addWidget(self.rv_plot, 1)

        rv_tab_layout.addWidget(blank_widget)


        gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 1])
        self.ax_order = self.rv_plot.figure.add_subplot(gs[0])
        self.ax_order_norm = self.rv_plot.figure.add_subplot(gs[1])
        self.ax_ccf = self.rv_plot.figure.add_subplot(gs[2])

        # Ticks, etc
        self.ax_order.set_xticklabels([])
        self.ax_order_norm.set_yticks([0, 1])
        self.ax_order_norm.set_ylim(0, 1.2)

        self.ax_order_norm.set_xlabel(u"Wavelength (Å)")
        self.ax_order.set_ylabel("Flux")

        # Draw an initial line for data and continuum.
        self.ax_order.plot([], [], c='k', drawstyle='steps-mid')
        self.ax_order.plot([], [], c='r', zorder=2)


        self.ax_order_norm.axhline(1, linestyle=":", c="#666666", zorder=-1)
        self.ax_order_norm.plot([], [], c='k', drawstyle='steps-mid')
        self.ax_order_norm.plot([], [], c='b') # Template.


        self.ax_ccf.plot([0, 1], [0, 1], 'ro-')

        self.rv_plot.draw()

        

        # Add the tab to the application.


        # Create signals
        rv_select_template_btn.clicked.connect(self.select_template) 
        rv_cross_correlate_btn.clicked.connect(self.cross_correlate) 
        rv_correct_btn.clicked.connect(self.correct_radial_velocity)
        rv_ccc_btn.clicked.connect(self.cross_correlate_and_correct)

        self.wl_region.currentIndexChanged.connect(self.redraw_order)

        # Read in the default settings from the SMH session file.

        # TODO: Put the default I/O somewhere else since it will be common to many
        #       tabs.
        with open(Session._default_settings_path, "rb") as fp:
            defaults = yaml.load(fp)["rv"]

        # Template filename.
        self.template_path.setReadOnly(False)
        self.template_path.setText(defaults["template_path"])
        self.template_path.setReadOnly(True)

        # Wavelength regions.
        for each in defaults["wavelength_regions"]:
            self.wl_region.addItem(u"{0:.0f}-{1:.0f} Å".format(*each))

        # Normalization function.
        self.norm_function.setCurrentIndex(
            norm_functions.index(defaults["normalization"]["function"].lower()))

        # Normalization order.
        self.norm_order.setCurrentIndex(
            norm_orders.index(defaults["normalization"]["order"]))

        # Normalization maximum iterations.
        self.norm_max_iter.setCurrentIndex(
            norm_max_iters.index(defaults["normalization"]["max_iterations"]))

        # Normalization low and high sigma clip:
        low, high = defaults["normalization"]["sigma_clip"]
        self.norm_low_sigma.setText(str(low))
        self.norm_high_sigma.setText(str(high))
        
        # Normalization knot spacing.
        self.norm_knot_spacing.setText(
            str(defaults["normalization"]["knot_spacing"]))

        return None


    def select_template(self):
        """
        Select the template spectrum filename from disk.
        """

        path, _ = QtGui.QFileDialog.getOpenFileName(self.parent,
            caption="Select template spectrum", dir="", filter="*.fits")
        if not path: return

        # Set the text.
        self.template_path.setReadOnly(False)
        self.template_path.setText(path)
        self.template_path.setReadOnly(True)
        
        return None


    def cross_correlate(self):
        """
        Normalize and cross-correlate the observed spectrum with the template.
        """

        # Get all the inputs.
        kwds = {
            "template_spectrum": self.template_path.text(),
            "wavelength_region": [float(_) \
                for _ in self.wl_region.currentText().split(" ")[0].split("-")],
            "resample": "template",
            "apodize": 0,
            "normalization_kwargs": {
                "function": self.norm_function.currentText(),
                "order": int(self.norm_order.currentText()),
                "max_iterations": int(self.norm_max_iter.currentText()),
                "sigma_clip": [
                    float(self.norm_low_sigma.text()),
                    float(self.norm_high_sigma.text())
                ],
                "knot_spacing": float(self.norm_knot_spacing.text())
            }
        }

        # Perform the cross-correlation.
        rv, rv_uncertainty = self.parent.session.rv_measure(**kwds)

        # Update the measured radial velocity in the GUI.
        self.rv_applied.setText("{0:+.1f}".format(rv))

        return None


    def correct_radial_velocity(self):
        """
        Correct the radial velocity of the observed spectra.
        """

        self.parent.session.rv["rv_applied"] = float(self.rv_applied.text())

        # Enable the next tab.
        self.parent.tabs.setTabEnabled(self.parent.tabs.indexOf(self) + 1, True)

        return None


    def cross_correlate_and_correct(self):
        """
        Normalize and cross-correlate the observed spectrum with a template,
        then correct the observed spectrum with that velocity.
        """

        self.cross_correlate()
        self.correct_radial_velocity()

        return None


    def redraw_order(self):
        """
        Re-draw the order selected and the continuum fit, as well as the preview
        of the normalized spectrum.
        """

        if self.parent.session is None: return

        print("redrawing order")

        wavelength_region = [float(_) \
            for _ in self.wl_region.currentText().split(" ")[0].split("-")]

        self._overlap_order, overlap_index, _ = \
            self.parent.session._get_overlap_order([wavelength_region]) 

        self.ax_order.lines[0].set_data([
            self._overlap_order.dispersion,
            self._overlap_order.flux,
        ])

        # Limits
        self.ax_order.set_xlim(wavelength_region)

        # Add 10% peak-to-peak.
        flux_limits = (
            np.nanmin(self._overlap_order.flux),
            np.nanmax(self._overlap_order.flux)
        )
        self.ax_order.set_ylim(
            flux_limits[0],
            flux_limits[1] + (np.ptp(flux_limits) * 0.10)
        )

        #self.redraw_continuum(self)

        self.rv_plot.draw()

        print("Assuming that the session has not just been loaded and it has a CCF/norm order, etc")

        # Normalize that order using the normalization settings supplied.
        #observed_spectrum = overlap_order.fit_continuum(**normalization_kwargs)

        self.redraw_continuum()

        return None


    def redraw_continuum(self):

        kwds = {
            "function": self.norm_function.currentText(),
            "order": int(self.norm_order.currentText()),
            "max_iterations": int(self.norm_max_iter.currentText()),
            "sigma_clip": [
                float(self.norm_low_sigma.text()),
                float(self.norm_high_sigma.text())
            ],
            "knot_spacing": float(self.norm_knot_spacing.text()),
            "full_output": True
        }
        normalized_order, continuum, _, __ = \
            self._overlap_order.fit_continuum(**kwds)

        self.ax_order.lines[1].set_data([
            self._overlap_order.dispersion,
            continuum
        ])

        self.ax_order_norm.lines[1].set_data([
            normalized_order.dispersion,
            normalized_order.flux,
        ])

        wavelength_region = [float(_) \
            for _ in self.wl_region.currentText().split(" ")[0].split("-")]
        self.ax_order_norm.set_xlim(wavelength_region)


        self.draw_template()


    def draw_template(self):

        template = specutils.Spectrum1D.read(self.template_path.text())
        self.ax_order_norm.lines[2].set_data([
            template.dispersion,
            template.flux
        ])

        self.rv_plot.draw()



def initialise_tab(tabs, parent):
    """
    Create a radial velocity tab and add it to the application tabs.

    :param tabs:
        The application tab widget to add a radial velocity tab to.

    :type tabs:
        :class:`QtGui.QTabWidget`
    """

    tab = RVTab(parent)

    tabs.addTab(tab, "Radial velocity")
    return tab
