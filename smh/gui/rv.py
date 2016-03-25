#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The radial velocity tab view for the Spectroscopy Made Hard GUI. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import sys

from PySide import QtCore, QtGui

import mpl

__all__ = ["initialise_tab"]



def initialise_tab(tabs):
    """
    Create a radial velocity tab and add it to the application tabs.

    :param tabs:
        The application tab widget to add a radial velocity tab to.

    :type tabs:
        :class:`QtGui.QTabWidget`
    """


    # Create the overall RV tab.
    rv_tab = QtGui.QWidget()
    sp = QtGui.QSizePolicy(
        QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding)
    sp.setHeightForWidth(rv_tab.sizePolicy().hasHeightForWidth())
    rv_tab.setSizePolicy(sp)

    # Create a top-level horizontal layout to contain a matplotlib figure and
    # a vertical layout of settings..
    rv_tab_layout = QtGui.QHBoxLayout(rv_tab)

    # This vertical layout will be for input settings.
    rv_settings_vbox = QtGui.QVBoxLayout()
    rv_settings_vbox.setSizeConstraint(QtGui.QLayout.SetMinimumSize)

    # A two-tab setting for 'template' and 'normalization' settings in the
    # radial velocity determination.
    rv_settings_tabs = QtGui.QTabWidget(rv_tab)

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
    rv_template_path_linedit = QtGui.QLineEdit(template_tab_widget)
    rv_template_path_linedit.setObjectName("rv_template_path_linedit")
    rv_template_path_linedit.setReadOnly(True)
    template_path_layout.addWidget(rv_template_path_linedit)

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
    rv_wl_region = QtGui.QComboBox(template_tab_widget)
    rv_wl_region.setSizeAdjustPolicy(QtGui.QComboBox.AdjustToContents)
    rv_wl_region.setObjectName("rv_wl_region")
    wl_region_layout.addWidget(rv_wl_region)

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
    norm_tab.setObjectName("normalization_tab")

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
    rv_norm_function = QtGui.QComboBox(norm_tab_widget)
    rv_norm_function.setObjectName("rv_norm_function")
    hbox.addWidget(rv_norm_function)
    norm_tab_grid_layout.addLayout(hbox, 0, 1, 1, 1)

    rv_norm_function.addItem("Polynomial")
    rv_norm_function.addItem("Spline")


    # Normalization function order.
    label = QtGui.QLabel(norm_tab_widget)
    label.setText("Order")
    norm_tab_grid_layout.addWidget(label, 1, 0, 1, 1)
    
    # Put the normalization order combo box in a horizontal layout with a spacer
    hbox = QtGui.QHBoxLayout()
    hbox.addItem(QtGui.QSpacerItem(
        40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
    rv_norm_order = QtGui.QComboBox(norm_tab_widget)
    rv_norm_order.setMaximumSize(QtCore.QSize(50, 16777215))
    rv_norm_order.setObjectName("rv_norm_order")
    hbox.addWidget(rv_norm_order)
    norm_tab_grid_layout.addLayout(hbox, 1, 1, 1, 1)

    for order in range(1, 10):
        rv_norm_order.addItem("{0:.0f}".format(order))


    # Maximum number of iterations.
    label = QtGui.QLabel(norm_tab_widget)
    label.setText("Maximum iterations")
    norm_tab_grid_layout.addWidget(label, 2, 0, 1, 1)

    # Put the maxium number of iterations in a horizontal layout with a spacer.
    hbox = QtGui.QHBoxLayout()
    hbox.addItem(QtGui.QSpacerItem(
        40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
    rv_norm_max_iter = QtGui.QComboBox(norm_tab_widget)
    rv_norm_max_iter.setMaximumSize(QtCore.QSize(50, 16777215))
    rv_norm_max_iter.setObjectName("rv_norm_max_iter")
    hbox.addWidget(rv_norm_max_iter)
    norm_tab_grid_layout.addLayout(hbox, 2, 1, 1, 1)

    for iteration in range(1, 10):
        rv_norm_max_iter.addItem("{0:.0f}".format(iteration))


    # Low sigma clipping.
    label = QtGui.QLabel(norm_tab_widget)
    label.setText("Low sigma clip")
    norm_tab_grid_layout.addWidget(label, 3, 0, 1, 1)

    # Put the low sigma line edit box in a horizontal layout with a spacer.
    hbox = QtGui.QHBoxLayout()
    hbox.setContentsMargins(-1, -1, 5, -1)
    hbox.addItem(QtGui.QSpacerItem(
        40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
    rv_norm_low_sigma = QtGui.QLineEdit(norm_tab_widget)
    rv_norm_low_sigma.setMaximumSize(QtCore.QSize(40, 16777215))
    rv_norm_low_sigma.setAlignment(QtCore.Qt.AlignCenter)
    rv_norm_low_sigma.setObjectName("rv_norm_low_sigma")
    hbox.addWidget(rv_norm_low_sigma)
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
    rv_norm_high_sigma = QtGui.QLineEdit(norm_tab_widget)
    rv_norm_high_sigma.setMaximumSize(QtCore.QSize(40, 16777215))
    rv_norm_high_sigma.setAlignment(QtCore.Qt.AlignCenter)
    rv_norm_high_sigma.setObjectName("rv_norm_high_sigma")
    hbox.addWidget(rv_norm_high_sigma)
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
    rv_norm_knot_spacing = QtGui.QLineEdit(norm_tab_widget)
    rv_norm_knot_spacing.setMaximumSize(QtCore.QSize(40, 16777215))
    rv_norm_knot_spacing.setAlignment(QtCore.Qt.AlignCenter)
    rv_norm_knot_spacing.setObjectName("rv_norm_knot_spacing")
    hbox.addWidget(rv_norm_knot_spacing)
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

    label = QtGui.QLabel(rv_tab)
    label.setText("Radial velocity:")
    hbox.addWidget(label)

    # Radial velocity measured.
    rv_measured = QtGui.QLineEdit(rv_tab)
    sp = QtGui.QSizePolicy(QtGui.QSizePolicy.Ignored, QtGui.QSizePolicy.Fixed)
    sp.setHorizontalStretch(0)
    sp.setVerticalStretch(0)
    sp.setHeightForWidth(rv_measured.sizePolicy().hasHeightForWidth())
    rv_measured.setSizePolicy(sp)
    rv_measured.setMaximumSize(QtCore.QSize(60, 16777215))
    rv_measured.setAlignment(QtCore.Qt.AlignCenter)
    rv_measured.setObjectName("rv_measured")
    hbox.addWidget(rv_measured)

    # Units/uncertainty label.
    rv_measured_units_label = QtGui.QLabel(rv_tab)
    rv_measured_units_label.setObjectName("rv_measured_units_label")
    rv_measured_units_label.setText("km/s")
    hbox.addWidget(rv_measured_units_label)

    # Correct for radial velocity button.
    rv_correct_btn = QtGui.QPushButton(rv_tab)
    rv_correct_btn.setObjectName("rv_correct_btn")
    rv_correct_btn.setText("Correct")
    hbox.addWidget(rv_correct_btn)
    rv_settings_vbox.addLayout(hbox)


    # Add a spacer until the big button.
    rv_settings_vbox.addItem(QtGui.QSpacerItem(
        20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding))


    # The cross-correlate and correct button.
    rv_ccc_btn = QtGui.QPushButton(rv_tab)
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
    blank_widget = QtGui.QWidget(rv_tab)
    sp = QtGui.QSizePolicy(
        QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
    sp.setHorizontalStretch(0)
    sp.setVerticalStretch(0)
    sp.setHeightForWidth(blank_widget.sizePolicy().hasHeightForWidth())
    blank_widget.setSizePolicy(sp)
    blank_widget.setObjectName("blank_widget")

    rv_plot = mpl.MPLWidget(blank_widget)
    layout = QtGui.QVBoxLayout(blank_widget)        
    layout.addWidget(rv_plot, 1)

    rv_tab_layout.addWidget(blank_widget)

    rv_plot.axes = rv_plot.figure.add_subplot(311, axisbg="#FFFFFF")
    rv_plot.axes.plot([0, 1], [0, 1],'bo-')
    rv_plot.draw()

    # Add the tab to the application.
    tabs.addTab(rv_tab, "Radial velocity")


    # HACK TODO:


    # Read in the default settings from the SMH session file.


    # Read in these values from the SMH session file.
    rv_wl_region.addItem("Test")
    rv_wl_region.addItem("Test2")    
    rv_wl_region.setItemText(0, "8450-8750 Angstroms")

    rv_norm_function.setCurrentIndex(1) # Indexed from the contents.

    #rv_norm_low_sigma.setText?

    return rv_tab
