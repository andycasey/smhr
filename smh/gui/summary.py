#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The summary tab view for the Spectroscopy Made Hard GUI. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

from PySide import QtCore, QtGui

import mpl


__all__ = ["SummaryTab"]


class SummaryTab(QtGui.QWidget):

    def __init__(self, parent=None):
        """
        Create a summary tab for a SMH analysis session.

        :param parent: [optional]
            The parent widget.
        """

        super(SummaryTab, self).__init__(parent)
        self.parent = parent
        self.setObjectName("summary_tab")

        # Right pane: A MPL figure with two axes, vertically aligned.
        # Left pane:
        # - Name, RA/DEC
        # - Link to search on Vizier/Simbad
        # - Comments/notes (user info?)

        # Set a size policy.
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding,
            QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)

        # Create a top-level horizontal layout to contain a MPL figure and
        # a vertical layout of settings..
        tab_layout = QtGui.QHBoxLayout(self)
        tab_layout.setContentsMargins(10, 10, 10, 10)

        # Create the left hand pane.
        summary_widget = QtGui.QWidget()
        summary_layout = QtGui.QVBoxLayout(summary_widget)

        # Star name label.
        self.star_label = QtGui.QLabel(self)
        self.star_label.setText("(No star)")
        summary_layout.addWidget(self.star_label)

        # Coordinates label.
        self.coordinates_label = QtGui.QLabel(self)
        self.coordinates_label.setText("")
        summary_layout.addWidget(self.coordinates_label)

        # Notes.
        self.summary_notes = QtGui.QPlainTextEdit(self)
        self.summary_notes.setObjectName("summary_notes")
        self.summary_notes.setPlainText("")

        summary_layout.addWidget(self.summary_notes)

        # External sources of information.
        hbox = QtGui.QHBoxLayout()

        # - Simbad
        self.btn_query_simbad = QtGui.QPushButton(self)
        self.btn_query_simbad.setObjectName("btn_query_simbad")
        self.btn_query_simbad.setText("Query Simbad..")
        self.btn_query_simbad.clicked.connect(self.query_simbad)

        hbox.addWidget(self.btn_query_simbad)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        summary_layout.addLayout(hbox)
        tab_layout.addWidget(summary_widget)


        # Create a matplotlib widget in the right hand pane.
        blank_widget = QtGui.QWidget(self)
        blank_widget.setFixedWidth(400)

        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(blank_widget.sizePolicy().hasHeightForWidth())
        blank_widget.setSizePolicy(sp)


        self.figure = mpl.MPLWidget(blank_widget, tight_layout=True)

        layout = QtGui.QVBoxLayout(blank_widget)
        layout.addWidget(self.figure, 1)
        tab_layout.addWidget(blank_widget)

        self.ax_top_comparison = self.figure.figure.add_subplot(211)
        self.ax_bottom_comparison = self.figure.figure.add_subplot(212)

        return None



    def redraw_literature_comparisons(self):
        """
        Update the literature comparison axes.
        """
        raise NotImplementedError


    def redraw_stellar_property(self):
        """
        The stellar parameters have been updated, so update the axes to reflect
        that change.
        """
        raise NotImplementedError



    def query_simbad(self):
        """
        Query Simbad for the star in the current session.
        """

        # Get the positional information.

        # Execute a system call to open a URL with the default browser.

        raise NotImplementedError




