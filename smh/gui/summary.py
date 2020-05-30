#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The summary tab view for the Spectroscopy Made Hard GUI. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
from os import system
from PySide2 import (QtCore, QtGui as QtGui2, QtWidgets as QtGui)
from urllib.parse import quote

import mpl


__all__ = ["SummaryTab"]

logger = logging.getLogger(__name__)

class SummaryTab(QtGui.QWidget):

    def __init__(self, parent):
        """
        Create a summary tab for a SMH analysis session.

        :param parent:
            The parent widget.
        """

        super(SummaryTab, self).__init__(parent)
        self.parent = parent

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
        tab_layout.setContentsMargins(20, 20, 20, 0)

        # Create the left hand pane.
        summary_widget = QtGui.QWidget()
        summary_layout = QtGui.QVBoxLayout(summary_widget)

        # Star name label.
        self.star_label = QtGui.QLabel(self)
        bold_monospaced = QtGui2.QFont("Courier", 14)
        bold_monospaced.setBold(True)
        self.star_label.setFont(bold_monospaced)
        summary_layout.addWidget(self.star_label)

        # Coordinates label.
        self.coordinates_label = QtGui.QLabel(self)
        self.coordinates_label.setFont(bold_monospaced)
        summary_layout.addWidget(self.coordinates_label)

        # Notes.
        self.summary_notes = QtGui.QPlainTextEdit(self)
        summary_layout.addWidget(self.summary_notes)

        # External sources of information.
        hbox = QtGui.QHBoxLayout()

        # - Simbad
        self.btn_query_simbad = QtGui.QPushButton(self)
        self.btn_query_simbad.setText("Query Simbad..")
        self.btn_query_simbad.clicked.connect(self.query_simbad)

        hbox.addWidget(self.btn_query_simbad)
        hbox.addItem(QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        summary_layout.addLayout(hbox)
        tab_layout.addWidget(summary_widget)


        # Create a matplotlib widget in the right hand pane.


        self.figure = mpl.MPLWidget(None, tight_layout=True)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(self.figure.sizePolicy().hasHeightForWidth())
        self.figure.setSizePolicy(sp)
        self.figure.setFixedWidth(400)
        tab_layout.addWidget(self.figure)

        self.ax_top_comparison = self.figure.figure.add_subplot(211)
        self.ax_bottom_comparison = self.figure.figure.add_subplot(212)

        # Initialize the widgets.
        self._populate_widgets()

        # Connect the widgets.
        self.summary_notes.textChanged.connect(self.update_summary_notes)

        return None


    def _populate_widgets(self):
        """
        Populate the widgets with values, either from the parent session or the
        default session file.
        """

        # Star label.
        try:
            star_label = str(self.parent.session.metadata["OBJECT"]).strip()
        except (AttributeError, KeyError):
            star_label = "(No star)"
        self.star_label.setText(star_label)

        # Star coordinates.
        try:
            ra = str(self.parent.session.metadata["RA"]).strip()
            dec = str(self.parent.session.metadata["DEC"]).strip()
        except (AttributeError, KeyError):
            ra, dec = ("", "")
            # Disable 'Query Simbad..' button since we have no position info.
            self.btn_query_simbad.setEnabled(False)

        else:
            self.btn_query_simbad.setEnabled(True)
            if ":" not in ra and ":" not in dec:
                # Convert to sexagesimal?
                # TODO
                None
                
        self.coordinates_label.setText(" ".join([ra, dec]))

        # Summary notes.
        try:
            summary_notes = self.parent.session.metadata["NOTES"]
        except (AttributeError, KeyError):
            summary_notes = ""    
        self.summary_notes.setPlainText(summary_notes)

        # Blank axes unless there are data.
        for ax in (self.ax_top_comparison, self.ax_bottom_comparison):
            ax.set_xticks([])
            ax.set_yticks([])
        
        return None


    def update_summary_notes(self):
        """
        Update the summary notes in the session metadata to reflect the changes
        made to the summary notes widget.
        """

        if self.parent.session is not None:
            self.parent.session.metadata["NOTES"] \
                = self.summary_notes.toPlainText()
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
        ra = str(self.parent.session.metadata["RA"]).strip()
        dec = str(self.parent.session.metadata["DEC"]).strip()

        # Execute a system call to open a URL with the default browser.
        url = "http://simbad.u-strasbg.fr/simbad/sim-coo?Coord={ra}%20{dec}"\
              "&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames"\
              "=none&Radius=5&Radius.unit=arcsec&submit=submit+query&CoordList="
        url = url.format(ra=quote(ra), dec=quote(dec))

        logger.info("Opening {}".format(url))

        system('open "{}"'.format(url))

        return None


