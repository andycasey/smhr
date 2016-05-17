#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The summary tab view for the Spectroscopy Made Hard GUI. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

from PySide import QtCore, QtGui

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

        text_edit = QtGui.QPlainTextEdit(self)
        text_edit.setObjectName("summary_text")
        text_edit.setGeometry(QtCore.QRect(80, 80, 401, 181))
        text_edit.setPlainText("""
- Name, RA/DEC
- Comments/session log
- HR diagram with isochrones, showing a comparison sample read from disk.
- Link to vizier/simbad search for this object""")

        return None

