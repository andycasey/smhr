#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The summary tab view for the Spectroscopy Made Hard GUI. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

from PySide import QtCore, QtGui

__all__ = ["SummaryTab"]


class SummaryTab(QtGui.QWidget):

    def __init__(self, parent=None):
        super(SummaryTab, self).__init__(parent)

        self.setObjectName("summary_tab")

        text_edit = QtGui.QPlainTextEdit(self)
        text_edit.setObjectName("summary_text")
        text_edit.setGeometry(QtCore.QRect(80, 80, 401, 181))
        text_edit.setPlainText(
"""
- position
- UTdate
- stellar params compared to some default comparison sample?
- stellar params on an isochrone?
- simbad/DSS image cutout?
- link to vizier for this object? (or retrieve info from vizier?)""")

        return None

