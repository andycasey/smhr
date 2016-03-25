#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" The summary tab view for the Spectroscopy Made Hard GUI. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

from PySide import QtCore, QtGui

__all__ = ["initialise_tab"]


def initialise_tab(tabs):
    """
    Create a summary tab and add it to the application tabs.

    :param tabs:
        The application tab widget to add a summary tab to.

    :type tabs:
        :class:`QtGui.QTabWidget`
    """

    # Create a tab and add a single text edit object.
    summary_tab = QtGui.QWidget()
    summary_tab.setObjectName("summary_tab")

    text_edit = QtGui.QPlainTextEdit(summary_tab)
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

    tabs.addTab(summary_tab, "Summary")

    return summary_tab