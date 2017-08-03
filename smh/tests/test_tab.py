#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
from PySide import QtCore, QtGui
import sys, os
from smh import Session
from smh.gui.review import ReviewTab

datadir = os.path.dirname(os.path.abspath(__file__))+'/test_data'
try:
    app = QtGui.QApplication(sys.argv)
except RuntimeError:
    pass
def create_blank_window(width=800, height=400):
    window = QtGui.QMainWindow()
    window.resize(width, height)
    cw = QtGui.QWidget(window)
    window.setCentralWidget(cw)
    return window, cw

def test_ReviewTab(session=None):
    if session is None: session = Session.load(datadir+"/test_G64-12_v02.smh")
    session.measure_abundances()
    window, cw = create_blank_window(width=1000,height=600)
    window.session = session
    tab = ReviewTab(window)

    vbox = QtGui.QVBoxLayout(cw)
    vbox.addWidget(tab)
    return window

if __name__=="__main__":
    session = Session.load(datadir+"/test_G64-12_v02.smh")
    w1 = test_ReviewTab(session)
    ws = [w1]

    for w in ws:
        w.show()
        w.raise_()
    sys.exit(app.exec_())

