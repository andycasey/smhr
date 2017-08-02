#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
from PySide import QtCore, QtGui
import sys, os

import matplotlib

from smh.gui import mpl
from smh.gui.base import SMHSpecDisplay, MeasurementTableModelBase
from smh import Session

datadir = os.path.dirname(os.path.abspath(__file__))+'/test_data'

app = QtGui.QApplication(sys.argv)

def create_blank_window(width=800, height=400):
    window = QtGui.QMainWindow()
    window.resize(width, height)
    cw = QtGui.QWidget(window)
    window.setCentralWidget(cw)
    return window, cw

def test_MPLWidget():
    window, cw = create_blank_window()
    fig = mpl.MPLWidget(None, toolbar=False, autofocus=True)
    ax = fig.figure.add_subplot(1,1,1)
    ax.plot([1,2,3],[1,0,1],'k')
    label = QtGui.QLabel("Test")

    vbox = QtGui.QVBoxLayout(cw)
    vbox.addWidget(fig)
    vbox.addWidget(label)
    return window

def test_SMHSpecDisplay():
    session = Session.load(datadir+"/test_G64-12.smh")
    window, cw = create_blank_window()
    
    vbox = QtGui.QVBoxLayout(cw)
    
    fig = SMHSpecDisplay(None, session,
                         enable_masks=True)
    
    vbox.addWidget(fig)
    vbox.setSpacing(0)
    #vbox.setContentsMargins(0,0,0,0)

    models = session.metadata["spectral_models"]
    model = models[0]
    model.fit()
    fig.set_selected_model(model)
    fig.update_spectrum_figure()
    return window

if __name__=="__main__":
    w1 = test_MPLWidget()
    w2 = test_SMHSpecDisplay()

    ws = [w1,w2]
    ws = [w2]

    for w in ws:
        w.show()
        w.raise_()
    sys.exit(app.exec_())
