#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
from PySide import QtCore, QtGui
import sys, os

import matplotlib

from smh.gui import mpl
from smh.gui.base import SMHSpecDisplay
from smh.gui.base import MeasurementTableView, MeasurementTableModelBase, MeasurementTableModelProxy
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
    session = Session.load(datadir+"/test_G64-12_v02.smh")
    window, cw = create_blank_window()
    
    vbox = QtGui.QVBoxLayout(cw)
    
    def get_selected_model():
        return session.spectral_models[0]
    fig = SMHSpecDisplay(None, session,
                         enable_masks=True,
                         get_selected_model=get_selected_model)
    
    
    vbox.addWidget(fig)
    vbox.setSpacing(0)
    #vbox.setContentsMargins(0,0,0,0)

    models = session.metadata["spectral_models"]
    model = models[0]
    model.fit()
    fig.set_selected_model(model)
    fig.update_spectrum_figure()
    return window

def test_SMHSpecDisplay_new_session():
    window, cw = create_blank_window()
    vbox = QtGui.QVBoxLayout(cw)
    fig = SMHSpecDisplay(None, None,
                         enable_masks=True)
    vbox.addWidget(fig)
    vbox.setSpacing(0)
    #vbox.setContentsMargins(0,0,0,0)

    session = Session.load(datadir+"/test_G64-12_v02.smh")
    m = session.metadata["spectral_models"].pop(0)
    fig.new_session(session)
    def get_selected_model():
        return session.spectral_models[0]
    # A bit of a hack
    fig.get_selected_model = get_selected_model

    models = session.metadata["spectral_models"]
    model = models[0]
    model.fit()
    fig.set_selected_model(model)
    fig.update_spectrum_figure()
    return window

def test_MeasurementTableView():
    session = Session.load(datadir+"/test_G64-12_v02.smh")
    session.measure_abundances()
    window, cw = create_blank_window()
    vbox = QtGui.QVBoxLayout(cw)

    columns = ["is_acceptable","wavelength", "expot", "species", "loggf", "elements",
               "equivalent_width","equivalent_width_uncertainty",
               "reduced_equivalent_width",
               "abundances","abundances_to_solar","abundance_uncertainties",
               "is_upper_limit","user_flag"]
    measurement_table = MeasurementTableModelBase(None, session, columns)
    view = MeasurementTableView(None, session)
    view.setModel(measurement_table)

    vbox.addWidget(view)
    return window

#### The proxy model is super slow if you emit dataChanged.
#### You can avoid this by explicitly connecting the view to update it.
#### It also causes segfaults for no apparent reason...so I'm commenting it out
#def test_MeasurementTableView_Sortable():
#    session = Session.load(datadir+"/test_G64-12_v02.smh")
#    session.measure_abundances()
#    window, cw = create_blank_window()
#    vbox = QtGui.QVBoxLayout(cw)
#
#    columns = ["is_acceptable","wavelength", "expot", "species", "loggf", "elements",
#               "equivalent_width","equivalent_width_uncertainty",
#               "reduced_equivalent_width",
#               "abundances","abundances_to_solar","abundance_uncertainties",
#               "is_upper_limit","user_flag"]
#    measurement_table = MeasurementTableModelBase(None, session, columns)
#    sort_and_filter = MeasurementTableModelProxy(None)
#    sort_and_filter.setSourceModel(measurement_table)
#    sort_and_filter.add_filter_function("test", lambda x: x.wavelength > 5000)
#    
#    view = MeasurementTableView(None)
#    view.setModel(sort_and_filter)
#    sort_and_filter.add_view_to_update(view)
#
#    vbox.addWidget(view)
#    return window

if __name__=="__main__":
    #w1 = test_MPLWidget()
    #w2 = test_SMHSpecDisplay()
    #w2b = test_SMHSpecDisplay_new_session()
    w3 = test_MeasurementTableView()
    w4 = test_MeasurementTableView_Sortable()

    #ws = [w1,w2,w2b,w3]
    ws = [w3, w4]

    for w in ws:
        w.show()
        w.raise_()
    sys.exit(app.exec_())
