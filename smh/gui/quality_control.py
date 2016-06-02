#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Quality control widget for setting spectral models as unacceptable. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["QualityControlDialog"]

import numpy as np
from PySide import QtCore, QtGui


class QualityControlDialog(QtGui.QDialog):

    _default_explanation = \
        "Set any lower and upper constraints on line properties. Any lines "\
        "that do not meet the specified constraints will be marked as "\
        "unacceptable."

    def __init__(self, session, filter_spectral_models=None, callbacks=None,
        explanation=None, **kwargs):
        """
        A widget to identify spectral models that do not meet some quality
        constraints, and to mark them as unacceptable.

        :param session:
            A session with spectral models to filter.

        :param filter_spectral_models: [optional]
            A function that returns True if a given spectral model should be
            checked to see whether it meets the specified quality criteria.
            If None is supplied, then all spectral models will be checked.

        :param callbacks: [optional]
            Callback functions to execute when this dialog is closed.

        :param explanation: [optional]
            Set the text in the explanatory header of the dialog.
        """

        super(QualityControlDialog, self).__init__(**kwargs)

        self.session = session
        self.callbacks = callbacks or []
        self.filter_spectral_models = filter_spectral_models

        # Display dialog in center and set size policy.
        self.setGeometry(400, 300, 400, 300)
        self.move(QtGui.QApplication.desktop().screen().rect().center() \
            - self.rect().center())
        self.setWindowTitle("Quality criteria for spectral models")

        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, 
            QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)

        vbox = QtGui.QVBoxLayout(self)

        # Explanation header.
        explanatory_label = QtGui.QLabel(self)
        explanatory_label.setWordWrap(True)
        explanatory_label.setText(explanation or self._default_explanation)
        vbox.addWidget(explanatory_label)

        # Grid layout.
        grid = QtGui.QGridLayout()

        # Wavelength.
        label_wavelength = QtGui.QLabel(self)
        label_wavelength.setText(u"Wavelength (Å)")

        self.edit_wavelength_lower = QtGui.QLineEdit(self)
        self.edit_wavelength_lower.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_wavelength_lower.setMaximumSize(QtCore.QSize(60, 16777215))

        self.edit_wavelength_upper = QtGui.QLineEdit(self)
        self.edit_wavelength_upper.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_wavelength_upper.setMaximumSize(QtCore.QSize(60, 16777215))

        grid.addWidget(label_wavelength, 0, 0, 1, 1)
        grid.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum), 0, 1, 1, 1)
        grid.addWidget(self.edit_wavelength_lower, 0, 2, 1, 1)
        grid.addWidget(self.edit_wavelength_upper, 0, 3, 1, 1)

        line_list = session.metadata.get("line_list", None)
        if line_list is not None:
            for item in (self.edit_wavelength_lower, self.edit_wavelength_upper):
                item.setValidator(QtGui.QDoubleValidator(
                    min(line_list["wavelength"]),
                    max(line_list["wavelength"]),
                    3, item))
                item.textChanged.connect(self.check_lineedit_state)


        # Excitation potential.
        label_expot = QtGui.QLabel(self)
        label_expot.setText("Excitation potential (eV)")

        self.edit_expot_lower = QtGui.QLineEdit(self)
        self.edit_expot_lower.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_expot_lower.setMaximumSize(QtCore.QSize(60, 16777215))

        self.edit_expot_upper = QtGui.QLineEdit(self)
        self.edit_expot_upper.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_expot_upper.setMaximumSize(QtCore.QSize(60, 16777215))

        grid.addWidget(label_expot, 1, 0, 1, 1)
        grid.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum), 1, 1, 1, 1)
        grid.addWidget(self.edit_expot_lower, 1, 2, 1, 1)
        grid.addWidget(self.edit_expot_upper, 1, 3, 1, 1)

        try:
            expot_range = (min(line_list["expot"]), max(line_list["expot"]))

        except (TypeError, ValueError):
            expot_range = (-np.inf, +np.inf)

        for item in (self.edit_expot_lower, self.edit_expot_upper):
            item.setValidator(QtGui.QDoubleValidator(
                expot_range[0], expot_range[1], 3, item))
            item.textChanged.connect(self.check_lineedit_state)


        # Equivalent width
        label_ew = QtGui.QLabel(self)
        label_ew.setText(u"Equivalent width (mÅ)")

        self.edit_ew_lower = QtGui.QLineEdit(self)
        self.edit_ew_lower.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_ew_lower.setMaximumSize(QtCore.QSize(60, 16777215))

        self.edit_ew_upper = QtGui.QLineEdit(self)
        self.edit_ew_upper.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_ew_upper.setMaximumSize(QtCore.QSize(60, 16777215))

        grid.addWidget(label_ew, 2, 0, 1, 1)
        grid.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum), 2, 1, 1, 1)
        grid.addWidget(self.edit_ew_lower, 2, 2, 1, 1)
        grid.addWidget(self.edit_ew_upper, 2, 3, 1, 1)

        for item in (self.edit_ew_lower, self.edit_ew_upper):
            item.setValidator(QtGui.QDoubleValidator(0, np.inf, 2, item))
            item.textChanged.connect(self.check_lineedit_state)


        # Reduced equivalent width
        label_rew = QtGui.QLabel(self)
        label_rew.setText(u"Reduced equivalent width")

        self.edit_rew_lower = QtGui.QLineEdit(self)
        self.edit_rew_lower.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_rew_lower.setMaximumSize(QtCore.QSize(60, 16777215))

        self.edit_rew_upper = QtGui.QLineEdit(self)
        self.edit_rew_upper.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_rew_upper.setMaximumSize(QtCore.QSize(60, 16777215))

        grid.addWidget(label_rew, 3, 0, 1, 1)
        grid.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum), 3, 1, 1, 1)
        grid.addWidget(self.edit_rew_lower, 3, 2, 1, 1)
        grid.addWidget(self.edit_rew_upper, 3, 3, 1, 1)

        for item in (self.edit_rew_lower, self.edit_rew_upper):
            item.setValidator(QtGui.QDoubleValidator(-np.inf, np.inf, 2, item))
            item.textChanged.connect(self.check_lineedit_state)


        # Abundance
        label_abundance = QtGui.QLabel(self)
        label_abundance.setText(u"Abundance (log ε, dex)")

        self.edit_abundance_lower = QtGui.QLineEdit(self)
        self.edit_abundance_lower.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_abundance_lower.setMaximumSize(QtCore.QSize(60, 16777215))

        self.edit_abundance_upper = QtGui.QLineEdit(self)
        self.edit_abundance_upper.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_abundance_upper.setMaximumSize(QtCore.QSize(60, 16777215))

        grid.addWidget(label_abundance, 4, 0, 1, 1)
        grid.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum), 4, 1, 1, 1)
        grid.addWidget(self.edit_abundance_lower, 4, 2, 1, 1)
        grid.addWidget(self.edit_abundance_upper, 4, 3, 1, 1)

        for item in (self.edit_abundance_lower, self.edit_abundance_upper):
            item.setValidator(QtGui.QDoubleValidator(-np.inf, np.inf, 2, item))
            item.textChanged.connect(self.check_lineedit_state)

        vbox.addLayout(grid)

        hbox = QtGui.QHBoxLayout()
        hbox.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, 
            QtGui.QSizePolicy.Minimum))

        self.btn_apply = QtGui.QPushButton(self)
        self.btn_apply.setText("Apply")
        self.btn_apply.setDefault(True)
        self.btn_apply.clicked.connect(self.apply)
        hbox.addWidget(self.btn_apply)

        self.btn_cancel = QtGui.QPushButton(self)
        self.btn_cancel.setText("Cancel")
        self.btn_cancel.clicked.connect(self.close)
        hbox.addWidget(self.btn_cancel)

        vbox.addLayout(hbox)

        return None


    def check_lineedit_state(self, *args, **kwargs):
        """
        Update the background color of a QLineEdit object based on whether the
        input is valid.
        """

        sender = self.sender()
        state = sender.validator().validate(sender.text(), 0)[0]

        color = {
            QtGui.QValidator.Acceptable: 'none',        # Normal background
            QtGui.QValidator.Intermediate: "#FFF79A",   # Yellow
        }.get(state, "#F6989D")                         # Red

        sender.setStyleSheet("QLineEdit {{ background-color: {} }}".format(color))
    
        return None


    def apply(self):
        """
        Apply the specified quality constraints to the parent session.
        """

        def safe_float(lineedit_widget):
            try:
                return float(lineedit_widget.getText())
            except:
                return None

        constraints = {
            "wavelength": [
                safe_float(self.edit_wavelength_lower),
                safe_float(self.edit_wavelength_upper)
            ],
        }

        raise NotImplementedError("requires a thinko..")
        
        affected = self.session.apply_quality_constraints(
            constraints, only=self.filter_spectral_models)

        # Show how many were affected.
        return self.show_affected(affected)


    def show_affected(self, N):
        """
        Show a dialog window outlining how many spectral models were affected
        by the specified quality constraints, then close this dialog.

        :param N:
            The number of spectral models affected.
        """

        QtGui.QMessageBox.information(self, "Information",
            "There {was_were} {N} spectral model{s} affected by these quality "
            "constraints.".format(N=N, was_were="was" if N == 1 else "were",
                s="" if N == 1 else "s"))

        return self.close()


    def closeEvent(self, event):
        """
        Perform any requested callbacks before letting the widget close.

        :param event:
            The close event.
        """

        for callback in self.callbacks:
            callback()
        event.accept()
        return None



if __name__ == "__main__":

    # This is just for development testing.
    import sys
    try:
        app = QtGui.QApplication(sys.argv)
    except RuntimeError:
        None
    window = QualityControlDialog(None)
    window.exec_()

    

