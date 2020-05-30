#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Quality control widget for setting spectral models as unacceptable. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["QualityControlDialog"]

import logging
import numpy as np

from PySide2 import (QtCore, QtGui as QtGui2, QtWidgets as QtGui)

logger = logging.getLogger(__name__)

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

        self.affected_indices = []
        self.session = session
        self.callbacks = callbacks or []
        self.filter_spectral_models = filter_spectral_models

        # Display dialog in center and set size policy.
        self.setGeometry(400, 400, 400, 400)
        desktop = QtGui.QApplication.desktop()
        self.move(desktop.screen().rect().center() \
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

        grid.addWidget(self.make_label("min"), 0, 2, 1, 1)
        grid.addWidget(self.make_label("max"), 0, 3, 1, 1)
        
        index = 1
        grid.addWidget(label_wavelength, index, 0, 1, 1)
        grid.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum), index, 1, 1, 1)
        grid.addWidget(self.edit_wavelength_lower, index, 2, 1, 1)
        grid.addWidget(self.edit_wavelength_upper, index, 3, 1, 1)

        line_list = session.metadata.get("line_list", None)
        if line_list is not None:
            for item in (self.edit_wavelength_lower, self.edit_wavelength_upper):
                item.setValidator(QtGui2.QDoubleValidator(
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

        index += 1
        grid.addWidget(label_expot, index, 0, 1, 1)
        grid.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum), index, 1, 1, 1)
        grid.addWidget(self.edit_expot_lower, index, 2, 1, 1)
        grid.addWidget(self.edit_expot_upper, index, 3, 1, 1)

        try:
            expot_range = (min(line_list["expot"]), max(line_list["expot"]))

        except (TypeError, ValueError):
            expot_range = (-np.inf, +np.inf)

        for item in (self.edit_expot_lower, self.edit_expot_upper):
            item.setValidator(QtGui2.QDoubleValidator(
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

        index += 1
        grid.addWidget(label_ew, index, 0, 1, 1)
        grid.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum), index, 1, 1, 1)
        grid.addWidget(self.edit_ew_lower, index, 2, 1, 1)
        grid.addWidget(self.edit_ew_upper, index, 3, 1, 1)

        for item in (self.edit_ew_lower, self.edit_ew_upper):
            item.setValidator(QtGui2.QDoubleValidator(0, np.inf, 2, item))
            item.textChanged.connect(self.check_lineedit_state)


        # Uncertainty in equivalent width
        label_ew = QtGui.QLabel(self)
        label_ew.setText(u"Equivalent width uncertainty (mÅ)")

        self.edit_ew_u_lower = QtGui.QLineEdit(self)
        self.edit_ew_u_lower.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_ew_u_lower.setMaximumSize(QtCore.QSize(60, 16777215))

        self.edit_ew_u_upper = QtGui.QLineEdit(self)
        self.edit_ew_u_upper.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_ew_u_upper.setMaximumSize(QtCore.QSize(60, 16777215))

        index += 1
        grid.addWidget(label_ew, index, 0, 1, 1)
        grid.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum), index, 1, 1, 1)
        grid.addWidget(self.edit_ew_u_lower, index, 2, 1, 1)
        grid.addWidget(self.edit_ew_u_upper, index, 3, 1, 1)

        for item in (self.edit_ew_u_lower, self.edit_ew_u_upper):
            item.setValidator(QtGui2.QDoubleValidator(0, np.inf, 1, item))
            item.textChanged.connect(self.check_lineedit_state)


        # Uncertainty in equivalent width
        label_ew = QtGui.QLabel(self)
        label_ew.setText(u"Equivalent width uncertainty (%)")

        self.edit_ew_u_percent_lower = QtGui.QLineEdit(self)
        self.edit_ew_u_percent_lower.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_ew_u_percent_lower.setMaximumSize(QtCore.QSize(60, 16777215))

        self.edit_ew_u_percent_upper = QtGui.QLineEdit(self)
        self.edit_ew_u_percent_upper.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_ew_u_percent_upper.setMaximumSize(QtCore.QSize(60, 16777215))

        index += 1
        grid.addWidget(label_ew, index, 0, 1, 1)
        grid.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum), index, 1, 1, 1)
        grid.addWidget(self.edit_ew_u_percent_lower, index, 2, 1, 1)
        grid.addWidget(self.edit_ew_u_percent_upper, index, 3, 1, 1)

        for item in (self.edit_ew_u_percent_lower, self.edit_ew_u_percent_upper):
            item.setValidator(QtGui2.QDoubleValidator(0, 100, 0, item))
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

        index += 1
        grid.addWidget(label_rew, index, 0, 1, 1)
        grid.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum), index, 1, 1, 1)
        grid.addWidget(self.edit_rew_lower, index, 2, 1, 1)
        grid.addWidget(self.edit_rew_upper, index, 3, 1, 1)

        for item in (self.edit_rew_lower, self.edit_rew_upper):
            item.setValidator(QtGui2.QDoubleValidator(-np.inf, np.inf, 2, item))
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

        index += 1
        grid.addWidget(label_abundance, index, 0, 1, 1)
        grid.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum), index, 1, 1, 1)
        grid.addWidget(self.edit_abundance_lower, index, 2, 1, 1)
        grid.addWidget(self.edit_abundance_upper, index, 3, 1, 1)

        for item in (self.edit_abundance_lower, self.edit_abundance_upper):
            item.setValidator(QtGui2.QDoubleValidator(-np.inf, np.inf, 2, item))
            item.textChanged.connect(self.check_lineedit_state)


        # Abundance uncertainty
        label_abundance = QtGui.QLabel(self)
        label_abundance.setText(u"Abundance uncertainty (log ε, dex)")

        self.edit_abundance_u_lower = QtGui.QLineEdit(self)
        self.edit_abundance_u_lower.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_abundance_u_lower.setMaximumSize(QtCore.QSize(60, 16777215))

        self.edit_abundance_u_upper = QtGui.QLineEdit(self)
        self.edit_abundance_u_upper.setMinimumSize(QtCore.QSize(60, 0))
        self.edit_abundance_u_upper.setMaximumSize(QtCore.QSize(60, 16777215))

        index += 1
        grid.addWidget(label_abundance, index, 0, 1, 1)
        grid.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum), index, 1, 1, 1)
        grid.addWidget(self.edit_abundance_u_lower, index, 2, 1, 1)
        grid.addWidget(self.edit_abundance_u_upper, index, 3, 1, 1)

        for item in (self.edit_abundance_u_lower, self.edit_abundance_u_upper):
            item.setValidator(QtGui2.QDoubleValidator(-np.inf, np.inf, 2, item))
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
            QtGui2.QValidator.Acceptable: 'none',        # Normal background
            QtGui2.QValidator.Intermediate: "#FFF79A",   # Yellow
        }.get(state, "#F6989D")                         # Red

        sender.setStyleSheet("QLineEdit {{ background-color: {} }}".format(color))
    
        return None


    def apply(self):
        """
        Apply the specified quality constraints to the parent session.
        """

        def safe_float(lineedit_widget):
            try:
                return float(lineedit_widget.text())
            except ValueError:
                return None

        constraints = {
            "wavelength": [
                safe_float(self.edit_wavelength_lower),
                safe_float(self.edit_wavelength_upper)
            ],
            "abundance": [
                safe_float(self.edit_abundance_lower),
                safe_float(self.edit_abundance_upper)
            ],
            "abundance_uncertainty": [
                safe_float(self.edit_abundance_u_lower),
                safe_float(self.edit_abundance_u_upper)
            ],
            "equivalent_width": [
                safe_float(self.edit_ew_lower),
                safe_float(self.edit_ew_upper)
            ],
            "equivalent_width_uncertainty": [
                safe_float(self.edit_ew_u_lower),
                safe_float(self.edit_ew_u_upper)
            ],
            "equivalent_width_percentage_uncertainty": [
                safe_float(self.edit_ew_u_percent_lower),
                safe_float(self.edit_ew_u_percent_upper)
            ],
            "reduced_equivalent_width": [
                safe_float(self.edit_rew_lower),
                safe_float(self.edit_rew_upper)
            ],
            "excitation_potential": [
                safe_float(self.edit_expot_lower),
                safe_float(self.edit_expot_upper)
            ]
        }

        logger.info("Supplying constraints: {}".format(constraints))
        affected, self.affected_indices \
            = self.session.apply_spectral_model_quality_constraints(
                constraints, only=self.filter_spectral_models, full_output=True)

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

    
    def make_label(self, text):
        label = QtGui.QLabel(self)
        label.setText(text)
        return label


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

    import smh
    a = smh.Session.load("hd122563.smh")
    window = QualityControlDialog(a)
    window.exec_()

    

