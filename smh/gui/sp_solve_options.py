
#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Dialog to show options for solving stellar parameters. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["SolveOptionsDialog"]

import logging

from PySide2 import (QtCore, QtWidgets as QtGui)

logger = logging.getLogger(__name__)


class SolveOptionsDialog(QtGui.QDialog):

    def __init__(self, session, **kwargs):
        """
        A widget to show solve options when determining the stellar parameters.

        :param session:
            A session.
        """

        super(SolveOptionsDialog, self).__init__(**kwargs)

        self.session = session

        # Display dialog in center and set size policy.
        self.setGeometry(640, 480, 640, 480)
        desktop = QtGui.QApplication.desktop()
        self.move(desktop.screen().rect().center() \
            - self.rect().center())
        self.setWindowTitle("Solve options for stellar parameter determination")

        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, 
            QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)

        vbox = QtGui.QVBoxLayout(self)

        self.check_use_uncertainties_in_line_fits = QtGui.QCheckBox(self)
        self.check_use_uncertainties_in_line_fits.setText(
            "Use abundance uncertainties in line fits")
        vbox.addWidget(self.check_use_uncertainties_in_line_fits)

        vbox.addItem(QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, 
            QtGui.QSizePolicy.Expanding))

        hbox = QtGui.QHBoxLayout()
        self.btn_save_as_default = QtGui.QPushButton(self)
        self.btn_save_as_default.setText("Save settings as default")
        self.btn_save_as_default.clicked.connect(self.save_as_default)
        hbox.addWidget(self.btn_save_as_default)


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

        self.populate_widgets()

        return None


    def populate_widgets(self):
        """ Populate the dialog widgets with information from the session. """

        # Use abundance uncertainties in line fit.
        self.check_use_uncertainties_in_line_fits.setCheckState(
            QtCore.Qt.Checked   if self.session.setting((
                                    "stellar_parameter_inference",
                                    "use_abundance_uncertainties_in_line_fits"),
                                    True) \
                                else QtCore.Qt.Unchecked)

        return None


    def apply(self):
        """ Apply these settings to the parent session. """

        meta = self.session.metadata
        meta.setdefault("stellar_parameter_inference", {})

        meta["stellar_parameter_inference"]\
            ["use_abundance_uncertainties_in_line_fits"] \
                = self.check_use_uncertainties_in_line_fits.isChecked()

        return self.close()


    def save_as_default(self):
        """ Save the dialog settings to the local default session file. """ 

        self.session.update_default_setting((
            "stellar_parameter_inference",
            "use_abundance_uncertainties_in_line_fits"),
        self.check_use_uncertainties_in_line_fits.isChecked())


        QtGui.QMessageBox.information(self, "Saved", 
            "Settings saved to local session default file.")

        return True
