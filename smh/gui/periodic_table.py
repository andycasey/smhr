#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A widget to select element(s). """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["PeriodicTableDialog"]

from numpy import ceil
from PySide2 import (QtCore, QtWidgets as QtGui)
from textwrap import dedent

class PeriodicTableDialog(QtGui.QDialog):

    elements = """
        H                                                  He
        Li Be                               B  C  N  O  F  Ne
        Na Mg                               Al Si P  S  Cl Ar
        K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
        Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
        Cs Ba Lu Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
        Fr Ra Lr Rf

              La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb
              Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No"""

    element_length = 3 # Includes spaces between elements.


    def __init__(self, selectable_elements=None, explanation=None, 
        multiple_select=False, close_on_first_select=True, **kwargs):
        """
        Show a widget that will let the user select element(s) from the
        specified list.

        :param elements: [None] 
            A list of elements that the user can select from. If None is given
            then all elements can be selected.

        :param explanation: [optional]
            Provide explanatory text at the top of the widget.

        :param multiple_select: [optional]
            Allow multiple elements to be selected.

        :param close_on_first_select: [optional]
            If `multiple_select` is set to False and the this option is set to
            True, then the dialog will be closed the first time an element is
            selected.
        """

        super(PeriodicTableDialog, self).__init__(**kwargs)

        self._offset = 1
        self._selected_elements = []
        self.multiple_select = multiple_select
        self.close_on_first_select = close_on_first_select

        if selectable_elements is None:
            selectable_elements = self.elements.split()

        else:
            selectable_elements \
                = [se.title().strip() for se in selectable_elements]

        # Display dialog in center and set size policy.
        self.setGeometry(800, 500, 800, 500)
        desktop = QtGui.QApplication.desktop()
        self.move(desktop.screen().rect().center() \
            - self.rect().center())
        self.setWindowTitle(
            "Select element{s}".format(s="s" if multiple_select else ""))

        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, 
            QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)

        vbox = QtGui.QVBoxLayout(self)

        if explanation is not None:
            self._offset += 1
            explanatory_label = QtGui.QLabel(self)
            explanatory_label.setWordWrap(True)
            explanatory_label.setText(explanation)
            vbox.addWidget(explanatory_label)

        # Add grid layout.
        grid = QtGui.QGridLayout()

        # How many columns and rows of axes?
        N_rows = len(self.elements.split("\n"))
        N_cols = max([int(ceil(float(len(row))/self.element_length)) \
            for row in dedent(self.elements).split("\n")])

        size = 30

        # We go from 1: so that the periodic_table variable can be formatted in
        # a human-readable way.
        for i, row in enumerate(dedent(self.elements).split("\n")[1:]):

            print(i, row)
            for j in range(N_cols):
                element = row[j*self.element_length:(j + 1)*self.element_length]
                element = element.strip()

                if element:    
                    label = QtGui.QPushButton(self)
                    label.setMinimumSize(QtCore.QSize(size, size))
                    label.setMaximumSize(QtCore.QSize(size, size))
                    label.setObjectName("element_{}".format(element))
                    label.setFlat(True)
                    label.setText(element)

                    if element in selectable_elements:
                        label.clicked.connect(self.toggle_element)
                    else:
                        label.setEnabled(False)

                    grid.addWidget(label, i, j, 1, 1)

            if not row.strip(): # Blank row.
                grid.addItem(
                    QtGui.QSpacerItem(size, size, 
                        QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum),
                    i, 0, 1, 1)

        vbox.addLayout(grid)

        # Add an 'OK' button if this is a multiple selection dialog.
        if multiple_select:
            hbox = QtGui.QHBoxLayout()
            hbox.addItem(QtGui.QSpacerItem(10, 10, 
                QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Minimum))

            self.btn_ok = QtGui.QPushButton(self)
            self.btn_ok.setText("OK")
            self.btn_ok.setEnabled(False)
            self.btn_ok.clicked.connect(self.close)
            hbox.addWidget(self.btn_ok)

            vbox.addLayout(hbox)

        return None


    @property
    def selected_elements(self):
        
        # This could/should be a property, but I am following Qt coding
        # practices and making it a function.
        return self._selected_elements


    def toggle_element(self):
        """ An element was selected in the GUI. """

        element = self.sender().text()

        if element not in self.selected_elements:
            self.select_element(element)
        else:
            self.deselect_element(element)

        return None


    def _get_element_widget(self, element):
        index = self.elements.split().index(element)
        widget = self.children()[self._offset + index]

        assert widget.text() == element
        return widget


    def select_element(self, element):
        """
        Show an element as selected.
        """

        widget = self._get_element_widget(element)
        widget.setStyleSheet("QPushButton { font-weight: bold; }")

        if self.multiple_select:
            self._selected_elements.append(element)
            self.btn_ok.setEnabled(True)

        else:
            self._selected_elements = [element]
            if self.close_on_first_select:
                self.close()

        return None


    def deselect_element(self, element):
        """
        Show an element as being de-selected.
        """

        widget = self._get_element_widget(element)
        widget.setStyleSheet("QPushButton { font-weight: normal; }")

        try:
            self._selected_elements.remove(element)

        except ValueError:
            None

        if len(self._selected_elements) == 0 and self.multiple_select:
            self.btn_ok.setEnabled(False)

        return None


    def closeEvent(self, event):
        """
        The widget cannot be closed until at least one element is selected.

        :param event:
            The close event.
        """

        if not self.selected_elements:
            event.ignore()
            return False

        event.accept()
        return None


    


if __name__ == "__main__":

    # This is just for development testing.
    import sys
    try:
        app = QtGui.QApplication(sys.argv)
    except RuntimeError:
        None
    window = PeriodicTableDialog(["Fe", "Ti"], multiple_select=True,
        explanation="Do something please")
    window.exec_()

