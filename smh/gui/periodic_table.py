#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A widget to select element(s). """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["PeriodicTableDialog"]

import numpy as np
from PySide import QtCore, QtGui
from textwrap import dedent

import mpl

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
        multiple_select=False, **kwargs):
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
        """

        super(PeriodicTableDialog, self).__init__(**kwargs)

        if selectable_elements is None:
            selectable_elements = self.elements.split()

        else:
            selectable_elements \
                = [se.title().strip() for se in selectable_elements]

        # Display dialog in center and set size policy.
        self.setGeometry(800, 500, 800, 500)
        self.move(QtGui.QApplication.desktop().screen().rect().center() \
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
            explanatory_label = QtGui.QLabel(self)
            explanatory_label.setWordWrap(True)
            explanatory_label.setText(explanation)
            vbox.addWidget(explanatory_label)

        # Add matplotlib widget.
        self.figure = mpl.MPLWidget(None)
        vbox.addWidget(self.figure)

        # How many columns and rows of axes?
        N_rows = len(self.elements.split("\n")) - 1
        N_cols = max([int(np.ceil(float(len(row))/self.element_length)) \
            for row in dedent(self.elements).split("\n")])

        axes = []
        for i, row in enumerate(dedent(self.elements).split("\n")[1:]):

            for j in range(N_cols):
                element = row[j*self.element_length:(j + 1)*self.element_length]
                if element.strip():
                    k = (j + 1) + (i * N_cols)
                    axes.append([
                        element.strip(),
                        self.figure.figure.add_subplot(N_rows, N_cols, k)
                    ])


        for element, ax in axes:
            ax.set_xticks([])
            ax.set_yticks([])

            ax.text(0.5, 0.5, element, verticalalignment="center",
                horizontalalignment="center", transform=ax.transAxes,
                color="#000000" if element in selectable_elements else "#666666")

        self.figure.figure.subplots_adjust(hspace=0, wspace=0)

        self.figure.draw()
        self.figure.mpl_connect("button_press_event", self.figure_mouse_press)

        return None


    def figure_mouse_press(self, event):
        """
        The mouse button has been pressed in the figure.

        :param event:
            The matplotlib event.
        """

        print("event", event.__dict__)



if __name__ == "__main__":

    # This is just for development testing.
    import sys
    try:
        app = QtGui.QApplication(sys.argv)
    except RuntimeError:
        None
    window = PeriodicTableDialog(None)
    window.exec_()

