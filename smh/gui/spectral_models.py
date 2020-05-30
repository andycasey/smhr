#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A widget for fitting spectral models. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)


import operator
import numpy as np
from PySide2 import (QtCore, QtWidgets as QtGui)
from time import time

import mpl

DOUBLE_CLICK_INTERVAL = 0.1 # MAGIC HACK


class SpectralModelsTableModel(QtCore.QAbstractTableModel):

    header = [" ", "Wavelength", "Element", "EW"]
    attrs = ("is_acceptable", "_repr_wavelength", "_repr_element", 
        "equivalent_width")

    def __init__(self, parent, data, *args):
        super(SpectralModelsTableModel, self).__init__(parent, *args)
        self._data = data

    def rowCount(self, parent):
        return len(self._data)

    def columnCount(self, parent):
        return len(self.header)

    def data(self, index, role):
        if not index.isValid():
            return None

        attr = self.attrs[index.column()]
        try:
            value = getattr(self._data[index.row()], attr)
                
        except AttributeError:
            try:
                result = self._data[index.row()].metadata["fitted_result"] 
                value = result[2][attr]
            except (AttributeError, KeyError):
                value = np.nan

            else:
                if attr == "equivalent_width":
                    value = "{0:.1f} +/- {1:.1f}".format(
                        1000 * value[0],
                        1000 * np.max(np.abs(value[1:])))

        if index.column() == 0:
            if role == QtCore.Qt.CheckStateRole:
                return QtCore.Qt.Checked if value else QtCore.Qt.Unchecked
            else:
                return None

        elif role != QtCore.Qt.DisplayRole:
            return None

        return value

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal \
        and role == QtCore.Qt.DisplayRole:
            return self.header[col]
        return None


    def setData(self, index, value, role=QtCore.Qt.DisplayRole):
        attr = self.attrs[index.column()]
        if attr != "is_acceptable":
            return False

        self._data[index.row()].metadata[attr] = value
        self.dataChanged.emit(index, index)
        return value


    def sort(self, column, order):

        self.emit(QtCore.SIGNAL("layoutAboutToBeChanged()"))
        self._data = sorted(self._data,
            key=lambda sm: getattr(sm, self.attrs[column]))
        
        if order == QtCore.Qt.DescendingOrder:
            self._data.reverse()

        self.dataChanged.emit(self.createIndex(0, 0),
            self.createIndex(self.rowCount(0), self.columnCount(0)))
        self.emit(QtCore.SIGNAL("layoutChanged()"))


    def flags(self, index):
        if not index.isValid(): return
        return QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsUserCheckable




class SpectralModelsWidget(QtGui.QWidget):
    def __init__(self, spectral_models, *args):
        """
        Initialize a widget for inspecting a list of spectral models.

        :param spectral_models:
            A list of spectral models that sub-class from 
            `smh.spectral_models.BaseSpectralModel`.
        """

        super(SpectralModelsWidget, self).__init__(*args)
        self.spectral_models = spectral_models

        self.setGeometry(300, 200, 570, 450)
        self.setWindowTitle("Spectral models")

        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, 
            QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)

        self.table = QtGui.QTableView()
        self.table.setModel(SpectralModelsTableModel(self, spectral_models))
        self.table.resizeColumnsToContents()
        self.table.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.table.setSortingEnabled(True)

        _ = self.table.selectionModel()
        _.selectionChanged.connect(self.row_selected)

        # Create a matplotlib widget.
        blank_widget = QtGui.QWidget(self)
        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sp.setHorizontalStretch(0)
        sp.setVerticalStretch(0)
        sp.setHeightForWidth(blank_widget.sizePolicy().hasHeightForWidth())
        blank_widget.setSizePolicy(sp)
        
        self.mpl_figure = mpl.MPLWidget(blank_widget, 
            tight_layout=True, autofocus=True)

        mpl_layout = QtGui.QVBoxLayout(blank_widget)
        mpl_layout.addWidget(self.mpl_figure, 1)

        layout = QtGui.QHBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.addWidget(self.table)
        layout.addWidget(blank_widget)

        self.setLayout(layout)

        # Initiate the matplotlib figure.
        self.mpl_axis = self.mpl_figure.figure.add_subplot(111)
        
        # Draw the spectrum first. 
        normalized_spectrum = self.spectral_models[0].session.normalized_spectrum
        self.mpl_axis.plot(
            normalized_spectrum.dispersion,
            normalized_spectrum.flux,
            c="k", drawstyle="steps-mid")
        self.mpl_axis.fill_between([], [], [], facecolor="k", alpha=0.5,
            edgecolor=None, zorder=1)

        self.mpl_axis.set_ylim(0, 1.2)
        self.mpl_axis.set_xlabel(r"Wavelength (${\rm \AA}$)")
        self.mpl_axis.set_ylabel(r"Normalized flux")
        self.mpl_figure.draw()

        # Plot the model spectrum.
        self.mpl_axis.plot([], [], c="r", zorder=2)
        self.mpl_axis.fill_between([], [], [], facecolor="r", alpha=0.5,
            edgecolor=None, zorder=1)

        self.mpl_figure.mpl_connect(
            "button_press_event", self.figure_mouse_press)
        self.mpl_figure.mpl_connect(
            "button_release_event", self.figure_mouse_release)

        self._mpl_masked_regions = []
        self._mpl_nearby_lines_masked_regions = []

        return None


    def figure_mouse_press(self, event):
        """
        Mouse button was clicked on the matplotlib figure.

        :param event:
            The matplotlib event.
        """

        if event.dblclick:
            # Double click.

            # If the user has double-clicked in a masked region for the current
            # model, then remove the most recently added match.

            # Otherwise, add a point and wait for a second continuum point?

            show_index = self.table.selectionModel().selectedRows()[0]
            table_model = self.table.model()
            spectral_model = self.table.model()._data[show_index.row()]

            for i, (s, e) in enumerate(spectral_model.metadata["mask"][::-1]):
                if e >= event.xdata >= s:

                    mask = spectral_model.metadata["mask"]
                    index = len(mask) - 1 - i
                    del mask[index]

                    # Re-fit the current spectral_model.
                    spectral_model.fit()

                    # Update the table view for this row.
                    table_model.dataChanged.emit(
                        table_model.createIndex(show_index.row(), 0),
                        table_model.createIndex(
                            show_index.row(), table_model.columnCount(0)))

                    # Update the view of the current model.
                    self.row_selected()

                    break

            else:
                # No match with a masked region. Add a point and wait for a
                # second continuum point?
                None

        else:
            # Single click.

            # Set up/update the excluded region.
            xmin, xmax, ymin, ymax = (event.xdata, np.nan, -1e8, +1e8)
            try:
                self._mask_region
            except AttributeError:
                self._mask_region = self.mpl_axis.axvspan(**{
                    "xmin": xmin,
                    "xmax": xmax,
                    "ymin": ymin,
                    "ymax": ymax,
                    "facecolor": "r",
                    "edgecolor": "none",
                    "alpha": 0.25,
                    "zorder": -1
                })

            else:
                self._mask_region.set_xy([
                    [xmin, ymin],
                    [xmin, ymax],
                    [xmax, ymax],
                    [xmax, ymin],
                    [xmin, ymin]
                ])

            # Set the signal and the time.
            self._mask_region_signal = (
                time(),
                self.mpl_figure.mpl_connect(
                    "motion_notify_event", self.update_mask_region)
                )
        return None


    def figure_mouse_release(self, event):
        """
        Mouse button was released on the matploltib figure.

        :param event:
            The matplotlib event.
        """

        try:
            signal_time, signal_cid = self._mask_region_signal
        except AttributeError:
            return None
        
        xy = self._mask_region.get_xy()
        
        if event.xdata is None:
            # Out of axis; exclude based on the closest axis limit.
            xdata = xy[2, 0]
        else:
            xdata = event.xdata

        # If the two mouse events were within some time interval,
        # then we should not add a mask because those signals were probably
        # part of a double-click event.
        if  time() - signal_time > DOUBLE_CLICK_INTERVAL \
        and np.abs(xy[0,0] - xdata) > 0:
            
            # Get current spectral model.
            show_index = self.table.selectionModel().selectedRows()[0]
            table_model = self.table.model()
            spectral_model = table_model._data[show_index.row()]

            # Add mask metadata.
            spectral_model.metadata["mask"].append([xy[0,0], xy[2, 0]])

            # Re-fit the spectral model.
            spectral_model.fit()

            # Update the table view for this row.
            table_model.dataChanged.emit(
                table_model.createIndex(show_index.row(), 0),
                table_model.createIndex(
                    show_index.row(), table_model.columnCount(0)))

            # Update the view of the spectral model.
            self.row_selected()

        xy[:, 0] = np.nan

        self._mask_region.set_xy(xy)

        self.mpl_figure.mpl_disconnect(signal_cid)
        self.mpl_figure.draw()
        del self._mask_region_signal
        return None


    def update_mask_region(self, event):
        """
        Update the visible selected masked region for the selected spectral
        model. This function is linked to a callback for when the mouse 
        position moves.

        :param event:
            The matplotlib motion event to show the current mouse position.
        """
        if event.xdata is None:
            return
        
        signal_time, signal_cid = self._mask_region_signal
        if time() - signal_time > DOUBLE_CLICK_INTERVAL: 
            
            data = self._mask_region.get_xy()

            # Update xmax.
            data[2:4, 0] = event.xdata
            self._mask_region.set_xy(data)

            self.mpl_figure.draw()

        return None


    def row_selected(self, *args):
        """
        The row which is selected has been changed.
        """

        show_index = self.table.selectionModel().selectedRows()[0]
        model = self.table.model()._data[show_index.row()]

        window = model.metadata["window"]
        wavelengths = model.transitions["wavelength"]
        try:
            lower_wavelength = wavelengths[0]
            upper_wavelength = wavelengths[-1]

        except IndexError:
            # Single row.
            lower_wavelength, upper_wavelength = (wavelengths, wavelengths)

        lower_wavelength = lower_wavelength - window
        upper_wavelength = upper_wavelength + window

        pixel = np.mean(np.diff(model.session.normalized_spectrum.dispersion))
        self.mpl_axis.set_xlim(
            lower_wavelength - pixel,
            upper_wavelength + pixel
        )
        
        # Show masked regions.
        for i, (start, end) in enumerate(model.metadata["mask"]):
            try:
                patch = self._mpl_masked_regions[i]
            except IndexError:
                self._mpl_masked_regions.append(self.mpl_axis.axvspan(
                    np.nan, np.nan, facecolor="r", edgecolor="None",
                    alpha=0.25))
                patch = self._mpl_masked_regions[-1]

            patch.set_xy([
                [start, -1e8],
                [start, +1e8],
                [end,   +1e8],
                [end,   -1e8],
                [start, -1e8]
            ])
            patch.set_visible(True)

        # Hide the masks from other spectral models.
        for each in self._mpl_masked_regions[len(model.metadata["mask"]):]:
            each.set_visible(False)

        # Any result for the selected spectral model?
        try:
            opt, cov, meta = model.metadata["fitted_result"]
        except KeyError:
            # Hide the model data and any masked regions.
            self.mpl_axis.lines[1].set_data([], [])

            for each in self._mpl_nearby_lines_masked_regions:
                each.set_visible(False)
        else:

            # Set the model data.
            self.mpl_axis.lines[1].set_data(meta["model_x"], meta["model_y"])

            # Any regions masked because of enarby lines?
            if "nearby_lines" in meta:
                for i, (_, (start, end)) in enumerate(meta["nearby_lines"]):
                    try:
                        patch = self._mpl_nearby_lines_masked_regions[i]
                    except IndexError:
                        self._mpl_nearby_lines_masked_regions.append(
                            self.mpl_axis.axvspan(np.nan, np.nan,
                                facecolor="b", edgecolor=None, alpha=0.25))
                        patch = self._mpl_nearby_lines_masked_regions[-1]

                    patch.set_xy([
                        [start, -1e8],
                        [start, +1e8],
                        [end,   +1e8],
                        [end,   -1e8],
                        [start, -1e8]
                    ])
                    patch.set_visible(True)

            # Hide the nearby line regions from other spectral models.
            N = len(meta.get("nearby_lines", []))
            for each in self._mpl_nearby_lines_masked_regions[N:]:
                each.set_visible(False)


        self.mpl_figure.draw()

        return True


if __name__ == "__main__":

    import sys

    from smh import linelists, Session, specutils
    transitions = linelists.LineList.read("/Users/arc/research/ges/linelist/vilnius.ew.fe")

    session = Session([
        "/Users/arc/codes/smh/hd44007red_multi.fits",
        "/Users/arc/codes/smh/hd44007blue_multi.fits",
    ])

    session.normalized_spectrum = specutils.Spectrum1D.read(
        "../smh/hd44007-rest.fits")

    session.metadata["line_list"] = transitions

    from smh import spectral_models as sm
    foo = []
    for i in range(len(transitions)):
        if i % 2:
            foo.append(sm.ProfileFittingModel(session, transitions["hash"][[i]]))
        else:
            foo.append(sm.SpectralSynthesisModel(session, transitions["hash"][[i]],
                transitions["elem1"][i]))

    
    app = QtGui.QApplication(sys.argv)
    window = SpectralModelsWidget(foo)
    window.show()
    sys.exit(app.exec_())
    
