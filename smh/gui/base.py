from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import sys
import logging
import os
from PySide import QtCore, QtGui

class SMHWidgetBase(QtGui.QWidget):
    """
    Base class for interactive widgets
    """
    def __init__(self, parent, session=None, widgets_to_update = []):
        super(SMHWidgetBase, self).__init__(parent)
        self.parent = parent
        self.widgets_to_update = widgets_to_update
        self.session = session
        ## TODO set style options here
    def reload_from_session(self):
        """
        Rewrite any internal information with session data.
        """
        raise NotImplementedError
    def send_to_session(self):
        """
        Save any internal cached information to session (if any).
        """
        raise NotImplementedError
    def update_widgets_after_selection(self, selected_models):
        """
        Call update_after_selection for all widgets_to_update
        """
        for widget in self.widgets_to_update:
            widget.update_after_selection(selected_models)
    def update_widgets_after_measurement_change(self, changed_model):
        """
        Call update_after_measurement_change for all widgets_to_update
        """
        for widget in self.widgets_to_update:
            widget.update_after_measurement_change(changed_model)

class SMHSpecDisplay(SMHWidgetBase):
    def __init__(self, parent, session=None, widgets_to_update = []):
        super(SMHSpecDisplay, self).__init__(parent, session, widgets_to_update)
    def update_after_selection(self, selected_models):
        raise NotImplementedError
    def update_after_measurement_change(self, changed_model):
        raise NotImplementedError


class SMHFittingOptions(SMHWidgetBase):
    def __init__(self, parent, session=None, widgets_to_update = []):
        super(SMHFittingOptions, self).__init__(parent, session, widgets_to_update)
    def update_after_selection(self, selected_models):
        raise NotImplementedError
    def update_after_measurement_change(self, changed_model):
        raise NotImplementedError
    

class SMHMeasurementList(SMHWidgetBase):
    def __init__(self, parent, session=None, widgets_to_update = []):
        super(SMHMeasurementList, self).__init__(parent, session, widgets_to_update)
    def update_after_selection(self, selected_models):
        raise NotImplementedError
    def update_after_measurement_change(self, changed_model):
        raise NotImplementedError
    

class SMHScatterplot(SMHWidgetBase):
    def __init__(self, parent, session=None, widgets_to_update = []):
        super(SMHScatterplot, self).__init__(parent, session, widgets_to_update)
    def update_after_selection(self, selected_models):
        raise NotImplementedError
    def update_after_measurement_change(self, changed_model):
        raise NotImplementedError
