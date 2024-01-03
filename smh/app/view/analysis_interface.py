# coding:utf-8
from typing import List

# coding:utf-8
# coding:utf-8
from PyQt5.QtCore import QEvent, Qt, QUrl
from PyQt5.QtGui import QColor, QDesktopServices, QIcon, QPainter, QPen
from PyQt5.QtWidgets import (QApplication, QCompleter, QFrame, QHBoxLayout, QFileDialog,
                             QLabel, QSizePolicy, QStackedWidget, QVBoxLayout,
                             QWidget)
from qfluentwidgets import (BodyLabel, BreadcrumbBar, CaptionLabel, CheckBox,
                            ComboBox, EditableComboBox, FlowLayout, FluentIcon,
                            IconWidget, Pivot, PushButton, ScrollArea,
                            StrongBodyLabel, TabBar, TabCloseButtonDisplayMode,
                            LineEdit, PrimaryPushButton,
                            Theme, TitleLabel, ToolButton, ToolTipFilter,
                            applyThemeColor, isDarkTheme, qrouter, toggleTheme)

from ..common.config import EXAMPLE_URL, FEEDBACK_URL, HELP_URL, cfg
from ..common.icon import Icon
from ..common.signal_bus import signalBus
from ..common.style_sheet import StyleSheet
from ..common.translator import Translator
from ..common.trie import Trie
from .gallery_interface import SeparatorWidget

import time

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import \
    NavigationToolbar2QT as BaseNavigationToolbar
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.figure import Figure
from qfluentwidgets import CardWidget, FlowLayout, IconWidget, TextWrap

from .dialog_interface import DialogInterface
import matplotlib.pyplot as plt
from matplotlib.backend_bases import NavigationToolbar2

from matplotlib.backends.backend_qtagg import FigureCanvas
# Remember to import matplotlib after Qt.
from matplotlib.figure import Figure
from PyQt5.QtWidgets import (QApplication, QFrame, QHBoxLayout, QLabel,
                             QSizePolicy, QStackedWidget, QVBoxLayout, QWidget)
from qfluentwidgets import Action, BodyLabel, CommandBar, FlowLayout
from qfluentwidgets import FluentIcon
from qfluentwidgets import FluentIcon as FIF
from qfluentwidgets import (IconWidget, RoundMenu, SearchLineEdit,
                            SmoothScrollArea, StrongBodyLabel, Theme,
                            applyThemeColor, isDarkTheme)



class StarBar(QWidget):
    """ Star bar """

    def __init__(self, title, subtitle, parent=None):
        super().__init__(parent=parent)
        self.titleLabel = TitleLabel(title, self)
        #self.titleLabel = LineEdit(self)
        #self.titleLabel.setText(title)
        

        self.subtitleLabel = CaptionLabel(subtitle, self)


        self.saveButton = PushButton(
            self.tr('Save'), 
            self, 
            FluentIcon.SAVE
        )
        # add 'save as'
        self.saveAsButton = PushButton(
            self.tr('Save as'), 
            self, 
            FluentIcon.SAVE_AS
        )
        # add 'export'
        self.exportButton = PushButton(
            self.tr('Export tables'),
            self,
            FluentIcon.DOCUMENT
        )
        self.labelButton = PushButton(
            "Label",
            self,
            FluentIcon.TAG
        )
        

        self.simbadButton = PushButton(
            self.tr('SIMBAD'), 
            self, FluentIcon.SEARCH
        )
        self.simbadButton.setToolTip("Search SIMBAD for this object")
        self.gaiaButton = PushButton(            
            self.tr('Gaia'), 
            self,
            FluentIcon.SEARCH
        )
        self.gaiaButton.setToolTip("Search the Gaia data archive for these coordinates")
        
        self.vBoxLayout = QVBoxLayout(self)
        self.buttonLayout = QHBoxLayout()

        self.__initWidget()

    def __initWidget(self):
        self.setFixedHeight(138)
        self.vBoxLayout.setSpacing(0)
        self.vBoxLayout.setContentsMargins(36, 22, 36, 12)
        self.vBoxLayout.addWidget(self.titleLabel)
        self.vBoxLayout.addSpacing(4)
        self.vBoxLayout.addWidget(self.subtitleLabel)
        self.vBoxLayout.addSpacing(4)
        self.vBoxLayout.addLayout(self.buttonLayout, 1)
        self.vBoxLayout.setAlignment(Qt.AlignTop)

        self.buttonLayout.setSpacing(4)
        self.buttonLayout.setContentsMargins(0, 0, 0, 0)
        self.buttonLayout.addWidget(self.saveButton, 0, Qt.AlignLeft)
        self.buttonLayout.addWidget(self.saveAsButton, 0, Qt.AlignLeft)
        self.buttonLayout.addWidget(self.labelButton, 0, Qt.AlignLeft)
        self.buttonLayout.addWidget(self.exportButton, 0, Qt.AlignLeft)
        # add separator
        self.buttonLayout.addWidget(SeparatorWidget(), 0, Qt.AlignLeft)
        self.buttonLayout.addStretch(1)
                
        self.buttonLayout.addWidget(self.simbadButton, 0, Qt.AlignRight)
        self.buttonLayout.addWidget(self.gaiaButton, 0, Qt.AlignRight)
        self.buttonLayout.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)

        self.simbadButton.clicked.connect(
            lambda: QDesktopServices.openUrl(QUrl(HELP_URL)))
        self.gaiaButton.clicked.connect(
            lambda: QDesktopServices.openUrl(QUrl(EXAMPLE_URL)))

        self.subtitleLabel.setTextColor(QColor(96, 96, 96), QColor(216, 216, 216))




class PlotWidget(QWidget):
    
    def __init__(self, parent=None):
        super().__init__(parent=parent)
    
        fig = Figure(figsize=(5, 3))
        canvas = FigureCanvas(fig)

        ax = fig.add_subplot(111)
        ax.plot([0, 1, 2, 3, 4, 5], [10, 1, 20, 3, 40, 5])


        
        layout = QVBoxLayout()

        #layout.addWidget(toolbar)

        canvas_layout = QVBoxLayout()
        canvas_layout.addWidget(canvas)
        layout.addChildLayout(canvas_layout)
        

        #widget = QWidget()
        #widget.setLayout(layout)
        canvas.draw()
        
        self.setFixedHeight(200)


class NavigationToolbar(BaseNavigationToolbar):
    
    #def __init__(self, canvas, parent=None, coordinates=False):
    #    super().__init__(canvas, parent, coordinates)
    
    def __init__(self, canvas, parent=None, coordinates=False):
        """coordinates: should we show the coordinates on the right?"""
        QtWidgets.QToolBar.__init__(self, parent)
        #self.setAllowedAreas(QtCore.Qt.ToolBarArea(
        #    _to_int(QtCore.Qt.ToolBarArea.TopToolBarArea) |
        #    _to_int(QtCore.Qt.ToolBarArea.BottomToolBarArea)))
        self.coordinates = coordinates
        self._actions = {}  # mapping of toolitem method names to QActions.
        self._subplot_dialog = None

        for text, tooltip_text, icon, callback in self.toolitems:
            #if text is None:
            #    self.addSeparator()
            #else:
            a = self.addAction(icon.qicon(), text, getattr(self, callback))
            self._actions[callback] = a
            if callback in ['zoom', 'pan']:
                a.setCheckable(True)
            if tooltip_text is not None:
                a.setToolTip(tooltip_text)

        # Add the (x, y) location widget at the right side of the toolbar
        # The stretch factor is 1 which means any resizing of the toolbar
        # will resize this label instead of the buttons.

        NavigationToolbar2.__init__(self, canvas)
        
    @property
    def toolitems(self):
        return [
            ("Home", "Reset original view", FluentIcon.HOME, "home"),
            ('Pan', 'Left button pans, Right button zooms\nx/y fixes axis, CTRL fixes aspect', FluentIcon.MOVE, 'pan'),
            ('Zoom', 'Zoom to rectangle\nx/y fixes axis', FluentIcon.ZOOM_IN, 'zoom'),            
        ]









class ExampleCard(QWidget):
    """ Example card """

    def __init__(self, title, widget: QWidget, leftButtonText=None, rightButtonText=None, stretch=0, parent=None):
        super().__init__(parent=parent)
        self.widget = widget
        self.stretch = stretch

        self.titleLabel = StrongBodyLabel(title, self)
        self.card = QFrame(self)

        self.sourceWidget = QFrame(self.card)
        
        self.leftButtonText = leftButtonText
        self.rightButtonText = rightButtonText
        #self.sourcePathLabel = BodyLabel(
        #    self.tr('Source code'), 
        #    self.sourceWidget
        #)
        #self.linkIcon = IconWidget(FluentIcon.LINK, self.sourceWidget)

        self.vBoxLayout = QVBoxLayout(self)
        self.cardLayout = QVBoxLayout(self.card)
        self.topLayout = QHBoxLayout()
        self.bottomLayout = QHBoxLayout(self.sourceWidget)

        self.__initWidget()

    def __initWidget(self):
        #self.linkIcon.setFixedSize(16, 16)
        self.__initLayout()

        #self.sourceWidget.setCursor(Qt.PointingHandCursor)
        self.sourceWidget.installEventFilter(self)

        self.card.setObjectName('card')
        self.sourceWidget.setObjectName('sourceWidget')

    def __initLayout(self):
        self.vBoxLayout.setSizeConstraint(QVBoxLayout.SetMinimumSize)
        self.cardLayout.setSizeConstraint(QVBoxLayout.SetMinimumSize)
        self.topLayout.setSizeConstraint(QHBoxLayout.SetMinimumSize)

        self.vBoxLayout.setSpacing(12)
        self.vBoxLayout.setContentsMargins(0, 0, 0, 0)
        self.topLayout.setContentsMargins(12, 12, 12, 12)
        self.bottomLayout.setContentsMargins(18, 18, 18, 18)
        self.cardLayout.setContentsMargins(0, 0, 0, 0)

        self.vBoxLayout.addWidget(self.titleLabel, 0, Qt.AlignTop)
        self.vBoxLayout.addWidget(self.card, 0, Qt.AlignTop)
        self.vBoxLayout.setAlignment(Qt.AlignTop)

        self.cardLayout.setSpacing(0)
        self.cardLayout.setAlignment(Qt.AlignTop)
        self.cardLayout.addLayout(self.topLayout, 0)
        self.cardLayout.addWidget(self.sourceWidget, 0, Qt.AlignBottom)

        self.widget.setParent(self.card)
        self.topLayout.addWidget(self.widget)
        if self.stretch == 0:
            self.topLayout.addStretch(1)

        self.widget.show()

        #self.bottomLayout.addWidget(self.sourcePathLabel, 0, Qt.AlignLeft)
        #self.bottomLayout.addStretch(1)
        #self.bottomLayout.addWidget(self.linkIcon, 0, Qt.AlignRight)
        
        if self.leftButtonText is not None:
            leftButton = PushButton(self.leftButtonText)
            #leftButton.clicked.connect(self.enable_norm)
            self.bottomLayout.addWidget(leftButton, 0, Qt.AlignLeft)
            
        self.bottomLayout.addStretch(1)
        if self.rightButtonText is not None:
            rightButton = PrimaryPushButton(self.rightButtonText, self)
            #rightButton.clicked.connect(self.enable_norm)
            
            self.bottomLayout.addWidget(rightButton, 0, Qt.AlignRight)

        self.bottomLayout.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)

    #def enable_norm(self):
    #    self.parent.vBoxLayout.addWidget(self.parent.continuum_card)
                

    #def eventFilter(self, obj, e):
    #    if obj is self.sourceWidget:
    #        if e.type() == QEvent.MouseButtonRelease:
    #            QDesktopServices.openUrl(QUrl(self.sourcePath))

    #        return super().eventFilter(obj, e)
    
#from style_sheet import StyleSheet


class SinglePlotWidget(QWidget):
    
    def __init__(
        self, 
        parent=None, 
        figsize=(8, 2),
        size_policy=(QSizePolicy.Expanding, QSizePolicy.Fixed),
        toolbar=True,
    ):
        super().__init__(parent)
        self.parent = parent
        
        self.figure = Figure(figsize=figsize)
        canvas = FigureCanvas(self.figure)
        canvas.setSizePolicy(*size_policy)
        
        self.ax = canvas.figure.subplots()
        self.ax.plot([1, 2, 3, 4, 5], [1, 2, 3, 4, 5])
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.figure.tight_layout()
        self.figure.canvas.draw()
        
        self.layout = QVBoxLayout(self)
        self.layout.addWidget(canvas)
        
        #StyleSheet.GALLERY_INTERFACE.apply(self)
        
        if toolbar:
    
            
            toolbar = CommandBar(parent=self)
            
            action_home = Action(FluentIcon.HOME, "Home", self)
            action_zoom = Action(FluentIcon.ZOOM_IN, "Zoom", self)
            action_zoom.setCheckable(True)
            action_pan = Action(FluentIcon.MOVE, "Pan", self)
            action_pan.setCheckable(True)
            action_home.triggered.connect(lambda: print("home"))
            action_zoom.triggered.connect(lambda: print("zoom"))
            action_pan.triggered.connect(lambda: print("pan"))
                        
            action_left = Action(FluentIcon.PAGE_LEFT, "Left", self)
            action_left.setDisabled(True)
            action_right = Action(FluentIcon.PAGE_RIGHT, "Right", self)
            action_left.triggered.connect(lambda: print("left triggered"))
            action_left.hovered.connect(lambda: print("left hovered"))
            
            toolbar.addAction(action_home)
            toolbar.addAction(action_zoom)
            toolbar.addAction(action_pan)
            toolbar.addSeparator()
            toolbar.addAction(action_left)
            toolbar.addAction(action_right)
            self.layout.addWidget(toolbar)
        
        self.layout.addStretch(1)    
        return None
    
    


class MySubAnalysisInterface(ScrollArea):
    """ Dialog interface """

    def __init__(self, text, parent=None):
        super().__init__(parent=parent)
        self.view = QWidget(self)
        self.toolBar = StarBar(text, "15:43:03 -10:56:01", self)
        self.vBoxLayout = QVBoxLayout(self.view)

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setViewportMargins(0, self.toolBar.height(), 0, 0)
        self.setWidget(self.view)
        self.setWidgetResizable(True)

        self.vBoxLayout.setSpacing(30)
        self.vBoxLayout.setAlignment(Qt.AlignTop)
        self.vBoxLayout.setContentsMargins(36, 20, 36, 36)

        self.view.setObjectName('view')
        StyleSheet.GALLERY_INTERFACE.apply(self)
        
        self.setObjectName('fooInterface')

        button1 = PushButton(self.tr('Show dialog'))
        button1.clicked.connect(self.showDialog)


        # Radial velocity stuff
        
        
        rv_view = QWidget(self)
        #rv_view.setFixedHeight(200)
        # set rv_view to span to fill available horizontal space
        rv_layout = QHBoxLayout(rv_view)
        rv_layout.setContentsMargins(0, 0, 0, 0)
                
        plot_layout = QVBoxLayout()
        #plot_layout.setSpacing(30)
        
        plt.style.use({
            "font.family": "sans-serif",
            #"font.sans-serif": ["Segoe UI"],
            "axes.labelsize": "10",
            "lines.linewidth": "2",
            "lines.markersize": "10",
            "xtick.labelsize": "10",
            "ytick.labelsize": "10",

            #"figure.facecolor": "#272727f2", #-> main background color
            #"figure.facecolor": "#202020",
            "axes.facecolor": "#292929", # rgba(41,41,41,255)            
            "figure.facecolor": "#292929", # rgba(41,41,41,255)
            "axes.edgecolor": "white",
            "axes.labelcolor": "white",
            "xtick.color": "white",
            "xtick.labelcolor": "white",
            "ytick.color": "white",
            "ytick.labelcolor": "white"
        })        
        
        plot_layout.addWidget(SinglePlotWidget())
        #plot_layout.addWidget(NavigationToolbar(dynamic_canvas, self)) # Too hard to get looking beautiful right now

        t = np.linspace(0, 10, 501)
        
        rv_layout.addLayout(plot_layout)
        
        button_select_template = PushButton("Select template")
        rv_nav_layout = QVBoxLayout()
        rv_nav_layout.addWidget(button_select_template)
        # align to left:
        rv_nav_layout.addWidget(QLabel("Wavelength:"), 0, Qt.AlignLeft)
        
        self.comboBox = EditableComboBox(self)
        self.comboBox.setPlaceholderText(cfg.get(cfg.RVWavelengthRange))
        
        #items = ['shoko 🥰', '西宫硝子', '宝多六花', '小鸟游六花']
        #self.comboBox.addItems(items)
        self.comboBox.setCurrentIndex(-1)
        self.comboBox.currentTextChanged.connect(print)
        self.completer = QCompleter([], self)
        self.comboBox.setCompleter(self.completer)
                
        rv_nav_layout.addWidget(self.comboBox)
        rv_nav_layout.addStretch(1)
        
        #rv_nav_layout.setFixedWidth(100)
        rv_layout.addLayout(rv_nav_layout)
        
        '''
        self.addExampleCard(
            self.tr("Radial Velocity"),
            rv_view,
            "Do not shift",
            "Measure and shift"
        )
        '''

        card = ExampleCard("Radial velocity", rv_view, "Do not shift", "Measure and shift", self.view)
        self.vBoxLayout.addWidget(card, 0, Qt.AlignTop)
        
        # Continuum normalization
        
        continuum_view = QWidget(self)
        continuum_layout = QHBoxLayout(continuum_view)
        continuum_layout.setContentsMargins(0, 0, 0, 0)

        continuum_layout.addWidget(SinglePlotWidget())
        self.continuum_card = ExampleCard(
            self.tr("Continuum normalization"),
            continuum_view,
            "Do not normalize",
            "Normalize"
        )
        
        
        '''
        self.addExampleCard(
            self.tr('A frameless message box'),
            button,
            'https://github.com/zhiyiYo/PyQt-Fluent-Widgets/blob/master/examples/dialog_flyout/dialog/demo.py'
        )
        
        button = PushButton(self.tr('Show dialog'))
        button.clicked.connect(self.showDialog)
        self.addExampleCard(
            self.tr('A frameless message box'),
            button,
            'https://github.com/zhiyiYo/PyQt-Fluent-Widgets/blob/master/examples/dialog_flyout/dialog/demo.py'
        )
        
        button = PushButton(self.tr('Show dialog'))
        button.clicked.connect(self.showDialog)
        self.addExampleCard(
            self.tr('A frameless message box'),
            button,
            'https://github.com/zhiyiYo/PyQt-Fluent-Widgets/blob/master/examples/dialog_flyout/dialog/demo.py'
        )
        
                        
        button = PushButton(self.tr('Show dialog'))
        button.clicked.connect(self.showDialog)
        self.addExampleCard(
            self.tr('A frameless message box'),
            button,
            'https://github.com/zhiyiYo/PyQt-Fluent-Widgets/blob/master/examples/dialog_flyout/dialog/demo.py'
        )
        
                
        button = PushButton(self.tr('Show dialog'))
        button.clicked.connect(self.showDialog)
        self.addExampleCard(
            self.tr('A frameless message box'),
            button,
            'https://github.com/zhiyiYo/PyQt-Fluent-Widgets/blob/master/examples/dialog_flyout/dialog/demo.py'
        )
        
                
        button = PushButton(self.tr('Show dialog'))
        button.clicked.connect(self.showDialog)
        self.addExampleCard(
            self.tr('A frameless message box'),
            button,
            'https://github.com/zhiyiYo/PyQt-Fluent-Widgets/blob/master/examples/dialog_flyout/dialog/demo.py'
        )
        '''
        
    def _update_canvas(self):
        t = np.linspace(0, 10, 101)
        # Shift the sinusoid as a function of time.
        self._line.set_data(t, np.sin(t + time.time()))
        self._line.figure.tight_layout()
        self._line.figure.canvas.draw()

                                                        

    def addExampleCard(self, title, widget, sourcePath: str, stretch=0):
        card = ExampleCard(title, widget, sourcePath, stretch, self.view)
        self.vBoxLayout.addWidget(card, 0, Qt.AlignTop)
        return card
        
    def showDialog(self):
        print("NAH")


class TabInterface(QWidget):
    """ Tab interface """

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.tabCount = 0

        self.tabBar = TabBar(self)
        self.tabBar.setMovable(True)
        self.stackedWidget = QStackedWidget(self)
        
        self.tabView = QWidget(self)

        self.hBoxLayout = QHBoxLayout(self)
        self.vBoxLayout = QVBoxLayout(self.tabView)

        # add items to pivot
        self.tabBar.setTabMaximumWidth(200)

        self.hBoxLayout.addWidget(self.tabView, 1)
        self.hBoxLayout.setContentsMargins(0, 0, 0, 0)

        self.vBoxLayout.addWidget(self.tabBar)
        self.vBoxLayout.addWidget(self.stackedWidget)
        self.vBoxLayout.setContentsMargins(0, 0, 0, 0)
        StyleSheet.GALLERY_INTERFACE.apply(self)

        self.connectSignalToSlot()

        self.tabBar.setCloseButtonDisplayMode(TabCloseButtonDisplayMode.ON_HOVER)

        self.myInterface = MySubAnalysisInterface("HD 140283", self)        
        self.addMySubInterface(self.myInterface, 'myInterface', 'HD140283')

        qrouter.setDefaultRouteKey(
            self.stackedWidget, self.myInterface.objectName())
        
        
        
    def connectSignalToSlot(self):
        self.tabBar.tabAddRequested.connect(self.addTab)
        self.tabBar.tabCloseRequested.connect(self.removeTab)

        self.stackedWidget.currentChanged.connect(self.onCurrentIndexChanged)


    def addSubInterface(self, widget: QLabel, objectName, text, icon):
        widget.setObjectName(objectName)
        widget.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        self.stackedWidget.addWidget(widget)
        self.tabBar.addTab(
            routeKey=objectName,
            text=text,
            icon=icon,
            onClick=lambda: self.stackedWidget.setCurrentWidget(widget)
        )

    def addMySubInterface(self, thing, objectName, text):
        thing.setObjectName(objectName)
        thing.setAlignment(Qt.AlignTop | Qt.AlignLeft)

        self.stackedWidget.addWidget(thing)
        self.tabBar.addTab(
            routeKey=objectName,
            text=text,
            icon=None,
            onClick=lambda: self.stackedWidget.setCurrentWidget(thing)            
        )


    def onCurrentIndexChanged(self, index):
        widget = self.stackedWidget.widget(index)
        if not widget:
            return

        self.tabBar.setCurrentTab(widget.objectName())
        qrouter.push(self.stackedWidget, widget.objectName())


    def addTab(self):
        # load new file to select
        filenames, selected_filter = QFileDialog.getOpenFileNames(
            self, 
            caption="Select input spectra", 
            directory="", 
            filter="*"
        )
        if filenames:            
            
            # load the file, parse a name from it.        
            text = f"Untitled-{self.tabCount}"
            widget = MySubAnalysisInterface(text, self)

            self.addMySubInterface(widget, text, text)        
            self.tabCount += 1
            
            # now set the new sub interface widget as the current one
            #print(self.tabBar.setCurrentTab(widget.objectName()))
            #print(self.tabBar.setCurrentIndex(self.tabCount - 1))
            self.stackedWidget.setCurrentIndex(self.tabCount - 1)
            
            
            
            
            #self.tabBar.setCurrentTab(widget.objectName())
            #qrouter.push(self.stackedWidget, widget.objectName())
            #self.onCurrentIndexChanged(self.tabBar.currentIndex())
            


    def removeTab(self, index):
        # ask to save before quit
        item = self.tabBar.tabItem(index)
        print(index, item, item.routeKey)
        widget = self.findChild(QLabel, item.routeKey())

        self.stackedWidget.removeWidget(widget)
        self.tabBar.removeTab(index)
        self.tabCount -= 1
        try:
            widget.deleteLater()
        except:
            print("cannot delete later")


class AnalysisInterface(ScrollArea):
    """ Analysis Interface """

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.view = TabInterface(self)        
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setViewportMargins(0, 0, 0, 0)
        self.setWidget(self.view)
        self.setWidgetResizable(True)
        self.view.setObjectName('view')
        #StyleSheet.GALLERY_INTERFACE.apply(self)
        self.setObjectName('AnalysisInterface')

