# coding:utf-8
from typing import List

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import QApplication, QFrame, QVBoxLayout, QLabel, QWidget, QHBoxLayout
from qfluentwidgets import (FluentIcon, IconWidget, FlowLayout, isDarkTheme,
                            Theme, applyThemeColor, SmoothScrollArea, SearchLineEdit, StrongBodyLabel,
                            BodyLabel)

from .gallery_interface import GalleryInterface, ToolBar, SeparatorWidget, ExampleCard
from ..common.translator import Translator
from ..common.config import cfg
from ..common.style_sheet import StyleSheet
from ..common.trie import Trie


class LineEdit(SearchLineEdit):
    """ Search line edit """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setPlaceholderText(self.tr('Search icons'))
        self.setFixedWidth(304)
        self.textChanged.connect(self.search)


class IconCard(QFrame):
    """ Icon card """

    clicked = pyqtSignal(FluentIcon)

    def __init__(self, icon: FluentIcon, parent=None):
        super().__init__(parent=parent)
        self.icon = icon
        self.isSelected = False

        self.iconWidget = IconWidget(icon, self)
        self.nameLabel = QLabel(self)
        self.vBoxLayout = QVBoxLayout(self)

        self.setFixedSize(96, 96)
        self.vBoxLayout.setSpacing(0)
        self.vBoxLayout.setContentsMargins(8, 28, 8, 0)
        self.vBoxLayout.setAlignment(Qt.AlignTop)
        self.iconWidget.setFixedSize(28, 28)
        self.vBoxLayout.addWidget(self.iconWidget, 0, Qt.AlignHCenter)
        self.vBoxLayout.addSpacing(14)
        self.vBoxLayout.addWidget(self.nameLabel, 0, Qt.AlignHCenter)

        text = self.nameLabel.fontMetrics().elidedText(icon.value, Qt.ElideRight, 90)
        self.nameLabel.setText(text)

    def mouseReleaseEvent(self, e):
        if self.isSelected:
            return

        self.clicked.emit(self.icon)

    def setSelected(self, isSelected: bool, force=False):
        if isSelected == self.isSelected and not force:
            return

        self.isSelected = isSelected

        if not isSelected:
            self.iconWidget.setIcon(self.icon)
        else:
            icon = self.icon.icon(Theme.LIGHT if isDarkTheme() else Theme.DARK)
            self.iconWidget.setIcon(icon)

        self.setProperty('isSelected', isSelected)
        self.setStyle(QApplication.style())


class IconInfoPanel(QFrame):
    """ Icon info panel """

    def __init__(self, icon: FluentIcon, parent=None):
        super().__init__(parent=parent)
        self.nameLabel = QLabel(icon.value, self)
        self.iconWidget = IconWidget(icon, self)
        self.iconNameTitleLabel = QLabel(self.tr('Icon name'), self)
        self.iconNameLabel = QLabel(icon.value, self)
        self.enumNameTitleLabel = QLabel(self.tr('Enum member'), self)
        self.enumNameLabel = QLabel("FluentIcon." + icon.name, self)

        self.vBoxLayout = QVBoxLayout(self)
        self.vBoxLayout.setContentsMargins(16, 20, 16, 20)
        self.vBoxLayout.setSpacing(0)
        self.vBoxLayout.setAlignment(Qt.AlignTop)

        self.vBoxLayout.addWidget(self.nameLabel)
        self.vBoxLayout.addSpacing(16)
        self.vBoxLayout.addWidget(self.iconWidget)
        self.vBoxLayout.addSpacing(45)
        self.vBoxLayout.addWidget(self.iconNameTitleLabel)
        self.vBoxLayout.addSpacing(5)
        self.vBoxLayout.addWidget(self.iconNameLabel)
        self.vBoxLayout.addSpacing(34)
        self.vBoxLayout.addWidget(self.enumNameTitleLabel)
        self.vBoxLayout.addSpacing(5)
        self.vBoxLayout.addWidget(self.enumNameLabel)

        self.iconWidget.setFixedSize(48, 48)
        self.setFixedWidth(216)

        self.nameLabel.setObjectName('nameLabel')
        self.iconNameTitleLabel.setObjectName('subTitleLabel')
        self.enumNameTitleLabel.setObjectName('subTitleLabel')

    def setIcon(self, icon: FluentIcon):
        self.iconWidget.setIcon(icon)
        self.nameLabel.setText(icon.value)
        self.iconNameLabel.setText(icon.value)
        self.enumNameLabel.setText("FluentIcon."+icon.name)


class IconCardView(QWidget):
    """ Icon card view """

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.trie = Trie()
        self.iconLibraryLabel = StrongBodyLabel(self.tr('Fluent Icons Library'), self)
        self.searchLineEdit = LineEdit(self)

        self.view = QFrame(self)
        self.scrollArea = SmoothScrollArea(self.view)
        self.scrollWidget = QWidget(self.scrollArea)
        self.infoPanel = IconInfoPanel(FluentIcon.MENU, self)

        self.vBoxLayout = QVBoxLayout(self)
        self.hBoxLayout = QHBoxLayout(self.view)
        self.flowLayout = FlowLayout(self.scrollWidget, isTight=True)

        self.cards = []     # type:List[IconCard]
        self.icons = []
        self.currentIndex = -1

        self.__initWidget()

    def __initWidget(self):
        self.scrollArea.setWidget(self.scrollWidget)
        self.scrollArea.setViewportMargins(0, 5, 0, 5)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.vBoxLayout.setContentsMargins(0, 0, 0, 0)
        self.vBoxLayout.setSpacing(12)
        self.vBoxLayout.addWidget(self.iconLibraryLabel)
        self.vBoxLayout.addWidget(self.searchLineEdit)
        self.vBoxLayout.addWidget(self.view)

        self.hBoxLayout.setSpacing(0)
        self.hBoxLayout.setContentsMargins(0, 0, 0, 0)
        self.hBoxLayout.addWidget(self.scrollArea)
        self.hBoxLayout.addWidget(self.infoPanel, 0, Qt.AlignRight)

        self.flowLayout.setVerticalSpacing(8)
        self.flowLayout.setHorizontalSpacing(8)
        self.flowLayout.setContentsMargins(8, 3, 8, 8)

        self.__setQss()
        cfg.themeChanged.connect(self.__setQss)
        self.searchLineEdit.clearSignal.connect(self.showAllIcons)
        self.searchLineEdit.searchSignal.connect(self.search)

        for icon in FluentIcon._member_map_.values():
            self.addIcon(icon)

        self.setSelectedIcon(self.icons[0])

    def addIcon(self, icon: FluentIcon):
        """ add icon to view """
        card = IconCard(icon, self)
        card.clicked.connect(self.setSelectedIcon)

        self.trie.insert(icon.value, len(self.cards))
        self.cards.append(card)
        self.icons.append(icon)
        self.flowLayout.addWidget(card)

    def setSelectedIcon(self, icon: FluentIcon):
        """ set selected icon """
        index = self.icons.index(icon)

        if self.currentIndex >= 0:
            self.cards[self.currentIndex].setSelected(False)

        self.currentIndex = index
        self.cards[index].setSelected(True)
        self.infoPanel.setIcon(icon)

    def __setQss(self):
        self.view.setObjectName('iconView')
        self.scrollWidget.setObjectName('scrollWidget')

        StyleSheet.ICON_INTERFACE.apply(self)
        StyleSheet.ICON_INTERFACE.apply(self.scrollWidget)

        if self.currentIndex >= 0:
            self.cards[self.currentIndex].setSelected(True, True)

    def search(self, keyWord: str):
        """ search icons """
        items = self.trie.items(keyWord.lower())
        indexes = {i[1] for i in items}
        self.flowLayout.removeAllWidgets()

        for i, card in enumerate(self.cards):
            isVisible = i in indexes
            card.setVisible(isVisible)
            if isVisible:
                self.flowLayout.addWidget(card)

    def showAllIcons(self):
        self.flowLayout.removeAllWidgets()
        for card in self.cards:
            card.show()
            self.flowLayout.addWidget(card)


from PyQt5.QtWidgets import QWidget, QStackedWidget, QVBoxLayout, QLabel, QHBoxLayout, QFrame, QSizePolicy
from qfluentwidgets import (Pivot, qrouter, SegmentedWidget, TabBar, CheckBox, ComboBox, EditableComboBox,
                            TabCloseButtonDisplayMode, BodyLabel, SpinBox, BreadcrumbBar,
                            SegmentedToggleToolWidget, FluentIcon)

# coding:utf-8
from PyQt5.QtCore import Qt, pyqtSignal, QUrl, QEvent
from PyQt5.QtGui import QDesktopServices, QPainter, QPen, QColor, QIcon
from PyQt5.QtWidgets import QWidget, QLabel, QVBoxLayout, QHBoxLayout, QFrame

from qfluentwidgets import (ScrollArea, PushButton, ToolButton, FluentIcon,
                            isDarkTheme, IconWidget, Theme, ToolTipFilter, TitleLabel, CaptionLabel,
                            StrongBodyLabel, BodyLabel, toggleTheme)
from ..common.config import cfg, FEEDBACK_URL, HELP_URL, EXAMPLE_URL
from ..common.icon import Icon
from ..common.style_sheet import StyleSheet
from ..common.signal_bus import signalBus

class StarBar(QWidget):
    """ Tool bar """

    def __init__(self, title, subtitle, parent=None):
        super().__init__(parent=parent)
        self.titleLabel = TitleLabel(title, self)
        self.subtitleLabel = CaptionLabel(subtitle, self)

        self.simbadButton = PushButton(
            self.tr('SIMBAD'), 
            self, FluentIcon.SEARCH
        )
        self.simbadButton.setToolTip("Search SIMBAD for this object")
        self.gaiaButton = PushButton(self.tr('Gaia archive'), self, None)
        self.gaiaButton.setToolTip("Search the Gaia data archive for these coordinates")
        
        #self.themeButton = ToolButton(FluentIcon.CONSTRACT, self)
        self.separator = SeparatorWidget(self)
        self.supportButton = ToolButton(FluentIcon.HEART, self)
        self.feedbackButton = ToolButton(FluentIcon.FEEDBACK, self)

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
        self.buttonLayout.addWidget(self.simbadButton, 0, Qt.AlignLeft)
        self.buttonLayout.addWidget(self.gaiaButton, 0, Qt.AlignLeft)
        self.buttonLayout.addStretch(1)
        #self.buttonLayout.addWidget(self.themeButton, 0, Qt.AlignRight)
        self.buttonLayout.addWidget(self.separator, 0, Qt.AlignRight)
        self.buttonLayout.addWidget(self.supportButton, 0, Qt.AlignRight)
        self.buttonLayout.addWidget(self.feedbackButton, 0, Qt.AlignRight)
        self.buttonLayout.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)

        #self.themeButton.installEventFilter(ToolTipFilter(self.themeButton))
        self.supportButton.installEventFilter(ToolTipFilter(self.supportButton))
        self.feedbackButton.installEventFilter(
            ToolTipFilter(self.feedbackButton))
        #self.themeButton.setToolTip(self.tr('Toggle theme'))
        self.supportButton.setToolTip(self.tr('Support me'))
        self.feedbackButton.setToolTip(self.tr('Send feedback'))

        #self.themeButton.clicked.connect(lambda: toggleTheme(True))
        self.supportButton.clicked.connect(signalBus.supportSignal)
        self.simbadButton.clicked.connect(
            lambda: QDesktopServices.openUrl(QUrl(HELP_URL)))
        self.gaiaButton.clicked.connect(
            lambda: QDesktopServices.openUrl(QUrl(EXAMPLE_URL)))
        self.feedbackButton.clicked.connect(
            lambda: QDesktopServices.openUrl(QUrl(FEEDBACK_URL)))

        self.subtitleLabel.setTextColor(QColor(96, 96, 96), QColor(216, 216, 216))


from qfluentwidgets import IconWidget, TextWrap, FlowLayout, CardWidget


from .dialog_interface import DialogInterface    

from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import \
    NavigationToolbar2QT as BaseNavigationToolbar
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.figure import Figure

import time

import numpy as np

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

import matplotlib.pyplot as plt
from matplotlib.backend_bases import NavigationToolbar2

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






from PyQt5.QtWidgets import QWidget, QStackedWidget, QVBoxLayout, QLabel, QHBoxLayout, QFrame, QSizePolicy

from PyQt5.QtWidgets import QApplication, QFrame, QVBoxLayout, QLabel, QWidget, QHBoxLayout
from qfluentwidgets import (FluentIcon, IconWidget, FlowLayout, isDarkTheme,
                            Theme, applyThemeColor, SmoothScrollArea, SearchLineEdit, StrongBodyLabel,
                            BodyLabel, CommandBar)
from qfluentwidgets import RoundMenu, Action
from qfluentwidgets import FluentIcon as FIF

# Remember to import matplotlib after Qt.
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvas
#from matplotlib.backends.backend_qtagg import \
#    NavigationToolbar2QT as BaseNavigationToolbar
#from matplotlib.backends.qt_compat import QtWidgets

class ExampleCard(QWidget):
    """ Example card """

    def __init__(self, title, widget: QWidget, sourcePath, stretch=0, parent=None):
        super().__init__(parent=parent)
        self.widget = widget
        self.stretch = stretch

        self.titleLabel = StrongBodyLabel(title, self)
        self.card = QFrame(self)

        self.sourceWidget = QFrame(self.card)
        self.sourcePath = sourcePath
        self.sourcePathLabel = BodyLabel(
            self.tr('Source code'), self.sourceWidget)
        self.linkIcon = IconWidget(FluentIcon.LINK, self.sourceWidget)

        self.vBoxLayout = QVBoxLayout(self)
        self.cardLayout = QVBoxLayout(self.card)
        self.topLayout = QHBoxLayout()
        self.bottomLayout = QHBoxLayout(self.sourceWidget)

        self.__initWidget()

    def __initWidget(self):
        self.linkIcon.setFixedSize(16, 16)
        self.__initLayout()

        self.sourceWidget.setCursor(Qt.PointingHandCursor)
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

        self.bottomLayout.addWidget(self.sourcePathLabel, 0, Qt.AlignLeft)
        self.bottomLayout.addStretch(1)
        self.bottomLayout.addWidget(self.linkIcon, 0, Qt.AlignRight)
        self.bottomLayout.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)

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

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.view = QWidget(self)
        self.toolBar = StarBar("HD140283", "subtitle", self)
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
        
        fig = Figure(figsize=(8, 2))
        dynamic_canvas = FigureCanvas(fig)
        dynamic_canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)

        dynamic_canvas.setStyle(QApplication.style())
        
        plot_layout.addWidget(dynamic_canvas)
        plot_layout.addWidget(SinglePlotWidget())
        #plot_layout.addWidget(NavigationToolbar(dynamic_canvas, self)) # Too hard to get looking beautiful right now

        t = np.linspace(0, 10, 501)

        self._dynamic_ax = dynamic_canvas.figure.subplots()
        self._dynamic_ax.set_xlabel("Wavelength [$\AA$]")
        t = np.linspace(0, 10, 101)
        # Set up a Line2D.
        self._line, = self._dynamic_ax.plot(t, np.sin(t + time.time()))
        self._timer = dynamic_canvas.new_timer(50)
        self._timer.add_callback(self._update_canvas)
        self._timer.start()
        
        rv_layout.addLayout(plot_layout)
        
        button_select_template = PushButton("Select template")
        rv_nav_layout = QVBoxLayout()
        rv_nav_layout.addWidget(button_select_template)
        # align to left:
        rv_nav_layout.addWidget(QLabel("Wavelength:"), 0, Qt.AlignLeft)
        
        self.comboBox = EditableComboBox(self)
        self.comboBox.setPlaceholderText("选择一个脑婆")

        items = ['shoko 🥰', '西宫硝子', '宝多六花', '小鸟游六花']
        self.comboBox.addItems(items)
        self.comboBox.setCurrentIndex(-1)
        self.comboBox.currentTextChanged.connect(print)
        self.completer = QCompleter(items, self)
        self.comboBox.setCompleter(self.completer)
                
        rv_nav_layout.addWidget(self.comboBox)
        rv_nav_layout.addStretch(1)
        
        #rv_nav_layout.setFixedWidth(100)
        rv_layout.addLayout(rv_nav_layout)
        
        
        self.addExampleCard(
            self.tr("Radial Velocity"),
            rv_view,
            ""
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
        
                
        button = PushButton(self.tr('Show dialog'))
        button.clicked.connect(self.showDialog)
        self.addExampleCard(
            self.tr('A frameless message box'),
            button,
            'https://github.com/zhiyiYo/PyQt-Fluent-Widgets/blob/master/examples/dialog_flyout/dialog/demo.py'
        )
        
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
        self.tabCount = 1

        self.tabBar = TabBar(self)
        self.stackedWidget = QStackedWidget(self)
        
        self.tabView = QWidget(self)
        #self.controlPanel = QFrame(self)

        #self.movableCheckBox = CheckBox(self.tr('IsTabMovable'), self)
        #self.scrollableCheckBox = CheckBox(self.tr('IsTabScrollable'), self)
        #self.shadowEnabledCheckBox = CheckBox(self.tr('IsTabShadowEnabled'), self)
        #self.tabMaxWidthLabel = BodyLabel(self.tr('TabMaximumWidth'), self)
        #self.tabMaxWidthSpinBox = SpinBox(self)
        #self.closeDisplayModeLabel = BodyLabel(self.tr('TabCloseButtonDisplayMode'), self)
        #self.closeDisplayModeComboBox = ComboBox(self)

        self.hBoxLayout = QHBoxLayout(self)
        self.vBoxLayout = QVBoxLayout(self.tabView)
        #self.panelLayout = QVBoxLayout(self.controlPanel)

        #self.songInterface = QLabel('Song Interface', self)
        #self.myInterface = StarBar("HD140283", "15:43:03.09 -10:56:00.59", self)
        self.myInterface = MySubAnalysisInterface(self)
        self.diagInterface = DialogInterface(self)

        # add items to pivot
        self.__initWidget()

    def __initWidget(self):
        self.initLayout()

        #self.shadowEnabledCheckBox.setChecked(True)

        #self.tabMaxWidthSpinBox.setRange(60, 400)
        #self.tabMaxWidthSpinBox.setValue(self.tabBar.tabMaximumWidth())

        #self.closeDisplayModeComboBox.addItem(self.tr('Always'), userData=TabCloseButtonDisplayMode.ALWAYS)
        #self.closeDisplayModeComboBox.addItem(self.tr('OnHover'), userData=TabCloseButtonDisplayMode.ON_HOVER)
        #self.closeDisplayModeComboBox.addItem(self.tr('Never'), userData=TabCloseButtonDisplayMode.NEVER)
        #self.closeDisplayModeComboBox.currentIndexChanged.connect(self.onDisplayModeChanged)

        #self.addSubInterface(self.songInterface,
        #                     'tabSongInterface', self.tr('Song'), ':/gallery/images/MusicNote.png')
        
        self.addMySubInterface(self.myInterface, 'myInterface', 'HD140283')
        self.addMySubInterface(self.diagInterface, "diagInterface", "HD122563")

        #self.controlPanel.setObjectName('controlPanel')
        StyleSheet.NAVIGATION_VIEW_INTERFACE.apply(self)

        #self.connectSignalToSlot()

        self.tabBar.setCloseButtonDisplayMode(TabCloseButtonDisplayMode.ON_HOVER)

        qrouter.setDefaultRouteKey(
            self.stackedWidget, self.myInterface.objectName())
        
        
    '''
    def connectSignalToSlot(self):
        self.movableCheckBox.stateChanged.connect(
            lambda: self.tabBar.setMovable(self.movableCheckBox.isChecked()))
        self.scrollableCheckBox.stateChanged.connect(
            lambda: self.tabBar.setScrollable(self.scrollableCheckBox.isChecked()))
        self.shadowEnabledCheckBox.stateChanged.connect(
            lambda: self.tabBar.setTabShadowEnabled(self.shadowEnabledCheckBox.isChecked()))

        self.tabMaxWidthSpinBox.valueChanged.connect(self.tabBar.setTabMaximumWidth)

        self.tabBar.tabAddRequested.connect(self.addTab)
        self.tabBar.tabCloseRequested.connect(self.removeTab)

        self.stackedWidget.currentChanged.connect(self.onCurrentIndexChanged)
    '''
    def initLayout(self):
        self.tabBar.setTabMaximumWidth(200)

        #self.controlPanel.setFixedWidth(220)
        self.hBoxLayout.addWidget(self.tabView, 1)
        #self.hBoxLayout.addWidget(self.controlPanel, 0, Qt.AlignRight)
        self.hBoxLayout.setContentsMargins(0, 0, 0, 0)

        self.vBoxLayout.addWidget(self.tabBar)
        self.vBoxLayout.addWidget(self.stackedWidget)
        self.vBoxLayout.setContentsMargins(0, 0, 0, 0)

        #self.panelLayout.setSpacing(8)
        #self.panelLayout.setContentsMargins(14, 16, 14, 14)
        #self.panelLayout.setAlignment(Qt.AlignTop)

        #self.panelLayout.addWidget(self.movableCheckBox)
        #self.panelLayout.addWidget(self.scrollableCheckBox)
        #self.panelLayout.addWidget(self.shadowEnabledCheckBox)

        #self.panelLayout.addSpacing(4)
        #self.panelLayout.addWidget(self.tabMaxWidthLabel)
        #self.panelLayout.addWidget(self.tabMaxWidthSpinBox)

        #self.panelLayout.addSpacing(4)
        #self.panelLayout.addWidget(self.closeDisplayModeLabel)
        #self.panelLayout.addWidget(self.closeDisplayModeComboBox)

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


    def onDisplayModeChanged(self, index):
        mode = self.closeDisplayModeComboBox.itemData(index)
        self.tabBar.setCloseButtonDisplayMode(mode)

    def onCurrentIndexChanged(self, index):
        widget = self.stackedWidget.widget(index)
        if not widget:
            return

        self.tabBar.setCurrentTab(widget.objectName())
        qrouter.push(self.stackedWidget, widget.objectName())

    def addTab(self):
        # load new file to select
        text = f'硝子酱一级棒卡哇伊×{self.tabCount}'
        self.addSubInterface(QLabel('🥰 ' + text), text, text, ':/gallery/images/Smiling_with_heart.png')
        self.tabCount += 1

    def removeTab(self, index):
        # ask to save before quit
        item = self.tabBar.tabItem(index)
        widget = self.findChild(QLabel, item.routeKey())

        self.stackedWidget.removeWidget(widget)
        self.tabBar.removeTab(index)
        widget.deleteLater()


# coding:utf-8
from PyQt5.QtCore import Qt, pyqtSignal, QUrl, QEvent
from PyQt5.QtGui import QDesktopServices, QPainter, QPen, QColor
from PyQt5.QtWidgets import QWidget, QLabel, QVBoxLayout, QHBoxLayout, QFrame, QCompleter

from qfluentwidgets import (ScrollArea, PushButton, ToolButton, FluentIcon,
                            isDarkTheme, IconWidget, Theme, ToolTipFilter, TitleLabel, CaptionLabel,
                            StrongBodyLabel, BodyLabel, toggleTheme)
from ..common.config import cfg, FEEDBACK_URL, HELP_URL, EXAMPLE_URL
from ..common.icon import Icon
from ..common.style_sheet import StyleSheet
from ..common.signal_bus import signalBus

class GalleryInterface(ScrollArea):
    """ Gallery interface """

    def __init__(self, title: str, subtitle: str, parent=None):
        """
        Parameters
        ----------
        title: str
            The title of gallery

        subtitle: str
            The subtitle of gallery

        parent: QWidget
            parent widget
        """
        super().__init__(parent=parent)
        self.view = QWidget(self)
        #self.toolBar = ToolBar(title, subtitle, self)
        self.vBoxLayout = QVBoxLayout(self.view)

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setViewportMargins(0, 0, 0, 0)
        self.setWidget(self.view)
        self.setWidgetResizable(True)

        #self.vBoxLayout.setSpacing(30)
        self.vBoxLayout.setAlignment(Qt.AlignTop)
        #self.vBoxLayout.setContentsMargins(36, 20, 36, 36)

        self.view.setObjectName('view')
        StyleSheet.GALLERY_INTERFACE.apply(self)

    def addExampleCard(self, title, widget, sourcePath: str, stretch=0):
        card = ExampleCard(title, widget, sourcePath, stretch, self.view)
        self.vBoxLayout.addWidget(card, 0, Qt.AlignTop)
        return card

    def scrollToCard(self, index: int):
        """ scroll to example card """
        w = self.vBoxLayout.itemAt(index).widget()
        self.verticalScrollBar().setValue(w.y())

    #def resizeEvent(self, e):
    #    super().resizeEvent(e)
    #    self.toolBar.resize(self.width(), self.toolBar.height())



class AnalysisInterface(ScrollArea):
    """ Analysis Interface """

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.view = TabInterface(self)
        
        #self.toolBar = ToolBar(title, subtitle, self)
        #self.vBoxLayout = QVBoxLayout(self.view)

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setViewportMargins(0, 0, 0, 0)
        self.setWidget(self.view)
        self.setWidgetResizable(True)

        #self.vBoxLayout.setSpacing(30)
        #self.vBoxLayout.setAlignment(Qt.AlignTop)
        #self.vBoxLayout.setContentsMargins(36, 20, 36, 36)

        self.view.setObjectName('view')
        StyleSheet.GALLERY_INTERFACE.apply(self)
        self.setObjectName('AnalysisInterface')

        '''
        card = self.addExampleCard(
            title=self.tr('A tab bar'),
            widget=TabInterface(self),
            sourcePath='https://github.com/zhiyiYo/PyQt-Fluent-Widgets/blob/master/examples/navigation/tab_view/demo.py',
            stretch=1
        )
        self.vBoxLayout.addWidget(card, 0, Qt.AlignTop)        
        card.topLayout.setContentsMargins(12, 0, 0, 0)
        '''
        #self.vBoxLayout.addWidget(, 0, Qt.AlignTop)
        
        #self.addWidget(TabInterface(self))
        

        #self.iconView = IconCardView(self)
        #self.vBoxLayout.addWidget(self.iconView)

class RadialVelocityCard(QWidget):
    """ Radial Velocity card """

    def __init__(self, title, widget: QWidget, sourcePath, stretch=0, parent=None):
        super().__init__(parent=parent)
        self.widget = widget
        self.stretch = stretch

        self.titleLabel = StrongBodyLabel(title, self)
        self.card = QFrame(self)

        self.sourceWidget = QFrame(self.card)
        self.sourcePath = sourcePath
        self.sourcePathLabel = BodyLabel(
            self.tr('Source code'), self.sourceWidget)
        self.linkIcon = IconWidget(FluentIcon.LINK, self.sourceWidget)

        self.vBoxLayout = QVBoxLayout(self)
        self.cardLayout = QVBoxLayout(self.card)
        self.topLayout = QHBoxLayout()
        self.bottomLayout = QHBoxLayout(self.sourceWidget)

        self.__initWidget()

    def __initWidget(self):
        self.linkIcon.setFixedSize(16, 16)
        self.__initLayout()

        self.sourceWidget.setCursor(Qt.PointingHandCursor)
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

        self.bottomLayout.addWidget(self.sourcePathLabel, 0, Qt.AlignLeft)
        self.bottomLayout.addStretch(1)
        self.bottomLayout.addWidget(self.linkIcon, 0, Qt.AlignRight)
        self.bottomLayout.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)

    def eventFilter(self, obj, e):
        if obj is self.sourceWidget:
            if e.type() == QEvent.MouseButtonRelease:
                QDesktopServices.openUrl(QUrl(self.sourcePath))

        return super().eventFilter(obj, e)



'''
class ExampleCard(QWidget):
    """ Example card """

    def __init__(self, title, widget: QWidget, sourcePath, stretch=0, parent=None):
        super().__init__(parent=parent)
        self.widget = widget
        self.stretch = stretch

        self.titleLabel = StrongBodyLabel(title, self)
        self.card = QFrame(self)

        self.sourceWidget = QFrame(self.card)
        self.sourcePath = sourcePath
        self.sourcePathLabel = BodyLabel(
            self.tr('Source code'), self.sourceWidget)
        self.linkIcon = IconWidget(FluentIcon.LINK, self.sourceWidget)

        self.vBoxLayout = QVBoxLayout(self)
        self.cardLayout = QVBoxLayout(self.card)
        self.topLayout = QHBoxLayout()
        self.bottomLayout = QHBoxLayout(self.sourceWidget)

        self.__initWidget()

    def __initWidget(self):
        self.linkIcon.setFixedSize(16, 16)
        self.__initLayout()

        self.sourceWidget.setCursor(Qt.PointingHandCursor)
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

        self.bottomLayout.addWidget(self.sourcePathLabel, 0, Qt.AlignLeft)
        self.bottomLayout.addStretch(1)
        self.bottomLayout.addWidget(self.linkIcon, 0, Qt.AlignRight)
        self.bottomLayout.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)

    def eventFilter(self, obj, e):
        if obj is self.sourceWidget:
            if e.type() == QEvent.MouseButtonRelease:
                QDesktopServices.openUrl(QUrl(self.sourcePath))

        return super().eventFilter(obj, e)

'''