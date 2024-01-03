# coding:utf-8
import time
from typing import List
import numpy as np
import matplotlib.pyplot as plt
from PyQt5.QtCore import QEvent, Qt, QUrl, QObject
from PyQt5.QtGui import QColor, QDesktopServices, QIcon, QPainter, QPen
from PyQt5.QtWidgets import (QApplication, QCompleter, QFileDialog, QFrame,
                             QHBoxLayout, QLabel, QSizePolicy, QStackedWidget,
                             QVBoxLayout, QWidget)
from qfluentwidgets import (Action, BodyLabel, BreadcrumbBar, CaptionLabel,
                            CardWidget, CheckBox, ComboBox, CommandBar,
                            EditableComboBox, FlowLayout)
from qfluentwidgets import (IconWidget, LineEdit, Pivot, PrimaryPushButton,
                            PushButton, RoundMenu, ScrollArea, SearchLineEdit,
                            SmoothScrollArea, StrongBodyLabel, TabBar,
                            TabCloseButtonDisplayMode, TextWrap, Theme,
                            TitleLabel, ToolButton, ToolTipFilter,
                            applyThemeColor, isDarkTheme, qrouter, toggleTheme)

from ..common.config import EXAMPLE_URL, FEEDBACK_URL, HELP_URL, cfg
from ..common.icon import Icon
from ..common.style_sheet import StyleSheet
from ..components.tool_bar import ToolBar

from ..components.plot import SinglePlotWidget
from ..view.gallery_interface import SeparatorWidget



class ExampleCard2(QWidget):
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


class KeyPressFilter(QObject):

    def eventFilter(self, widget, event):
        try:
            print(f"ok: {event.type()} {event.key()} {event.text()}")
        except:
            None
        '''
        if event.type() == QEvent.KeyPress:
            text = event.text()
            if event.modifiers():
                text = event.keyCombination().key().name.decode(encoding="utf-8")
            print(f"{event} {event.type}: {text}")
        '''            
        return False
    

class AnalysisWidget(QWidget):
    
    def __init__(self, session, title=None, parent=None, stretch=0):
        super().__init__(parent=parent)
        self.session = session
        self.stretch = stretch
        self.widget = QWidget(self)

        self.titleLabel = StrongBodyLabel(title, self)
        self.card = QFrame(self)

        self.sourceWidget = QFrame(self.card)
        
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

        self.topLayout.addWidget(self.widget)
        if self.stretch == 0:
            self.topLayout.addStretch(1)


        '''
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

        '''
        self.bottomLayout.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        self.widget.setParent(self.card)
        self.widget.show()
        
        
    def eventFilter(self, widget, event):
        try:
            print(f"ok: {event.type()} {event.key()} {event.text()}")
        except:
            None
        '''
        if event.type() == QEvent.KeyPress:
            text = event.text()
            if event.modifiers():
                text = event.keyCombination().key().name.decode(encoding="utf-8")
            print(f"{event} {event.type}: {text}")
        '''            
        return False        
    


class ContinuumRectificationWidget(AnalysisWidget):
    
    def __init__(self, session, parent=None, stretch=0):
        
        # Add buttons to bottom layout
        super().__init__(parent=parent, session=session, title="Continuum rectification")
        layout = QVBoxLayout(self.widget)
        layout.setContentsMargins(0, 0, 0, 0)  

        import os
        
        def expand_path(path):
            return os.path.abspath(os.path.expanduser(path))


        from specutils import Spectrum1D

        wavelength, flux, ivar, meta = Spectrum1D.read_fits_multispec(expand_path("~/Downloads/j174239-133332red_multi.fits"), flux_ext=6, ivar_ext=3)
        
        from continuum import load_basis_vectors, BaseContinuumModel
        
        '''
        model_wavelength = 10 * (10**(2.57671464 + np.arange(167283) * 2.25855074e-06))
        basis_vectors = load_basis_vectors("H.pkl.gz", wavelength.size)

        model = BaseContinuumModel(model_wavelength, basis_vectors, deg=3)

        index = 30
        #(v_rel, continuum, template_wavelength, template_flux, p_opt) = model.fit_relative_velocity(wavelength[index], flux[index], ivar[index])

        #print(v_rel)
        
        print(continuum)
        index = 28
        '''
        
        figure = SinglePlotWidget(
            xlabel=r"Wavelength [$\mathrm{\AA}$]",
            ylabel="Rectified flux",
            figsize=(100, 4), 
            toolbar=True, 
            toolbar_left_right=True
        )
        current_index = 0
        
        line, = figure.ax.plot([], [])
        
        def update_canvas(index):
            nonlocal current_index
            if index < 0 or index > (len(wavelength) - 1):
                return
            line.set_data(wavelength[index], flux[index] / np.nanmedian(flux[index]))
            figure.ax.set_xlim(wavelength[index][[0, -1]])
            figure.ax.set_ylim(0, 1.2)
            figure.canvas.draw()
            current_index = index
            figure.page_left.setEnabled(current_index > 0)
            figure.page_right.setEnabled(current_index < (len(wavelength) - 1))
            self.canvas.setFocus()
            
        
        figure.page_left.triggered.connect(lambda: update_canvas(current_index - 1))
        figure.page_right.triggered.connect(lambda: update_canvas(current_index + 1))
        
        update_canvas(current_index)

        layout.addWidget(figure)                
        
        self.button_do_not_rectify = PushButton("No continuum normalization")
        self.button_rectify = PrimaryPushButton("Rectify")
        self.bottomLayout.addWidget(self.button_do_not_rectify)
        self.bottomLayout.addStretch(1)
        self.bottomLayout.addWidget(self.button_rectify)

        self.widget.setParent(self.card)
        self.widget.show()
        
        return None
    
    
    
class RadialVelocityWidget(AnalysisWidget):
    
    
    def __init__(self, session, callback=None, parent=None, stretch=0):
        
        # Add buttons to bottom layout
        super().__init__(parent=parent, session=session, title="Radial velocity correction")
        self.callback = callback
        '''
        layout = QVBoxLayout(self.widget)
        layout.addWidget(QLabel("Test"))
        
        
        self.button_do_not_rectify = PushButton("No continuum normalization")
        self.button_rectify = PrimaryPushButton("Rectify")
        self.bottomLayout.addWidget(self.button_do_not_rectify)
        self.bottomLayout.addStretch(1)
        self.bottomLayout.addWidget(self.button_rectify)
        '''
            
        from specutils import Spectrum1D
        
        import os
        
        def expand_path(path):
            return os.path.abspath(os.path.expanduser(path))

        wavelength, flux, ivar, meta = Spectrum1D.read_fits_multispec(expand_path("~/Downloads/hd122563red_multi.fits"), flux_ext=6, ivar_ext=3)
        
        layout = QVBoxLayout(self.widget)
        layout.setContentsMargins(0, 0, 0, 0)  
        index = np.random.choice(len(wavelength))
        
        #index = 28
        plot = SinglePlotWidget(
            x=wavelength[index],
            y=flux[index] / np.nanmedian(flux[index]),            
            xlabel=r"Wavelength [$\mathrm{\AA}$]",
            ylabel="Rectified flux",
            figsize=(100, 2), toolbar=False, toolbar_left_right=False)
        plot.ax.set_xlim(wavelength[index][[0, -1]])
        
        layout.addWidget(plot)

        
        
        rv_opt_layout = QHBoxLayout()
        rv_opt_layout.setSpacing(30)
        
        lhs_layout = QVBoxLayout()
        middle_layout = QVBoxLayout()
        rhs_layout = QVBoxLayout()
        
        
        wavelength_layout = QHBoxLayout()
        wavelength_layout.addWidget(QLabel("Wavelength range:"))
        wavelength_combo = ComboBox()
        wavelength_combo.addItems([
            "8000 - 9000 Å"
        ])
        wavelength_layout.addWidget(wavelength_combo)


        template_layout = QHBoxLayout()
        template_combo = ComboBox()
        template_combo.setPlaceholderText("Auto")
        template_combo.addItems([
            "Auto",
            "Select from file",
            "Synthesise spectrum"
        ])
        for combo in (wavelength_combo, template_combo):
            combo.setFixedWidth(180)
        template_layout.addWidget(QLabel("Template:"))
        template_layout.addWidget(template_combo)

        lhs_layout.addLayout(wavelength_layout)
        lhs_layout.addLayout(template_layout)

    
        # Middle layout
        continuum_layout = QHBoxLayout()
        continuum_layout.setAlignment(Qt.AlignTop)
        continuum_layout.addWidget(QLabel("Continuum:"))
        continuum_combo = ComboBox()
        continuum_combo.setPlaceholderText("Sinusoids")
        continuum_combo.addItems([
            "Sinusoids",
            "Polynomial",
            "Spline",
        ])
        continuum_layout.addWidget(continuum_combo)        
        middle_layout.addLayout(continuum_layout)
        

        rv_opt_layout.addLayout(lhs_layout)
        rv_opt_layout.addLayout(middle_layout)
        rv_opt_layout.addStretch(1)
        rv_opt_layout.addLayout(rhs_layout)    
        
        layout.addLayout(rv_opt_layout)
                
        self.button_no_shift = PushButton("Do not shift")
        self.button_measure = PushButton("Measure")
        self.button_measure.clicked.connect(self.on_shift_button_clicked)
        self.button_measure_and_shift = PrimaryPushButton("Measure and shift")
        self.button_measure_and_shift.clicked.connect(self.on_measure_and_shift_button_clicked)
        self.bottomLayout.addWidget(self.button_no_shift)
        self.bottomLayout.addStretch(1)
        self.bottomLayout.addWidget(QLabel("Relative velocity:"))
        
        self.v_rel = LineEdit(self)
        self.v_rel.setText("0.0")
        self.v_rel.setFixedWidth(80)
        self.v_rel.setAlignment(Qt.AlignHCenter)
        self.v_rel.textEdited.connect(self.on_v_rel_changed)
        
        self.bottomLayout.addWidget(self.v_rel)
        self.bottomLayout.addWidget(QLabel("km/s"))
        self.bottomLayout.addWidget(SeparatorWidget())
        self.bottomLayout.addWidget(self.button_measure)
        self.bottomLayout.addWidget(self.button_measure_and_shift)
        
        self.button_measure.setVisible(False)
        #measure_layout = QHBoxLayout()
        #measure_layout.addStretch(1)
        #measure_layout.addWidget(PushButton("Measure relative velocity"), 0, Qt.AlignRight)
        #layout.addLayout(measure_layout)
                        
        return None    
    
    def on_no_shift_button_pushed(self):
        self.v_rel.setText("0.0")
        if self.callback is not None:
            self.callback()
            
    def on_measure_and_shift_button_clicked(self):
        print("on measure clicked")
        if self.callback is not None:
            self.callback()
            
    def on_shift_button_clicked(self):
        print("on shift clicked")
        if self.callback is not None:
            self.callback()
    
    def on_v_rel_changed(self):
        self.button_measure.setVisible(True)
        self.button_measure_and_shift.setText("Shift")
        print("edited")
        


class SessionInterface(ScrollArea):
    """ Dialog interface """

    def __init__(self, text, parent=None):
        super().__init__(parent=parent)
        self.view = QWidget(self)
        self.toolBar = ToolBar(text, "15:43:03 -10:56:01", self)
        self.vBoxLayout = QVBoxLayout(self.view)

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setViewportMargins(0, self.toolBar.height(), 0, 0)
        self.setWidget(self.view)
        self.setWidgetResizable(True)

        self.vBoxLayout.setSpacing(30)
        self.vBoxLayout.setAlignment(Qt.AlignTop)
        self.vBoxLayout.setContentsMargins(36, 20, 36, 36)

        self.view.setObjectName('view')
        StyleSheet.ANALYSIS_INTERFACE.apply(self)
        
        self.setObjectName('fooInterface')

        button1 = PushButton(self.tr('Show dialog'))
        button1.clicked.connect(self.showDialog)


        # Radial velocity stuff
        self.continuum = ContinuumRectificationWidget(None, self)
        self.continuum.setVisible(False)
        
        def rv_callback():
            self.continuum.setVisible(True)
        
        card = RadialVelocityWidget(None, callback=rv_callback, parent=self)
        self.vBoxLayout.addWidget(card, 0, Qt.AlignTop)
        
        # Continuum normalization

        self.vBoxLayout.addWidget(self.continuum, 0, Qt.AlignTop)
        
                    
        
    def _update_canvas(self):
        t = np.linspace(0, 10, 101)
        # Shift the sinusoid as a function of time.
        self._line.set_data(t, np.sin(t + time.time()))
        self._line.figure.tight_layout()
        self._line.figure.canvas.draw()

                                                        

    def addExampleCard2(self, title, widget, sourcePath: str, stretch=0):
        card = ExampleCard2(title, widget, sourcePath, stretch, self.view)
        self.vBoxLayout.addWidget(card, 0, Qt.AlignTop)
        return card
        
    def showDialog(self):
        print("NAH")


class SessionTabsInterface(QWidget):
    """ Session tabs interface """

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
        
        self.connectSignalToSlot()

        self.tabBar.setCloseButtonDisplayMode(TabCloseButtonDisplayMode.ON_HOVER)

        self.myInterface = SessionInterface("HD 140283", self)        
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
            widget = SessionInterface(text, self)

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
        self.view = SessionTabsInterface(self)        
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setViewportMargins(0, 0, 0, 0)
        self.setWidget(self.view)
        self.setWidgetResizable(True)
        self.view.setObjectName('view')
        #StyleSheet.GALLERY_INTERFACE.apply(self)
        self.setObjectName('AnalysisInterface')

