# coding: utf-8
from PyQt5.QtCore import QUrl, QSize, QObject, QEvent
from PyQt5.QtGui import QIcon, QDesktopServices, QKeyEvent
from PyQt5.QtWidgets import QApplication, QFileDialog

from qfluentwidgets import (Action, NavigationAvatarWidget, NavigationItemPosition, MessageBox, FluentWindow,
                            SplashScreen, FolderListDialog, NavigationWidget)
from qfluentwidgets import FluentIcon as FIF

from .home_interface import HomeInterface
from .setting_interface import SettingInterface
from .analysis_interface import AnalysisInterface

from ..common.config import NAME, SUPPORT_URL, cfg
from ..common.icon import Icon
from ..common.signal_bus import signalBus
from ..common.translator import Translator
from ..common import resource


class KeyPressFilter(QObject):

    def eventFilter(self, widget, event):
        try:
            print(f"__debug: {widget} {event.type()} {event.key()} {event.text()}")
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
    



class MainWindow(FluentWindow):

    def __init__(self):
        super().__init__()
        self.initWindow()

        # create sub interface
        self.homeInterface = HomeInterface(self)
        self.analysisInterface = AnalysisInterface(self)
        self.settingInterface = SettingInterface(self)

        # enable acrylic effect
        self.navigationInterface.setAcrylicEnabled(True)

        self.connectSignalToSlot()

        # add items to navigation interface
        self.initNavigation()
        self.splashScreen.finish()
                
        self.installEventFilter(KeyPressFilter(parent=self))        
        

    def connectSignalToSlot(self):
        signalBus.micaEnableChanged.connect(self.setMicaEffectEnabled)
        signalBus.switchToSampleCard.connect(self.switchToSample)
        #signalBus.supportSignal.connect(self.onSupport)

    def initNavigation(self):
        # add navigation items
        t = Translator()
        self.addSubInterface(self.homeInterface, FIF.HOME, self.tr('Home'))
        self.addSubInterface(self.analysisInterface, FIF.APPLICATION, self.tr('Analysis'))
        self.navigationInterface.addSeparator()
        
        # add an action for new analysis
        self.navigationInterface.addItem(
            "newAnalysis",
            FIF.ADD, 
            "New analysis",
            onClick=self.analysisInterface.view.addTab,
            tooltip="Start new analysis",
            selectable=False
        )


        # add custom widget to bottom
        self.navigationInterface.addWidget(
            routeKey='avatar',
            widget=NavigationAvatarWidget(cfg.get(cfg.userName), ':/gallery/images/shoko.png'),
            onClick=self.onAvatarClicked,
            position=NavigationItemPosition.BOTTOM
        )
        self.addSubInterface(
            self.settingInterface, FIF.SETTING, self.tr('Settings'), NavigationItemPosition.BOTTOM)

    def initWindow(self):
        self.resize(960, 780)
        self.setMinimumWidth(760)
        self.setWindowIcon(QIcon(':/gallery/images/logo.png'))
        self.setWindowTitle(NAME)

        self.setMicaEffectEnabled(cfg.get(cfg.micaEnabled))

        # create splash screen
        self.splashScreen = SplashScreen(self.windowIcon(), self)
        self.splashScreen.setIconSize(QSize(106, 106))
        self.splashScreen.raise_()

        desktop = QApplication.desktop().availableGeometry()
        w, h = desktop.width(), desktop.height()
        self.move(w//2 - self.width()//2, h//2 - self.height()//2)
        self.show()
        QApplication.processEvents()


    def onAvatarClicked(self):
        pass


    def resizeEvent(self, e):
        super().resizeEvent(e)
        if hasattr(self, 'splashScreen'):
            self.splashScreen.resize(self.size())

    def switchToSample(self, routeKey, index):
        """ switch to sample """
        interfaces = self.findChildren(HomeInterface)
        for w in interfaces:
            if w.objectName() == routeKey:
                self.stackedWidget.setCurrentWidget(w, False)
                w.scrollToCard(index)
