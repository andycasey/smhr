
import threading
from PyQt5.QtCore import Qt, QSize, QTimer, pyqtSignal, QObject, QEvent
from PyQt5.QtGui import QResizeEvent
from PyQt5.QtWidgets import QWidget, QStackedWidget, QVBoxLayout, QLabel, QHBoxLayout, QFrame, QSizePolicy

from PyQt5.QtWidgets import QApplication, QFrame, QVBoxLayout, QLabel, QWidget, QHBoxLayout
from qfluentwidgets import (FluentIcon, IconWidget, FlowLayout, isDarkTheme,
                            Theme, applyThemeColor, SmoothScrollArea, SearchLineEdit, StrongBodyLabel,
                            BodyLabel, CommandBar)
from qfluentwidgets import RoundMenu, Action
from time import time

# Remember to import matplotlib after Qt.
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvas

# TODO: Import qconfig from qfluentwidgets and when the theme changes, update the matplotlib style file
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
        
        
class ExcitationIonizationBalanceWidget(QWidget):    
    
    def __init__(
        self, 
        x=None,
        y=None,
        xlabel=None,
        ylabel=None,
        figsize=(8, 6),
        parent=None, 
        size_policy=(QSizePolicy.Expanding, QSizePolicy.Fixed),
        resize_interval=50
    ):
        super().__init__(parent)
        self.parent = parent
        self.resize_interval = resize_interval
        self.figure = Figure(figsize=figsize)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setParent(self)
        self.canvas.setFocusPolicy(Qt.WheelFocus)
        self.canvas.setFocus()
        self.canvas.setSizePolicy(*size_policy)
        
        self.canvas.mpl_connect("figure_enter_event", self._focus)
                    
        self.axes = self.canvas.figure.subplots(3, 1)
        self.figure.tight_layout()
        self.figure.canvas.draw()
        
        self.layout = QVBoxLayout(self)
        self.layout.addWidget(self.canvas)

        #self.layout.addStretch(1)    
        self.installEventFilter(self)
        return None
    
    def _focus(self, event):
        """ Set the focus of the canvas. """
        self.canvas.setFocus()    
        
    
    def resizeEvent(self, e):
        # Matplotlib wants to redraw the canvas while we are resizing, makes it all yucky
        try:
            self.resizeTimer
        except:
            
            self.resizeTimer = QTimer(self)
            self.resizeTimer.setSingleShot(True)
            self.resizeTimer.timeout.connect(self.afterResizeEvent)
        finally:
            self.resizeTimer.start(self.resize_interval)
                
        return None        

    def afterResizeEvent(self):        
        self.figure.tight_layout()
        self.figure.canvas.draw()
        try:
            self.resizeTimer.stop()
            del self.resizeTimer
        except:
            None
        return None
    
    def eventFilter(self, widget, event):
        try:
            print(f"plot widget {widget} {event.type()} {event.key()} {event.text()}")
        except:
            None
            
        if event.type() == 51:
            if event.key() == Qt.Key_Left:
                try:
                    self.page_left.trigger()
                except:
                    return False
                else:
                    return True
            elif event.key() == Qt.Key_Right:
                try:
                    self.page_right.trigger()
                except:
                    return False
                else:
                    return True
                
                                

        '''
        if event.type() == QEvent.KeyPress:
            text = event.text()
            if event.modifiers():
                text = event.keyCombination().key().name.decode(encoding="utf-8")
            print(f"{event} {event.type}: {text}")
        '''            
        return False
        

        
class SinglePlotWidget(QWidget):    
    
    def __init__(
        self, 
        x=None,
        y=None,
        xlabel=None,
        ylabel=None,
        figsize=(8, 2),
        parent=None, 
        size_policy=(QSizePolicy.Expanding, QSizePolicy.Fixed),
        toolbar=False,
        toolbar_left_right=True,
        resize_interval=50
    ):
        super().__init__(parent)
        self.parent = parent
        self.resize_interval = resize_interval
        self.figure = Figure(figsize=figsize)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setParent(self)
        self.canvas.setFocusPolicy(Qt.WheelFocus)
        self.canvas.setFocus()
        self.canvas.setSizePolicy(*size_policy)
        
        self.canvas.mpl_connect("figure_enter_event", self._focus)
                    
        self.ax = self.canvas.figure.subplots()
        if x is not None and y is not None:
            self.ax.plot(x, y)
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.figure.tight_layout()
        self.figure.canvas.draw()
        
        self.layout = QVBoxLayout(self)
        self.layout.addWidget(self.canvas)
        
                
        if toolbar:
            toolbar = CommandBar(parent=self)
            toolbar.setFocusPolicy(Qt.NoFocus)
            
            action_home = Action(FluentIcon.HOME, "Home", self)
            action_zoom = Action(FluentIcon.ZOOM_IN, "Zoom", self)
            action_zoom.setCheckable(True)
            action_pan = Action(FluentIcon.MOVE, "Pan", self)
            action_pan.setCheckable(True)
            action_home.triggered.connect(lambda: print("home"))
            action_zoom.triggered.connect(lambda: print("zoom"))
            action_pan.triggered.connect(lambda: print("pan"))

            toolbar.addAction(action_home)
            toolbar.addAction(action_zoom)
            toolbar.addAction(action_pan)
        
            if toolbar_left_right:
                self.page_left = Action(FluentIcon.PAGE_LEFT, "Left", self)
                self.page_right = Action(FluentIcon.PAGE_RIGHT, "Right", self)
                toolbar.addSeparator()
                toolbar.addAction(self.page_left)
                toolbar.addAction(self.page_right)
            self.layout.addWidget(toolbar)
        
        else:
            if toolbar_left_right:
                # Only left/right, set to align in middle
                toolbar = CommandBar(parent=self)
                toolbar.setFocusPolicy(Qt.NoFocus)
                self.page_left = Action(FluentIcon.PAGE_LEFT, "Left", self)
                self.page_right = Action(FluentIcon.PAGE_RIGHT, "Right", self)
                toolbar.addAction(self.page_left)
                toolbar.addAction(self.page_right)
                self.layout.addWidget(toolbar, alignment=Qt.AlignCenter)           
        
        #self.layout.addStretch(1)    
        self.installEventFilter(self)
        return None
    
    def _focus(self, event):
        """ Set the focus of the canvas. """
        self.canvas.setFocus()    
        
    
    def resizeEvent(self, e):
        # Matplotlib wants to redraw the canvas while we are resizing, makes it all yucky
        try:
            self.resizeTimer
        except:
            
            self.resizeTimer = QTimer(self)
            self.resizeTimer.setSingleShot(True)
            self.resizeTimer.timeout.connect(self.afterResizeEvent)
        finally:
            self.resizeTimer.start(self.resize_interval)
                
        return None        

    def afterResizeEvent(self):        
        self.figure.tight_layout()
        self.figure.canvas.draw()
        try:
            self.resizeTimer.stop()
            del self.resizeTimer
        except:
            None
        return None
    
    def eventFilter(self, widget, event):
        try:
            print(f"plot widget {widget} {event.type()} {event.key()} {event.text()}")
        except:
            None
            
                    
        if event.type() == 51:
            if event.key() == Qt.Key_Left:
                try:
                    self.page_left.trigger()
                except:
                    return False
                else:
                    return True
            elif event.key() == Qt.Key_Right:
                try:
                    self.page_right.trigger()
                except:
                    return False
                else:
                    return True
                
                                

        '''
        if event.type() == QEvent.KeyPress:
            text = event.text()
            if event.modifiers():
                text = event.keyCombination().key().name.decode(encoding="utf-8")
            print(f"{event} {event.type}: {text}")
        '''            
        return False
        
    
if __name__ == '__main__':
    # enable dpi scale
    from PyQt5.QtCore import QModelIndex, Qt, QRect, QSize
    import sys
    QApplication.setHighDpiScaleFactorRoundingPolicy(Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps)

    app = QApplication(sys.argv)
    w = SinglePlotWidget()
    w.resize(600, 600)
    w.show()
    app.exec_()    