


from PyQt5.QtWidgets import QWidget, QStackedWidget, QVBoxLayout, QLabel, QHBoxLayout, QFrame, QSizePolicy

from PyQt5.QtWidgets import QApplication, QFrame, QVBoxLayout, QLabel, QWidget, QHBoxLayout
from qfluentwidgets import (FluentIcon, IconWidget, FlowLayout, isDarkTheme,
                            Theme, applyThemeColor, SmoothScrollArea, SearchLineEdit, StrongBodyLabel,
                            BodyLabel, CommandBar)
from qfluentwidgets import RoundMenu, Action

# Remember to import matplotlib after Qt.
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvas

# TODO: Import qconfig from qfluentwidgets and when the theme changes, update the matplotlib style file

class SinglePlotWidget(QWidget):
    
    def __init__(
        self, 
        parent=None, 
        figsize=(8, 2),
        size_policy=(QSizePolicy.Expanding, QSizePolicy.Fixed),
        toolbar=False,
        toolbar_left_right=True,
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

            toolbar.addAction(action_home)
            toolbar.addAction(action_zoom)
            toolbar.addAction(action_pan)
        
            if toolbar_left_right:
                action_left = Action(FluentIcon.PAGE_LEFT, "Left", self)
                action_right = Action(FluentIcon.PAGE_RIGHT, "Right", self)
                action_left.triggered.connect(lambda: print("left triggered"))
                action_left.hovered.connect(lambda: print("left hovered"))
                toolbar.addSeparator()
                toolbar.addAction(action_left)
                toolbar.addAction(action_right)
            self.layout.addWidget(toolbar)
        
        else:
            if toolbar_left_right:
                # Only left/right, set to align in middle
                toolbar = CommandBar(parent=self)
                action_left = Action(FluentIcon.PAGE_LEFT, "Left", self)
                action_right = Action(FluentIcon.PAGE_RIGHT, "Right", self)
                action_left.triggered.connect(lambda: print("left triggered"))
                action_left.hovered.connect(lambda: print("left hovered"))
                toolbar.addAction(action_left)
                toolbar.addAction(action_right)
                self.layout.addWidget(toolbar, alignment=Qt.AlignCenter)           
        
        self.layout.addStretch(1)    
        return None
    

    
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