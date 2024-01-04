# coding:utf-8
import time
from typing import List
import numpy as np
import matplotlib.pyplot as plt
from PyQt5 import QtCore
from PyQt5.QtCore import QEvent, Qt, QUrl, QObject, pyqtSignal, QItemSelectionModel
from PyQt5.QtGui import QColor, QDesktopServices, QIcon, QPainter, QPen, QDoubleValidator, QIntValidator
from PyQt5.QtWidgets import (QApplication, QCompleter, QFileDialog, QFrame,
                             QGridLayout,
                             QHBoxLayout, QLabel, QSizePolicy, QStackedWidget, 
                             QVBoxLayout, QWidget, QTableWidgetItem, QHeaderView, QTableView)
from qfluentwidgets import (Action, BodyLabel, BreadcrumbBar, CaptionLabel,
                            CardWidget, CheckBox, ComboBox, CommandBar,
                            EditableComboBox, FlowLayout, Flyout)
from qfluentwidgets import (IconWidget, LineEdit, Pivot, PrimaryPushButton, InfoBarIcon,
                            PushButton, RoundMenu, ScrollArea, SearchLineEdit,
                            SmoothScrollArea, StrongBodyLabel, TabBar,
                            TabCloseButtonDisplayMode, TextWrap, Theme,
                            TitleLabel, ToolButton, ToolTipFilter, TableWidget, SmoothScrollDelegate,
                            applyThemeColor, isDarkTheme, qrouter, toggleTheme)
from qfluentwidgets import TableWidget, isDarkTheme, setTheme, Theme, TableView, TableItemDelegate, setCustomStyleSheet

from qfluentwidgets import TransparentDropDownPushButton, PrimaryDropDownPushButton, PrimarySplitPushButton

from ..common.config import EXAMPLE_URL, FEEDBACK_URL, HELP_URL, cfg
from ..common.icon import Icon
from ..common.style_sheet import StyleSheet
from ..components.tool_bar import ToolBar

from ..components.plot import SinglePlotWidget, ExcitationIonizationBalanceWidget
from ..view.gallery_interface import SeparatorWidget

from qfluentwidgets import RoundMenu, setTheme, Theme, Action, MenuAnimationType, MenuItemDelegate, CheckableMenu, MenuIndicatorType
from qfluentwidgets import FluentIcon as FIF


from session import Session
from astropy.io.fits import getval

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

class TopAndBottomAnalysisWidget(QWidget):
    
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
    
    def __init__(self, session, callback=None, parent=None, stretch=0):
        # Add buttons to bottom layout
        super().__init__(parent=parent, session=session, title="Continuum Rectification")
        self.callback = callback
        
        layout = QVBoxLayout(self.widget)
        layout.setContentsMargins(0, 0, 0, 0)  

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
        
        S = len(session.spectra)
        
        def update_canvas(index):
            nonlocal current_index
            if index < 0 or index > (S - 1):
                return
            wavelength, flux, ivar, meta = session.spectra[index]
            line.set_data(wavelength, flux / np.nanmedian(flux))
            figure.ax.set_xlim(wavelength[[0, -1]])
            figure.ax.set_ylim(0, 1.2)
            figure.canvas.draw()
            current_index = index
            figure.page_left.setEnabled(current_index > 0)
            figure.page_right.setEnabled(current_index < (S - 1))
            figure.canvas.setFocus()
            
        
        figure.page_left.triggered.connect(lambda: update_canvas(current_index - 1))
        figure.page_right.triggered.connect(lambda: update_canvas(current_index + 1))
        
        update_canvas(current_index)

        layout.addWidget(figure)                
        
        self.button_do_not_rectify = PushButton("No continuum normalization")
        self.button_rectify = PrimaryPushButton("Continuum normalize all spectra")
        
        self.button_do_not_rectify.clicked.connect(self.button_do_not_rectify_clicked)
        self.button_rectify.clicked.connect(self.button_rectify_clicked)
        self.bottomLayout.addWidget(self.button_do_not_rectify)
        self.bottomLayout.addStretch(1)
        self.bottomLayout.addWidget(self.button_rectify)

        self.widget.setParent(self.card)
        self.widget.show()
        
        return None
    
    def button_do_not_rectify_clicked(self):
        if self.callback is not None:
            self.callback()
        
    def button_rectify_clicked(self):
        if self.callback is not None:
            self.callback()
            

from PyQt5.QtCore import QModelIndex, Qt, QVariant

from PyQt5.QtGui import QPalette
    
from PyQt5.QtWidgets import QApplication, QStyleOptionViewItem, QTableWidget, QTableWidgetItem, QWidget, QHBoxLayout
    
class CustomTableItemDelegate(TableItemDelegate):
    """ Custom table item delegate """

    def initStyleOption(self, option: QStyleOptionViewItem, index: QModelIndex):
        super().initStyleOption(option, index)
        if index.column() != 1:
            return

        if isDarkTheme():
            option.palette.setColor(QPalette.Text, Qt.white)
            option.palette.setColor(QPalette.HighlightedText, Qt.white)
        else:
            option.palette.setColor(QPalette.Text, Qt.red)
            option.palette.setColor(QPalette.HighlightedText, Qt.red)    

class TransitionsTableModel(QtCore.QAbstractTableModel):

    header = ["", "λ", "Species", "χ", "log(gf)", "Γ", "Deg", "Tol", "Profile", "EW", "A(X)"]
    attrs = (
        "is_acceptable",
        "_repr_wavelength", 
        "_repr_element", 
        "chi",
        "loggf",
        "C6",
        "poly",
        "wl_tol",
        "option",
        "equivalent_width",
        "A(X)",
    )

    def __init__(self, parent, data, *args):
        super().__init__(parent, *args)
        self._data = data

    def rowCount(self, parent):
        return len(self._data)

    def columnCount(self, parent):
        return len(self.header)

    def data(self, index, role):
        if not index.isValid():
            return None

        value = self._data[index.row()][index.column()]

        if index.column() == 0:
            if role == QtCore.Qt.CheckStateRole:
                return value #QtCore.Qt.Checked if value else QtCore.Qt.Unchecked
            else:
                return None

        elif role != QtCore.Qt.DisplayRole:
            return None

        return QVariant(value)
    

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal \
        and role == QtCore.Qt.DisplayRole:
            return self.header[col]
        return None


    def setData(self, index, value, role=QtCore.Qt.DisplayRole):
        print(f"setData {index} {value}")
        attr = self.attrs[index.column()]
        if attr != "is_acceptable":
            return False

        old_value = self.data(index, role)
        
        # Enable tri-state for first column, where the order is :
        # false, true, upper limit
        print(f"setting {index} {role} to {value} from {old_value}")
        if old_value == 2:
            value = 1
        elif old_value == 1:
            value = 0
                        
        self._data[index.row()][index.column()] = value
        self.dataChanged.emit(index, index)
        return value

    def sort(self, column, order):
        #self.emit(QtCore.SIGNAL("layoutAboutToBeChanged()"))
        print(column, order)
        self._data = sorted(self._data, key=lambda sm: sm[column])
        
        if order == QtCore.Qt.DescendingOrder:
            self._data.reverse()

        self.dataChanged.emit(self.createIndex(0, 0),
            self.createIndex(self.rowCount(0), self.columnCount(0)))
        #self.emit(QtCore.SIGNAL("layoutChanged()"))

    def flags(self, index):
        #if not index.isValid(): return
        flags = QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable
        if index.column() == 0:
            flags |= QtCore.Qt.ItemIsUserCheckable 
        elif index.column() == 8:
            flags |= QtCore.Qt.ItemIsEditable
        return flags

# import requried things
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QItemDelegate, QComboBox, QListWidget, QListWidgetItem, QStyleOptionViewItem, QStyle, QHeaderView

from qfluentwidgets import FluentIcon as FIF
class ComboDelegate(QItemDelegate):
    editorItems=['Gaussian', 'Voight','Lorentzian']
    #height = 25
    #width = 200
    def createEditor(self, parent, option, index):
        print(f"creating editor with parent {parent} {option} {index}")
        editor = TransparentDropDownPushButton(parent)

        self.menu = RoundMenu()
        self.menu.addAction(Action('Gaussian'))
        self.menu.addAction(Action('Voight'))
        self.menu.addAction(Action('Lorentzian'))
        
        #editor = ComboBox(parent)
        # editor.addItems(self.editorItems)
        # editor.setEditable(True)
        editor.setMenu(self.menu)
        #editor.currentIndexChanged.connect(self.currentIndexChanged)
        
        #editor.setCurrentIndex(-1)#self.editorItems.index(index.data()))
        #editor.addItems(self.editorItems)
        return editor

    '''
    def setEditorData(self,editor,index):
        z = 0
        for item in self.editorItems:
            #ai = QListWidgetItem(item)
            editor.addItem(item)
            if item == index.data():
                editor.setCurrentItem(editor.item(z))
            z += 1
        #editor.setGeometry(0,index.row()*self.height,self.width,self.height*len(self.editorItems))
    '''
    def setModelData(self, editor, model, index):
        print(self, editor, model, index)
        return None
        editorIndex=editor.currentIndex()
        text=self.editorItems[editorIndex]
        model.setData(index, text)
        # print '\t\t\t ...setModelData() 1', text

    @pyqtSlot()
    def currentIndexChanged(self): 
        self.commitData.emit(self.sender())
            

class StellarParametersWidget(AnalysisWidget):
    
    def __init__(self, session, callback=None, parent=None):
        super().__init__(parent=parent, session=session, title="Stellar Parameters")
        layout = QVBoxLayout(self.widget)
        layout.setContentsMargins(0, 0, 0, 0)
        
        self.primaryButton = PrimaryPushButton("Measure all equivalent widths")
        top_layout = QHBoxLayout()
        top_layout.addStretch(1)
        top_layout.addWidget(self.primaryButton)
        layout.addLayout(top_layout)        
        
        lr_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        # stellar params
        left_layout.setContentsMargins(10, 10, 10, 10)
        left_layout.setSpacing(10)
        
        line_edit_teff = LineEdit(self)
        line_edit_teff.setValidator(QIntValidator(2_500, 10_000, self))
        line_edit_teff.setText(f"{cfg.get(cfg.initial_teff):,}")
        line_edit_teff.setFixedWidth(80)
        line_edit_teff.setAlignment(Qt.AlignHCenter)

        line_edit_logg = LineEdit(self)
        line_edit_logg.setValidator(QDoubleValidator(0, 5, 3, self))
        line_edit_logg.setText("4.5")
        line_edit_logg.setFixedWidth(80)
        line_edit_logg.setAlignment(Qt.AlignHCenter)        
        
        line_edit_feh = LineEdit(self)
        line_edit_feh.setValidator(QDoubleValidator(-5, 1, 3, self))
        line_edit_feh.setText("0")
        line_edit_feh.setFixedWidth(80)
        line_edit_feh.setAlignment(Qt.AlignHCenter)        

        line_edit_v_micro = LineEdit(self)
        line_edit_v_micro.setValidator(QDoubleValidator(0, 10, 3, self))
        line_edit_v_micro.setText("1")
        line_edit_v_micro.setFixedWidth(80)
        line_edit_v_micro.setAlignment(Qt.AlignHCenter)
                
        param_layout = QGridLayout()
        param_layout.setColumnStretch(0, 1)
        param_layout.addWidget(QLabel("Effective temperature:"), 0, 0)
        param_layout.addWidget(line_edit_teff, 0, 2)
        param_layout.addWidget(QLabel("K"), 0, 3)
        
        param_layout.addWidget(QLabel("Surface gravity:"), 1, 0)
        param_layout.addWidget(line_edit_logg, 1, 2)
        param_layout.addWidget(QLabel("dex"), 1, 3)
        
        param_layout.addWidget(QLabel("Metallicity ([M/H]):"), 2, 0)
        param_layout.addWidget(line_edit_feh, 2, 2)
        param_layout.addWidget(QLabel("dex"), 2, 3)
        
        param_layout.addWidget(QLabel("Microturbulence:"), 3, 0)
        param_layout.addWidget(line_edit_v_micro, 3, 2)
        param_layout.addWidget(QLabel("km/s"), 3, 3)
        
        left_layout.addLayout(param_layout)
        
        
        
        #right_layout = QVBoxLayout()
        #right_layout.addWidget(ExcitationIonizationBalanceWidget(parent=self), 0, Qt.AlignTop)
        
        
        self.tableView = TableView(self)        
        #header = ["", "Wavelength", "Species", "χ", "log(gf)", "Γ", "Deg", "Tol", "Profile", "Automask", "EW", "σ(EW)", "A(X)", "σ(A(X))"]
        self.table_model = TransitionsTableModel(self.tableView, [
            [0, 5160, "Fe I", 4.22, -1.01, 0.003, 3, 0.1, "Gaussian", "0.1 ± 0.1", "7.5 ± 0.1"],
            [1, 5172, "Fe II", 4.1, -3, None, -1, 0.1, "Voight", "0.1 ± 0.1", "7.5 ± 0.1"],         
        ])
        self.tableView.setModel(self.table_model)
        self.tableView.scrollDelagate = SmoothScrollDelegate(self.tableView)
        self.tableView.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        
        #    attrs = ("is_acceptable", "_repr_wavelength", "_repr_element", 
        #    "equivalent_width")
        
        #self.tableView.setFixedSize(400, 800)

        # NOTE: use custom item delegate
        #self.tableView.setItemDelegate(CustomTableItemDelegate(self.tableView))
        self.tableView.setItemDelegateForColumn(8, ComboDelegate(self))



        # select row on right-click
        # self.tableView.setSelectRightClickedRow(True)

        # enable border
        self.tableView.setBorderVisible(True)
        self.tableView.setBorderRadius(8)

        self.tableView.setWordWrap(False)

        self.tableView.verticalHeader().hide()
        self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)        
        self.tableView.horizontalHeader().setSectionResizeMode(0, QHeaderView.Fixed)
        self.tableView.setColumnWidth(0, 45)        
        self.tableView.resizeColumnsToContents()
        
        self.tableView.horizontalHeader().setSectionsMovable(True)
        self.tableView.setSortingEnabled(True)
        self.tableView.setDragEnabled(True) 
        self.tableView.setTabKeyNavigation(True)
        #self.tableView.selectionModel().currentChanged.connect(self.line_current_changed)
        self.tableView.selectionModel().selectionChanged.connect(self.line_selection_changed)
        
        #        self.tableView.selectionModel().
        self.tableView.setMinimumSize(800, 200)
        self.tableView.setSelectionMode(QTableView.SingleSelection)
        #self.tableView.setSelectionBehavior(..)

        self.setStyleSheet("Demo{background: rgb(255, 255, 255)} ")
        #layout.setContentsMargins(50, 30, 50, 30)
        lr_layout.addLayout(left_layout)
        lr_layout.addStretch(1)
        #lr_layout.addLayout(right_layout)
        layout.addWidget(self.tableView)
        layout.addWidget(ExcitationIonizationBalanceWidget(parent=self))
        layout.addLayout(lr_layout)
                
        
        #layout.addWidget(QLabel(" " g gglshd))  
        

        #self.primaryButton = PrimaryPushButton("Measure all equivalent widths")
        #self.bottomLayout.addStretch(1)
        #self.bottomLayout.addWidget(self.primaryButton)
        self.bottomLayout.addWidget(PushButton("Compute abundances with these stellar parameters"))
        self.bottomLayout.addStretch(1)
        push_button_solve = PrimarySplitPushButton("Solve stellar parameters")
        
        #menu = RoundMenu()
        menu = CheckableMenu(parent=self, indicatorType=MenuIndicatorType.CHECK)
        menu.addAction(Action("Hold effective temperature constant"))
        menu.addAction(Action("Hold surface gravity constant"))
        menu.addAction(Action("Hold metallicity constant"))
        menu.addAction(Action("Hold microturbulence constant"))
        for i in range(3):            
            menu.actions()[i].setCheckable(True)
                
        push_button_solve.setFlyout(menu)
        
        self.bottomLayout.addWidget(push_button_solve)
        
        self.widget.setParent(self.card)
        self.widget.show()
        
        return None
    
    def line_current_changed(self):
        """Force the first column (checkbox) to always be in focus."""
        model = self.tableView.selectionModel()
        current_index = model.currentIndex()
        if current_index.column() > 0:
            model.setCurrentIndex(
                self.table_model.index(current_index.row(), 0), 
                QItemSelectionModel.Select
            )
        return None
    
    def line_selection_changed(self):        
        try:
            index = self.tableView.selectionModel().selectedRows()[0]
        except IndexError:
            return None
        else:
            print(f"selection changed {index}")
        


class FloatLineEdit(LineEdit):

    valueChanged = pyqtSignal(str)

    def __init__(self, value, lower, upper, decimals, parent=None):
        super().__init__(parent)
        self.setText(str(value))
        #self.setFixedSize(136, 33)
        self.setClearButtonEnabled(True)
        self.setValidator(QDoubleValidator(lower, upper, decimals, self))

        self.textEdited.connect(self._onTextEdited)

    def _onTextEdited(self, text):
        """ text edited slot """
        
        state = self.validator().validate(text, 0)[0]
        if state == QDoubleValidator.Acceptable:
            self.valueChanged.emit(text)
            print(f"ok {text}")
        else:
            print("nope")
            
        return False

    
class RadialVelocityWidget(AnalysisWidget):
    
    
    def __init__(self, session, callback=None, parent=None, stretch=0):
        
        # Add buttons to bottom layout
        super().__init__(parent=parent, session=session, title="Radial Velocity")
        self.session = session
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
        
        # Get the order closest to the wavelength range.
        index = session._get_closest_spectrum_index(8500)                
        wavelength, flux, *_ = session.spectra[index]
        
        
        layout = QVBoxLayout(self.widget)
        layout.setContentsMargins(0, 0, 0, 0)  
        
        #index = 28
        plot = SinglePlotWidget(
            x=wavelength,
            y=flux / np.nanmedian(flux),            
            xlabel=r"Wavelength [$\mathrm{\AA}$]",
            ylabel="Rectified flux",
            figsize=(100, 2), toolbar=False, toolbar_left_right=False)
        plot.ax.set_xlim(wavelength[[0, -1]])
        
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
                
        self.button_no_shift = PushButton("Spectrum is already at rest")
        self.button_no_shift.clicked.connect(self.on_no_shift_button_pushed)
        self.button_measure = PushButton("Measure")
        self.button_measure.clicked.connect(self.on_shift_button_clicked)
        self.button_measure_and_shift = PrimaryPushButton("Measure and shift")
        self.button_measure_and_shift.clicked.connect(self.on_measure_and_shift_button_clicked)
        self.bottomLayout.addWidget(self.button_no_shift)
        self.bottomLayout.addStretch(1)
        self.bottomLayout.addWidget(QLabel("Relative velocity:"))
        
        self.v_rel = LineEdit(self)
        self.v_rel.setValidator(QDoubleValidator(-1000, 1000, 3, self))
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
        
        state = self.v_rel.validator().validate(self.v_rel.text(), 0)[0]
        if state == QDoubleValidator.Acceptable:    
            self.button_measure.setVisible(True)
            self.button_measure_and_shift.setText("Shift")
        else:
            print("warning")
        print("edited")
        


class SessionInterface(ScrollArea):
    """ Dialog interface """

    def __init__(self, session, parent=None):
        super().__init__(parent=parent)
        
        self.session = session
        
        self.view = QWidget(self)
        self.toolBar = ToolBar(parent=self)
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

        initial_visibility = False

        def continuum_callback():
            self.cog_eqw = StellarParametersWidget(session, self)
            self.vBoxLayout.addWidget(self.cog_eqw, 0, Qt.AlignTop)
            

        # Radial velocity stuff
        self.continuum = ContinuumRectificationWidget(session, callback=continuum_callback, parent=self)
        self.continuum.setVisible(initial_visibility)
        
        
        def rv_callback():
            self.continuum.setVisible(True)
        
        card = RadialVelocityWidget(session, callback=rv_callback, parent=self)
        self.vBoxLayout.addWidget(card, 0, Qt.AlignTop)
        
        # Continuum normalization

        self.vBoxLayout.addWidget(self.continuum, 0, Qt.AlignTop)
        
        
        # add stellar parameter analysis? --> differential, etc.
        
        
        
                    
        
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

        session = Session([
            "/Users/andycasey/Downloads/hd122563blue_multi.fits",
            "/Users/andycasey/Downloads/hd122563red_multi.fits"
        ])

        self.myInterface = SessionInterface(session, self)        
        self.addMySubInterface(self.myInterface, 'myInterface', 'HD 122563')

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
            
            session = Session(filenames)
            
            # Get a suggested name.
            try:
                name = getval(filenames[0], "OBJECT", 0)
            except:
                name = f"Untitled-{self.tabCount}"
                
            
            
            # load the file, parse a name from it.        
            widget = SessionInterface(session, self)

            self.addMySubInterface(widget, name, name)        
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

