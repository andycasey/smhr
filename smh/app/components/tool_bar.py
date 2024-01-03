from PyQt5.QtCore import Qt, QUrl
from PyQt5.QtGui import QColor, QDesktopServices
from PyQt5.QtWidgets import QHBoxLayout, QVBoxLayout, QWidget
from qfluentwidgets import CaptionLabel, FluentIcon, PushButton, TitleLabel

from ..common.config import EXAMPLE_URL, HELP_URL
from ..view.gallery_interface import SeparatorWidget


class ToolBar(QWidget):
    """ Tool bar """

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
