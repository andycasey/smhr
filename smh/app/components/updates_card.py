# coding:utf-8
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QFrame, QLabel, QVBoxLayout, QHBoxLayout, QSizePolicy

from qfluentwidgets import IconWidget, TextWrap, FlowLayout, CardWidget
from ..common.signal_bus import signalBus
from ..common.style_sheet import StyleSheet


class UpdatesCard(CardWidget):
    """ Updates card """

    def __init__(self, icon, title, content, routeKey, index, parent=None):
        super().__init__(parent=parent)
        self.index = index
        self.routekey = routeKey

        self.iconWidget = IconWidget(icon, self)
        self.titleLabel = QLabel(title, self)
        # TODO: 120 here is magic
        self.contentLabel = QLabel(TextWrap.wrap(content, 120, False)[0], self)

        self.hBoxLayout = QHBoxLayout(self)
        self.vBoxLayout = QVBoxLayout()
        self.iconLayout = QVBoxLayout()
        self.iconLayout.setSpacing(28)
        self.iconLayout.setContentsMargins(0, 20, 0, 20)

        # TODO: I don't know how to set this to expand to the full horizontal space, but it would be nicer if it did that.
        self.setFixedWidth(830)

        self.iconWidget.setFixedSize(48, 48)

        self.hBoxLayout.setSpacing(28)
        self.hBoxLayout.setContentsMargins(20, 0, 0, 0)
        self.vBoxLayout.setSpacing(2)
        self.vBoxLayout.setContentsMargins(0, 20, 0, 20)
        self.vBoxLayout.setAlignment(Qt.AlignVCenter)

        # align to vertical top
        self.hBoxLayout.setAlignment(Qt.AlignTop)
        
        self.iconLayout.addWidget(self.iconWidget, alignment=Qt.AlignTop)
        self.iconLayout.addStretch(1)
        #self.hBoxLayout.addWidget(self.iconWidget)
        self.hBoxLayout.addLayout(self.iconLayout)
        self.hBoxLayout.addLayout(self.vBoxLayout)
        
        self.vBoxLayout.addStretch(1)
        self.vBoxLayout.addWidget(self.titleLabel)
        self.vBoxLayout.addWidget(self.contentLabel)
        self.vBoxLayout.addStretch(1)

        self.titleLabel.setObjectName('titleLabel')
        self.contentLabel.setObjectName('contentLabel')

    def mouseReleaseEvent(self, e):
        super().mouseReleaseEvent(e)
        signalBus.switchToSampleCard.emit(self.routekey, self.index)


class UpdatesCardView(QWidget):
    """ Updates card view """

    def __init__(self, title: str, parent=None):
        super().__init__(parent=parent)
        self.titleLabel = QLabel(title, self)
        self.vBoxLayout = QVBoxLayout(self)
        self.flowLayout = FlowLayout()


        self.vBoxLayout.setContentsMargins(36, 0, 36, 0)
        self.vBoxLayout.setSpacing(10)
        self.flowLayout.setContentsMargins(0, 0, 0, 0)
        self.flowLayout.setHorizontalSpacing(12)
        self.flowLayout.setVerticalSpacing(12)

        self.vBoxLayout.addWidget(self.titleLabel)
        self.vBoxLayout.addLayout(self.flowLayout, 1)

        self.titleLabel.setObjectName('viewTitleLabel')
        StyleSheet.SAMPLE_CARD.apply(self)

    def addSampleCard(self, icon, title, content, routeKey, index):
        """ add sample card """
        card = UpdatesCard(icon, title, content, routeKey, index, self)        
        self.flowLayout.addWidget(card)
