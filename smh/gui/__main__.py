
#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Main GUI window for SMH """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

from PySide import QtCore, QtGui


class MainWindow(QtGui.QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        widget = QtGui.QWidget()
        self.setCentralWidget(widget)

        topFiller = QtGui.QWidget()
        topFiller.setSizePolicy(QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Expanding)

        self.context = QtGui.QLabel("What up?",
                alignment=QtCore.Qt.AlignCenter)
        self.context.setFrameStyle(
            QtGui.QFrame.StyledPanel | QtGui.QFrame.Sunken)

        bottomFiller = QtGui.QWidget()
        bottomFiller.setSizePolicy(QtGui.QSizePolicy.Expanding,
                QtGui.QSizePolicy.Expanding)

        vbox = QtGui.QVBoxLayout()
        vbox.setContentsMargins(5, 5, 5, 5)
        vbox.addWidget(topFiller)
        vbox.addWidget(self.context)
        vbox.addWidget(bottomFiller)
        widget.setLayout(vbox)

        # Create actions and menus.
        self.__init_menus__()

        self.statusBar().showMessage("")

        self.setWindowTitle("SMHr")
        self.setMinimumSize(160,160)
        self.resize(480,320)



    def new_session(self):
        """ Initialise new session. """
        print("New session")
        return None


    def open_session(self):
        """ Open existing session. """
        print("Open session")
        return None


    def save_session(self):
        """ Save session. """
        print("Save session")
        return None

    def save_session_as(self):
        """ Save session as new filename. """
        print("Save session as")
        return None


    def about(self):
        QtGui.QMessageBox.about(self, "About Menu",
                """
                SMHr
                
                Gotta pay back that tech debt.
                """)


    def __init_menus__(self):
        """ Initialize main window menus and associated actions. """

        # File menu actions.
        new_session = QtGui.QAction("&New", self,
                shortcut=QtGui.QKeySequence.New,
                statusTip="Create a new session",
                triggered=self.new_session)

        open_session = QtGui.QAction("&Open...", self,
                shortcut=QtGui.QKeySequence.Open,
                statusTip="Open an existing session",
                triggered=self.open_session)

        save_session = QtGui.QAction("&Save", self,
                shortcut=QtGui.QKeySequence.Save,
                statusTip="Save the session to disk",
                triggered=self.save_session)

        save_session_as = QtGui.QAction("Save &As", self,
            statusTip="Save the session to a new file",
            triggered=self.save_session_as)

        # Help menu actions.
        about = QtGui.QAction("&About", self,
                statusTip="Show the application's About box",
                triggered=self.about)

        # File menu.
        self.fileMenu = self.menuBar().addMenu("&File")
        self.fileMenu.addAction(new_session)
        self.fileMenu.addAction(open_session)
        self.fileMenu.addAction(save_session)
        self.fileMenu.addAction(save_session_as)

        # Help menu.
        self.helpMenu = self.menuBar().addMenu("&Help")
        self.helpMenu.addAction(about)




if __name__ == '__main__':

    import sys

    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())

