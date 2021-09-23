#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" GUIs related to exceptions and warnings. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
import os
import requests
import sys
import traceback as tb
#from tempfile import mkstemp
from smh.utils import mkstemp
from urllib.parse import quote
from PySide2 import (QtCore, QtGui as QtGui2, QtWidgets as QtGui)

#from smh import __git_status__

logger = logging.getLogger("smh")

def format(color, style=None):
    """
    Return a QTextCharFormat with the given attributes.

    :param color:
        The color to format the text.

    :param style: [optional]
        The styling for the text.

    """
    _color = QtGui2.QColor()
    _color.setNamedColor(color)

    _format = QtGui2.QTextCharFormat()
    _format.setForeground(_color)
    if style is not None and 'bold' in style:
        _format.setFontWeight(QtGui2.QFont.Bold)
    if style is not None and 'italic' in style:
        _format.setFontItalic(True)

    return _format


# Syntax styles that can be shared by all languages
STYLES = {
    'keyword': format('green'),
    'operator': format('red'),
    'brace': format('darkGray'),
    'defclass': format('black', 'bold'),
    'string': format('blue'),
    'string2': format('blue'),
    'comment': format('darkGreen', 'italic'),
    'self': format('black', 'italic'),
    'numbers': format('brown'),
}

class PythonHighlighter(QtGui2.QSyntaxHighlighter):
    """
    Syntax highlighter for the Python language.
    """
    # Python keywords
    keywords = [
        'and', 'assert', 'break', 'class', 'continue', 'def',
        'del', 'elif', 'else', 'except', 'exec', 'finally',
        'for', 'from', 'global', 'if', 'import', 'in',
        'is', 'lambda', 'not', 'or', 'pass', 'print',
        'raise', 'return', 'try', 'while', 'yield',
        'None', 'True', 'False',
    ]

    # Python operators
    operators = [
        '=',
        # Comparison
        '==', '!=', '<', '<=', '>', '>=',
        # Arithmetic
        '\+', '-', '\*', '/', '//', '\%', '\*\*',
        # In-place
        '\+=', '-=', '\*=', '/=', '\%=',
        # Bitwise
        '\^', '\|', '\&', '\~', '>>', '<<',
    ]

    # Python braces
    braces = [
        '\{', '\}', '\(', '\)', '\[', '\]',
    ]
    def __init__(self, document):
        QtGui2.QSyntaxHighlighter.__init__(self, document)

        # Multi-line strings (expression, flag, style)
        # FIXME: The triple-quotes in these two lines will mess up the
        # syntax highlighting from this point onward
        self.tri_single = (QtCore.QRegExp("'''"), 1, STYLES['string2'])
        self.tri_double = (QtCore.QRegExp('"""'), 2, STYLES['string2'])

        rules = []

        # Keyword, operator, and brace rules
        rules += [(r'\b%s\b' % w, 0, STYLES['keyword'])
            for w in PythonHighlighter.keywords]
        rules += [(r'%s' % o, 0, STYLES['operator'])
            for o in PythonHighlighter.operators]
        rules += [(r'%s' % b, 0, STYLES['brace'])
            for b in PythonHighlighter.braces]

        # All other rules
        rules += [
            # 'self'
            (r'\bself\b', 0, STYLES['self']),

            # Double-quoted string, possibly containing escape sequences
            (r'"[^"\\]*(\\.[^"\\]*)*"', 0, STYLES['string']),
            # Single-quoted string, possibly containing escape sequences
            (r"'[^'\\]*(\\.[^'\\]*)*'", 0, STYLES['string']),

            # 'def' followed by an identifier
            (r'\bdef\b\s*(\w+)', 1, STYLES['defclass']),
            # 'class' followed by an identifier
            (r'\bclass\b\s*(\w+)', 1, STYLES['defclass']),

            # From '#' until a newline
            (r'#[^\n]*', 0, STYLES['comment']),

            # Numeric literals
            (r'\b[+-]?[0-9]+[lL]?\b', 0, STYLES['numbers']),
            (r'\b[+-]?0[xX][0-9A-Fa-f]+[lL]?\b', 0, STYLES['numbers']),
            (r'\b[+-]?[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?\b', 0, STYLES['numbers']),
        ]

        # Build a QRegExp for each pattern
        self.rules = [(QtCore.QRegExp(pat), index, fmt)
            for (pat, index, fmt) in rules]


    def highlightBlock(self, text):
        """Apply syntax highlighting to the given block of text.
        """
        # Do other syntax formatting
        for expression, nth, format in self.rules:
            index = expression.indexIn(text, 0)

            while index >= 0:
                # We actually want the index of the nth match
                index = expression.pos(nth)
                length = len(expression.cap(nth))
                self.setFormat(index, length, format)
                index = expression.indexIn(text, index + length)

        self.setCurrentBlockState(0)

        # Do multi-line strings
        in_multiline = self.match_multiline(text, *self.tri_single)
        if not in_multiline:
            in_multiline = self.match_multiline(text, *self.tri_double)


    def match_multiline(self, text, delimiter, in_state, style):
        """Do highlighting of multi-line strings. ``delimiter`` should be a
        ``QRegExp`` for triple-single-quotes or triple-double-quotes, and
        ``in_state`` should be a unique integer to represent the corresponding
        state changes when inside those strings. Returns True if we're still
        inside a multi-line string when this function is finished.
        """
        # If inside triple-single quotes, start at 0
        if self.previousBlockState() == in_state:
            start = 0
            add = 0
        # Otherwise, look for the delimiter on this line
        else:
            start = delimiter.indexIn(text)
            # Move past this match
            add = delimiter.matchedLength()

        # As long as there's a delimiter match on this line...
        while start >= 0:
            # Look for the ending delimiter
            end = delimiter.indexIn(text, start + add)
            # Ending delimiter on this line?
            if end >= add:
                length = end - start + add + delimiter.matchedLength()
                self.setCurrentBlockState(0)
            # No; multi-line string
            else:
                self.setCurrentBlockState(in_state)
                length = len(text) - start + add
            # Apply formatting
            self.setFormat(start, length, style)
            # Look for the next match
            start = delimiter.indexIn(text, start + length)

        # Return True if still inside a multi-line string, False otherwise
        if self.currentBlockState() == in_state:
            return True
        else:
            return False



class ExceptionWidget(QtGui.QDialog):
    def __init__(self, exception_type, message, traceback, *args):
        """
        Initialize a widget to display details about an exception that was
        raised.


        :param exception_type:
            The type of exception that was raised.

        :param message:
            The exception message.

        :param traceback:
            The traceback object of the exception.
        """

        super(ExceptionWidget, self).__init__(*args)

        self.message = message
        self.exception_type = exception_type
        self.traceback = traceback

        self.setGeometry(600, 400, 600, 400)
        desktop = QtGui.QApplication.desktop()
        self.move(desktop.screen().rect().center() \
            - self.rect().center())
        self.setWindowTitle("Exception")

        sp = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, 
            QtGui.QSizePolicy.MinimumExpanding)
        sp.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sp)


        vertical_layout = QtGui.QVBoxLayout(self)

        # Label to explain what an exception is.
        label = QtGui.QLabel(self)
        label.setText(
            "An unexpected exception has occurred. The traceback and detailed "
            "information is available below.")
        vertical_layout.addWidget(label)

        # Text browser to display the exception details.
        self.show_traceback = QtGui.QPlainTextEdit(self)
        self.show_traceback.setFont(QtGui2.QFont("Courier", 12))
        self.show_traceback.setReadOnly(True)
        self.show_traceback.setPlainText("\n" + \
            "\n".join(tb.format_exception(exception_type, message, traceback)))

        self.highlighter = PythonHighlighter(self.show_traceback.document())
        

        vertical_layout.addWidget(self.show_traceback)

        # Horizontal layout with buttons.
        horizontal_layout = QtGui.QHBoxLayout()
        self.btn_ignore = QtGui.QPushButton(self)
        self.btn_ignore.setText("Ignore exceptions of this kind")
        horizontal_layout.addWidget(self.btn_ignore)

        self.btn_create_issue = QtGui.QPushButton(self)
        self.btn_create_issue.setText("Create GitHub issue..")
        horizontal_layout.addWidget(self.btn_create_issue)

        # Spacer with a minimum width.
        horizontal_layout.addItem(QtGui.QSpacerItem(40, 20, 
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))

        self.btn_ok = QtGui.QPushButton(self)
        self.btn_ok.setText("OK")
        self.btn_ok.setDefault(True)
        self.btn_ok.setFocus()
        horizontal_layout.addWidget(self.btn_ok)

        vertical_layout.addLayout(horizontal_layout)

        # Connect the buttons.
        self.btn_ok.clicked.connect(self.close)
        self.btn_create_issue.clicked.connect(self.create_github_issue)
        self.btn_ignore.clicked.connect(self.ignore_exceptions_of_this_kind)

        self.ignore_in_future = False
        return None


    def ignore_exceptions_of_this_kind(self):
        """
        Trigger for when the 'Ignore exceptions of this kind' button has been
        pressed.
        """

        self.ignore_in_future = True
        self.close()
        return None


    def create_github_issue(self):
        """
        Trigger for when the 'Create GitHub issue' button has been pushed.
        """

        self.close()

        # Get the parent and session etc.
        parent = QtCore.QCoreApplication.instance().window
        session = parent.session

        # Save and upload the session.
        # TODO: Once we can save sessions..

        # Get a screenshot of the main window and upload it.
        screenshot_url = None
        logger.info("Taking and uploading screenshot..")
        try:
            screenshot = QtGui.QPixmap.grabWindow(parent.winId())
            _, path = mkstemp(suffix=".png")
            screenshot.save(path, "png")

            response = requests.put(
                url="https://transfer.sh/screenshot.png",
                data=open(path, "rb"))

            if response.status_code == 200:
                screenshot_url = response.text.strip()  
            else:
                logger.warn("Screenshot could not be uploaded (response {})"\
                    .format(response.status_code))

            os.remove(path)

        except:
            logger.exception("No screenshot could be uploaded.")

        # Save session.
        include_session_str = ""
        logger.info("Uploading session..")
        try:
            _, path = mkstemp(suffix=".smh")
            session.save(path, overwrite=True)

            response = requests.put(
                url="https://transfer.sh/session.smh",
                data=open(path, "rb"))

            if response.status_code == 200:
                include_session_str = "in [this session]({}) ".format(
                    response.text.strip())

            else:
                logger.warn("Session could not be uploaded (response {})".format(
                        response.status_code))

            os.remove(path)

        except:
            logger.exception("Session could not be uploaded.")

        logger.info("Creating issue..")
        body_template = \
            "An exception was encountered {include_session_str}"\
            "using version {application_version} on Python {sys_version}:\n\n"\
            "````python\n"\
            "{formatted_exception}"\
            "````\n\n"\
            "{screenshot_str}"\
            "**Additional details on how this exception occurred**:\n"

        # Create text.
        title = quote("Exception raised: {}".format(self.message))
        body = quote(body_template.format(
            include_session_str=include_session_str,
            application_version=__git_status__,
            sys_version=sys.version.replace("\n", ""),
            formatted_exception="\n".join(tb.format_exception(
                self.exception_type, self.message, self.traceback)),
            screenshot_str="" if screenshot_url is None else 
                "Below is a screenshot from the application at the time that "\
                "the exception occurred:\n\n"\
                "![screenshot]({})\n\n".format(screenshot_url)))

        url = "https://github.com/andycasey/smhr/issues/new?title={}&body={}"\
            .format(title, body)

        os.system('python -m webbrowser "{}"'.format(url))
        return None



if sys.platform == "darwin":
        
    # See http://successfulsoftware.net/2013/10/23/fixing-qt-4-for-mac-os-x-10-9-mavericks/
    substitutes = [
        (".Lucida Grande UI", "Lucida Grande"),
        (".Helvetica Neue DeskInterface", "Helvetica Neue")
    ]
    for substitute in substitutes:
        QtGui2.QFont.insertSubstitution(*substitute)


if __name__ == "__main__":

    # This is just for development testing.
    #app = QtGui.QApplication(sys.argv)
    window = ExceptionWidget("something", "message", "traceback")
    window.exec_()

    
