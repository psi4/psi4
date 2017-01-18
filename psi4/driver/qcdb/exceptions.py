#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

"""Module with non-generic exceptions classes."""
from __future__ import print_function


class QcdbException(Exception):
    """Error class for QCDB."""
    pass


class FeatureNotImplemented(QcdbException):
    """Error called for functions defined but not yet implemented.
    Also for functions defined that will never be implemented.

    """
    def __init__(self, msg):
        QcdbException.__init__(self, msg)
        self.msg = msg
        print('\nQcdbException: Feature %s is not yet implemented.\n\n' % (msg))


class ValidationError(QcdbException):
    """Error called for problems with syntax input file. Prints
    error message *msg* to standard output stream.

    """
    def __init__(self, msg):
        QcdbException.__init__(self, msg)
        self.msg = msg
        print('\nQcdbException: %s\n\n' % (msg))


class IncompleteAtomError(QcdbException):
    """Error raised when not all variables in an atom specification
    have been defined at compute time. May be a temporary situation
    so message not printed but appears as traceback when error persists.

    """
    def __init__(self, msg):
        QcdbException.__init__(self, msg)
        self.msg = msg


class ParsingValidationError(QcdbException):
    """Error called for problems with syntax from a QC output file. Prints
    error message *msg* to standard output stream.

    """
    def __init__(self, msg):
        QcdbException.__init__(self, msg)
        self.msg = msg
        print('\nQcdbException: %s\n\n' % (msg))


class FragmentCountError(QcdbException):
    """Error called molecule has wrong number of fragments for method.
    Prints error message *msg* to standard output stream.

    """
    def __init__(self, msg):
        QcdbException.__init__(self, msg)
        self.msg = msg
        #print('\nQcdbException: %s\n\n' % (msg))


class BasisSetFileNotFound(QcdbException):
    """

    """
    def __init__(self, msg):
        QcdbException.__init__(self, msg)
        self.msg = msg
        print('\nQcdbException BasisSetFileNotFound: %s\n\n' % (msg))


class BasisSetNotFound(QcdbException):
    """

    """
    def __init__(self, msg, silent=False):
        QcdbException.__init__(self, msg)
        self.msg = msg
        if not silent:
            print('\nQcdbException BasisSetNotFound: %s\n\n' % (msg))


class BasisSetNotDefined(QcdbException):
    """

    """
    def __init__(self, msg):
        QcdbException.__init__(self, msg)
        self.msg = msg
        print('\nQcdbException BasisSetNotDefined: %s\n\n' % (msg))


class Dftd3Error(QcdbException):
    """

    """
    def __init__(self, msg):
        QcdbException.__init__(self, msg)
        self.msg = msg
        print('\nDftd3Error: %s\n\n' % (msg))
