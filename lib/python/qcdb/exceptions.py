#
#@BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
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
#@END LICENSE
#

"""Module with non-generic exceptions classes."""


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
