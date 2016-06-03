#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
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
from __future__ import absolute_import
import psi4


class PsiException(Exception):
    """Error class for Psi."""
    pass


class ValidationError(PsiException):
    """Error called for problems with the input file. Prints
    error message *msg* to standard output stream and output file.

    """
    def __init__(self, msg):
        PsiException.__init__(self, msg)
        self.message = msg
        psi4.print_out('\nPsiException: %s\n\n' % (msg))

class ParsingError(PsiException):
    """Error called for problems parsing a text file. Prints error message
    *msg* to standard output stream and output file.

    """
    def __init__(self, msg):
        PsiException.__init__(self,msg)
        self.message = msg
        psi4.print_out('\nPsiException: %s\n\n' % (msg))

class TestComparisonError(PsiException):
    """Error called when a test case fails due to a failed
    compare_values() call. Prints error message *msg* to standard
    output stream and output file.

    """
    def __init__(self, msg):
        PsiException.__init__(self, msg)
        self.message = msg
        psi4.print_out('\nPsiException: %s\n\n' % (msg))


class CSXError(PsiException):
    """Error called when CSX generation fails.

    """
    def __init__(self, msg):
        PsiException.__init__(self, msg)
        self.message = msg
        psi4.print_out('\nCSXException: %s\n\n' % (msg))


class ManagedMethodError(PsiException):
    def __init__(self, circs):
        if circs[5] == '':
            msg = """{0}: Method '{1}' with {2} '{3}' and REFERENCE '{4}' not available{5}""".format(*circs)
        else:
            msg = """{0}: Method '{1}' with {2} '{3}' and REFERENCE '{4}' not directable to QC_MODULE '{5}'""".format(*circs)
        PsiException.__init__(self, msg)
        self.message = msg
        psi4.print_out('\nPsiException: %s\n\n' % (msg))


class Dftd3Error(PsiException):
    """

    """
    def __init__(self, msg):
        PsiException.__init__(self, msg)
        self.message = msg
        print('\nDftd3Error: %s\n\n' % (msg))
