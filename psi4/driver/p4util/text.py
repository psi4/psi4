#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
"""Module with utility classes and functions related
to data tables and text.

"""
import sys
import warnings

from psi4 import core
from psi4.driver import constants


class Table(object):
    """Class defining a flexible Table object for storing data."""

    def __init__(self, rows=(), row_label_width=10, row_label_precision=4, cols=(), width=16, precision=10):
        warnings.warn(
            "Using `psi4.driver.p4util.Table` is deprecated, and in 1.4 it will stop working\n",
            category=FutureWarning,
            stacklevel=2)
        self.row_label_width = row_label_width
        self.row_label_precision = row_label_precision
        self.width = width
        self.precision = precision
        self.rows = rows

        if isinstance(cols, str):
            self.cols = (cols, )
        else:
            self.cols = cols

        self.labels = []
        self.data = []

    def format_label(self):
        """Function to pad the width of Table object labels."""
        #str = lambda x: (('%%%d.%df' % (self.row_label_width, self.row_label_precision)) % x)
        str = lambda x: (('%%%ds' % (self.row_label_width)) % x)
        return " ".join(map(str, self.labels))

    def format_values(self, values):
        """Function to pad the width of Table object data cells."""
        str = lambda x: (('%%%d.%df' % (self.width, self.precision)) % x)
        return " ".join(map(str, values))

    def __getitem__(self, value):
        self.labels.append(value)
        return self

    def __setitem__(self, name, value):
        self.labels.append(name)
        label = self.format_label()
        self.labels = []

        if isinstance(value, list):
            self.data.append((label, value))
        else:
            self.data.append((label, [value]))

    def save(self, file):
        """Function to save string of the Table object to *file*."""
        import pickle
        pickle_str = pickle.dumps(self)
        fileobj = open(file, "w")
        fileobj.write(str(self))
        fileobj.close()

    def __str__(self):
        rowstr = lambda x: '%%%ds' % self.row_label_width % x
        colstr = lambda x: '%%%ds' % self.width % x

        lines = []

        table_header = ""
        if isinstance(self.rows, str):
            table_header += "%%%ds" % self.row_label_width % self.rows
        else:
            table_header += " ".join(map(rowstr, self.rows))
        table_header += " ".join(map(colstr, self.cols))

        lines.append(table_header)

        for datarow in self.data:
            #print datarow
            row_data = datarow[0]
            row_data += self.format_values(datarow[1])
            lines.append(row_data)

        return "\n".join(lines) + "\n"

    def copy(self):
        """Function to return a copy of the Table object."""
        import copy
        return copy.deepcopy(self)

    def absolute_to_relative(self, Factor=constants.hartree2kcalmol):
        """Function to shift the data of each column of the Table object
        such that the lowest value is zero. A scaling factor of *Factor* is applied.

        """
        import copy

        if len(self.data) == 0:
            return

        current_min = list(copy.deepcopy(self.data[0][1]))
        for datarow in self.data:
            for col in range(0, len(datarow[1])):
                if current_min[col] > datarow[1][col]:
                    current_min[col] = datarow[1][col]

        for datarow in self.data:
            for col in range(0, len(datarow[1])):
                #print datarow[1][col]
                datarow[1][col] = (datarow[1][col] - current_min[col]) * Factor

    def scale(self, Factor=constants.hartree2kcalmol):
        """Function to apply a scaling factor *Factor* to the
        data of the Table object.

        """
        if len(self.data) == 0:
            return

        for datarow in self.data:
            for col in range(0, len(datarow[1])):
                #print datarow[1][col]
                datarow[1][col] = datarow[1][col] * Factor


def banner(text, type=1, width=35, strNotOutfile=False):
    """Function to print *text* to output file in a banner of
    minimum width *width* and minimum three-line height for
    *type* = 1 or one-line height for *type* = 2. If *strNotOutfile*
    is True, function returns string rather than printing it
    to output file.

    """
    lines = text.split('\n')
    max_length = 0
    for line in lines:
        max_length = max(len(line), max_length)

    max_length = max(width, max_length)

    null = ''
    if type == 1:
        banner = '  //' + null.center(max_length, '>') + '//\n'
        for line in lines:
            banner += '  //' + line.center(max_length) + '//\n'
        banner += '  //' + null.center(max_length, '<') + '//\n'

    if type == 2:
        banner = ''
        for line in lines:
            banner += (' ' + line + ' ').center(max_length, '=')

    if strNotOutfile:
        return banner
    else:
        core.print_out(banner)


def print_stdout(stuff):
    """Function to print *stuff* to standard output stream."""
    warnings.warn(
        "Using `psi4.driver.p4util.print_stdout` instead of `print` is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)

    print(stuff, file=sys.stdout)


def print_stderr(stuff):
    """Function to print *stuff* to standard error stream."""
    warnings.warn(
        "Using `psi4.driver.p4util.print_stderr` instead of `print(..., file=sys.stderr)` is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)

    print(stuff, file=sys.stderr)


def levenshtein(seq1, seq2):
    """Compute the Levenshtein distance between two strings."""

    oneago = None
    thisrow = list(range(1, len(seq2) + 1)) + [0]
    for x in range(len(seq1)):
        twoago, oneago, thisrow = oneago, thisrow, [0] * len(seq2) + [x + 1]
        for y in range(len(seq2)):
            delcost = oneago[y] + 1
            addcost = thisrow[y - 1] + 1
            subcost = oneago[y - 1] + (seq1[x] != seq2[y])
            thisrow[y] = min(delcost, addcost, subcost)
    return thisrow[len(seq2) - 1]


def find_approximate_string_matches(seq1, options, max_distance):
    """Find list of approximate (within `max_distance`) matches to string `seq1` among `options`."""

    return [seq2 for seq2 in options if (levenshtein(seq1, seq2) <= max_distance)]
