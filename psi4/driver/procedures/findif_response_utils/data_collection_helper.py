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

"""
Module of helper functions for distributed ccresponse computations.

Defines functions for retrieving data computed at displaced geometries.
"""
from __future__ import absolute_import
from __future__ import print_function
import collections
import shelve
import copy
import os

from psi4.driver import p4util
from psi4.driver.constants import *


def collect_displaced_matrix_data(db, signature, row_dim):
    """
        Gathers a list of tensors, one at each displaced geometry.

    db: (database) the database object for this property calculation
    signature: (string) The string that notifies the matrix reader that the
        targeted tensor data begins.
    row_dim: the expected number of rows that this value should be printed
        across in the file

    Returns a 2d list result[i][j]:
        i: indexes displacements
        j: indexes elements of the flattened tensor at some displacement
    Throws: none
    """
    result = []
    for job in db['job_status']:
        with open('{}/output.dat'.format(job)) as outfile:
            result.append(parse_geometry_matrix_data(outfile, signature, row_dim))

    return result

    # END collect_displaced_matrix_data()


def parse_geometry_matrix_data(outfile, matrix_name, row_tot):
    """
        Parses data from a 3 by n  matrix printed to a file

    outfile: ( file ) handle open in read mode, where the data should be found
    matrix_name: ( string ) that indicates the matrix data is found on the lines
        below
    row_tot: ( int ) indicates the number of lines that the matrix data should
        be printed across in the file

    Returns: matrix_data a list of matrix elements, len = 3*row_tot

    Throws: ParsingError (Collecting matrix data failed) if
            It can't find matrix_header in the file.
            It found matrix_header, but no data.
            It found matrix_header, and data but the number of elements is
            incorrect.

    """
    collect_matrix = False
    n_rows = 0
    n_tries = 0
    matrix_data = []
    for line in outfile:
        if matrix_name in line:
            collect_matrix = True
        if collect_matrix and (n_rows < row_tot):
            try:
                n_tries += 1
                if n_tries > (row_tot + 13):
                    raise ParsingError('{} Matrix was unreadable. Scanned {}'
                                    'lines.'.format(matrix_name, n_tries))
                else:
                    (index, x, y, z) = line.split()
                    matrix_data.append(float(x))
                    matrix_data.append(float(y))
                    matrix_data.append(float(z))
                    n_rows += 1
            except:
                pass
        if (n_rows == row_tot) and (len(matrix_data) != 3 * row_tot):
            raise p4util.ParsingError('Collecting {} data failed!'
                            '\nExpected {} elements but only captured {}'.format(
                                matrix_name, 3 * row_tot, len(matrix_data)))
        if len(matrix_data) == 3 * row_tot:
            return matrix_data

    raise p4util.ParsingError('data for {}  was not found in the output file, '
                    'but it was marked for collection. Check output files '
                    'in displacement sub-dirs!'.format(matrix_name))

    # END parse_geometry_matrix_data()
