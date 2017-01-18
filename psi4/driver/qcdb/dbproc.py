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

r"""File to

"""
from __future__ import absolute_import
from __future__ import print_function
import sys
import os
import glob
import ast


def useful():
    print("in qcdb.useful()")
    return 'qcdb successfully accessed'


def drop_duplicates(seq):
    """Function that given an array or array of arrays *seq*, returns an
    array without any duplicate entries. There is no guarantee of which
    duplicate entry is dropped.

    """
    noDupes = []
    seq2 = sum(seq, [])
    [noDupes.append(i) for i in seq2 if not noDupes.count(i)]
    return noDupes


def dictify_database_docstrings():
    """

    """
    db_path = os.path.dirname(__file__) + '/../databases'

    DSD = {}
    module_choices = []
    for module in glob.glob(db_path + '/*.py'):
        filename = os.path.split(module)[1]
        basename = os.path.splitext(filename)[0]
        div = '=' * len(basename)

        module_choices.append(basename)
        DSD[basename] = {}

        M = ast.parse(''.join(open(module)))
        DS = ast.get_docstring(M)
        if not DS:
            DS = ""
        DS = str.replace(DS, '|dl|', '-->')
        DS = str.replace(DS, '|dr|', '<--')
        DS = str.replace(DS, "``'", '')
        DS = str.replace(DS, "'``", '')

        lst = DS.split("\n- **")

        #DSD[basename]['general'] = str.replace(lst[0], '|', '')
        DSD[basename]['general'] = lst[0].split('\n')

        try:
            DSD[basename]['cp'] = [section for section in lst if section.startswith("cp")][0]
        except IndexError:
            DSD[basename]['cp'] = None

        try:
            DSD[basename]['rlxd'] = [section for section in lst if section.startswith("rlxd")][0]
        except IndexError:
            DSD[basename]['rlxd'] = None

        try:
            DSD[basename]['benchmark'] = [section for section in lst if section.startswith("benchmark")][0]
        except IndexError:
            DSD[basename]['benchmark'] = None

        try:
            #DSD[basename]['subset'] = [section for section in lst if section.startswith("subset")][0]
            temp = [section for section in lst if section.startswith("subset")][0].splitlines()
            temp = temp[2:]

            result = {}
            for item in temp:
                item = item.lstrip(" -")
                try:
                    key, val = item.split(" ", 1)
                    result[key] = val
                except ValueError:
                    result[item] = ""

            DSD[basename]['subset'] = result

        except IndexError:
            DSD[basename]['subset'] = {"": 'No subsets available'}

    return DSD

    #    print '\ngeneral\n\n', DSD[basename]['general']
    #    print '\ncp\n\n', DSD[basename]['cp']
    #    print '\nrlxd\n\n', DSD[basename]['rlxd']
    #    print '\nbenchmark\n\n', DSD[basename]['benchmark']
    #    print '\nsubset\n\n', DSD[basename]['subset']

        #print '  %-12s   %s' % ('[' + basename + ']', DSD[basename]['general'][0])

    #print 'DSD2\n', DSD['S22']['subset']
