#!/usr/bin/python

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

import sys
import os
import glob
import re


DriverPath = ''
if (len(sys.argv) == 2):
    DriverPath = sys.argv[1] + '/'
    sys.path.insert(0, os.path.abspath(os.getcwd()))


def pts(category, pyfile):
    print('Auto-documenting %s file %s' % (category, pyfile))


# Available databases in psi4/share/psi4/databases
fdriver = open('source/autodoc_available_databases.rst', 'w')
fdriver.write('\n\n')

for pyfile in glob.glob(DriverPath + '../../psi4/share/psi4/databases/*.py'):
    filename = os.path.split(pyfile)[1]
    basename = os.path.splitext(filename)[0]
    div = '=' * len(basename)

    if basename not in ['input']:

        pts('database', basename)

        fdriver.write(':srcdb:`%s`\n%s\n\n' % (basename, '"' * (9 + len(basename))))
        fdriver.write('.. automodule:: %s\n\n' % (basename))
        fdriver.write('----\n')

    fdriver.write('\n')
fdriver.close()
