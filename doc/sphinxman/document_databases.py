#!/usr/bin/python

#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2022 The Psi4 Developers.
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
        fdriver.write('.. automodule:: %s\n' % (basename))
        fdriver.write('   :noindex:\n\n')
        fdriver.write('----\n')

    fdriver.write('\n')
fdriver.close()
