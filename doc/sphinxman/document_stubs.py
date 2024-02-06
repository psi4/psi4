#!/usr/bin/python

#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
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

import glob
import os
import re
import sys
from pathlib import Path

DriverPath = ''
if (len(sys.argv) == 2):
    DriverPath = sys.argv[1] + '/'
    sys.path.insert(0, os.path.abspath(os.getcwd()))


def pts(category, pyfile):
    print('Auto-documenting %s file %s' % (category, pyfile))

Path("source/api/").mkdir(parents=True, exist_ok=True)

for stub in ["psi4.core.del_variable",
             "psi4.core.get_array_variable",
             "psi4.core.get_array_variables",
             "psi4.core.get_variable",
             "psi4.core.get_variables",
             "psi4.core.has_variable",
             "psi4.core.set_global_option_python",
             "psi4.core.set_variable",
             "psi4.core.variable",
             "psi4.core.variables",
    ]:
    curmod = ".".join(stub.split(".")[:-1])
    basename = stub.split(".")[-1]
    div = '=' * len(basename)

    pts('stub', basename)

    with open(f"source/api/{stub}.rst", "w") as fp:
        fp.write(basename + "\n")
        fp.write(div + "\n\n")
        fp.write('.. currentmodule:: %s\n\n' % (curmod))
        fp.write('.. autofunction:: %s\n\n' % (basename))

