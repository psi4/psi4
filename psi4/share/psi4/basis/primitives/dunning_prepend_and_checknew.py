#!/usr/bin/env python

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

from __future__ import print_function

import os
import subprocess
import sys

qcdb_module = os.path.normpath(os.path.dirname(os.path.abspath(__file__)) + '../../../../../driver')
sys.path.append(qcdb_module)
import qcdb
from qcdb.libmintsbasissetparser import Gaussian94BasisSetParser

output = subprocess.check_output("ls -1 *cc-*.gbs | grep -v 'autogen' | grep -v 'tight' | grep -v 'polarization' | grep -v 'molpro' | grep -v 'diffuse' | grep -v 'basis' | grep -v 'corevalence' | grep -v 'hold'", shell=True)
real_dunnings = output.decode().split('\n')

parser = Gaussian94BasisSetParser()
os.system("echo '#differing basis sets' > basisdunningfiles.txt")

for bfl in real_dunnings:
    if not bfl:
        continue

    with open(bfl, 'r') as basfile:
        bascontents = basfile.readlines()

    if bascontents[0] != "spherical\n":
        print('{:30} {:4}'.format(bfl, 'sph '))
        with open(bfl, 'w') as basfile:
            basfile.write("spherical\n\n")
            for ln in bascontents:
                basfile.write(ln)
    else:
        print('{:30} {:4}'.format(bfl, ''))

for bfl in sorted(real_dunnings, key=lambda v: v[::-1], reverse=True):
    if not bfl:
        continue

    os.system('./diff_gbs.py {} ../{}'.format(bfl, bfl))
