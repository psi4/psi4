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

"""Run on a single gbs file to get a periodic table printout or on two to compare gbs contents."""

from __future__ import print_function

import os
import subprocess
import sys

qcdb_module = os.path.normpath(os.path.dirname(os.path.abspath(__file__)) + '../../../../../driver')
sys.path.append(qcdb_module)
import qcdb
import qcelemental as qcel
from qcdb.libmintsbasissetparser import Gaussian94BasisSetParser


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def bas_sanitize(fl):
    if fl[-4] == '.gbs':
        fl = fl[:-4]
    return fl.lower().replace('+', 'p').replace('*', 's').replace('(', '_').replace(')', '_').replace(',', '_')


parser = Gaussian94BasisSetParser()
elements = qcel.periodictable.E
os.system("echo '#differing basis sets' > basisdunningdiffer.txt")

with open(sys.argv[1], 'r') as basfile:
    bascontents = basfile.readlines()
bname = bas_sanitize(sys.argv[1])

isdiff = False
if len(sys.argv) > 2:
    isdiff = True
    with open(sys.argv[2], 'r') as reffile:
        refcontents = reffile.readlines()
    rname = bas_sanitize(os.path.basename(sys.argv[2]))

if isdiff:
    if bname != rname:
        print('%s / %s' % (bname, rname), end='')
    else:
        print('%-40s' % (bname), end='')
else:
    print('%-40s' % (bname), end='')

anychange = False
forbiddenchange = False
postKr = False
for el in elements:
    if el.upper() == "RB":
        postKr = True

    shells, msg, ecp_shells, ecp_msg, ecp_ncore = parser.parse(el.upper(), bascontents)

    if isdiff:
        rshells, rmsg, recp_shells, recp_msg, recp_ncore = parser.parse(el.upper(), refcontents)
        if not shells and not rshells:
            print('%s' % ('' if postKr else '   '), end='')
            continue
        if shells and not rshells:
            print(bcolors.OKBLUE + '{:3}'.format(el.upper()) + bcolors.ENDC, end='')
            anychange = True
        if not shells and rshells:
            print(bcolors.FAIL + '{:3}'.format(el.upper()) + bcolors.ENDC, end='')
            anychange = True
            forbiddenchange = True
        if shells and rshells:
            mol = qcdb.Molecule("""\n{}\n""".format(el))
            mol.update_geometry()
            mol.set_basis_all_atoms(bname, role='BASIS')
            bdict = {bname: ''.join(bascontents)}
            rdict = {bname: ''.join(refcontents)}
            bs, msg, ecp = qcdb.BasisSet.construct(parser, mol, 'BASIS', None, bdict, False)
            rbs, rmsg, recp = qcdb.BasisSet.construct(parser, mol, 'BASIS', None, rdict, False)
            #if bs.allclose(rbs, verbose=2):   # see changed coeff/exp
            if bs.allclose(rbs):  # one line per BS
                print('{:3}'.format(el.lower()), end='')
            else:
                print(bcolors.WARNING + '{:3}'.format(el.upper()) + bcolors.ENDC, end='')
                anychange = True
                tbs = bs.print_detail(out='tmpB.txt')
                rtbs = rbs.print_detail(out='tmpR.txt')
                try:
                    outdiff = subprocess.check_output("diff -bwy -W 180 tmpB.txt tmpR.txt >> basisdunningdiffer.txt", shell=True)
                    #outdiff = subprocess.check_output("diff -bw --context=1 tmpB.txt tmpR.txt >> basisdunningdiffer.txt", shell=True)
                except subprocess.CalledProcessError:
                    pass
    else:
        if not shells:
            print('%s' % ('' if postKr else '   '), end='')
        else:
            print('{:3}'.format(el.lower()), end='')
print('')
if anychange and not forbiddenchange:
    os.system("echo 'mv {} ../' >> basisdunningfiles.txt".format(sys.argv[1]))
