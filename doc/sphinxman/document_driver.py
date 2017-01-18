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
InsertPath = '/../../../'
if (len(sys.argv) == 2):
    DriverPath = sys.argv[1] + '/'
    sys.path.insert(0, os.path.abspath(os.getcwd()))


def pts(category, pyfile):
    print('Auto-documenting %s file %s' % (category, pyfile))


# Main driver modules in psi4/driver
fdriver = open('source/autodoc_driver.rst', 'w')
fdriver.write('\n.. include:: /autodoc_abbr_options_c.rst\n\n')
fdriver.write('.. _`sec:driver`:\n\n')
fdriver.write('=============\n')
fdriver.write('Python Driver\n')
fdriver.write('=============\n\n')

for pyfile in glob.glob(DriverPath + '../../psi4/driver/*.py'):
    filename = os.path.split(pyfile)[1]
    basename = os.path.splitext(filename)[0]
    div = '=' * len(basename)

    if basename not in ['inpsight', 'pep8', 'diatomic_fits', 'pyparsing', 'computation_cache']:

        pts('driver', basename)
        fdriver.write(basename + '\n')
        fdriver.write(div + '\n\n')

        fdriver.write('.. automodule:: %s\n' % (basename))
        fdriver.write('   :members:\n')
        fdriver.write('   :undoc-members:\n')

        if basename == 'driver':
            fdriver.write('   :exclude-members: energy, optimize, opt, frequency, frequencies, freq, property, prop, molden, gdma, fchk, gradient, hessian\n')
        elif basename == 'wrapper_database':
            fdriver.write('   :exclude-members: db, database\n')
        elif basename == 'driver_nbody':
            fdriver.write('   :exclude-members: nbody_gufunc\n')
        elif basename == 'driver_cbs':
            fdriver.write('   :exclude-members: cbs, complete_basis_set, xtpl_highest_1,\n')
            fdriver.write('       scf_xtpl_helgaker_3, scf_xtpl_helgaker_2, corl_xtpl_helgaker_2, n_body\n')
#        elif basename == 'physconst':
#            fdriver.write('\n.. literalinclude:: %sdriver/%s\n' % (IncludePath, filename))
        elif basename == 'diatomic':
            fdriver.write('   :exclude-members: anharmonicity\n')
#        elif basename == 'interface_dftd3':
#            fdriver.write('   :exclude-members: run_dftd3\n')
#        elif basename == 'interface_cfour':
#            fdriver.write('   :exclude-members: run_cfour\n')
        elif basename == 'aliases':
            fdriver.write('   :exclude-members: sherrill_gold_standard, allen_focal_point\n')
        elif basename == 'p4util':
            fdriver.write('   :exclude-members: oeprop, cubeprop\n')
        elif basename == 'procedures':
            fdriver.write('   :exclude-members: interface_cfour\n')

    fdriver.write('\n')


# Python-only plugin modules in psi4/driver
for basename in os.walk(DriverPath + '../../psi4/driver').next()[1]:
    div = '=' * len(basename)

    if basename not in ['grendel']:

        pts('driver', basename)
        fdriver.write(basename + '\n')
        fdriver.write(div + '\n\n')

        fdriver.write('.. automodule:: %s\n' % (basename))
        fdriver.write('   :members:\n')
        fdriver.write('   :undoc-members:\n')

        for pyfile in glob.glob(DriverPath + '../../psi4/driver/' + basename + '/*py'):
            filename = os.path.split(pyfile)[1]
            basename2 = os.path.splitext(filename)[0]
            div = '=' * len(basename2)

            fdriver.write('.. automodule:: %s.%s\n' % (basename, basename2))
            fdriver.write('   :members:\n')
            fdriver.write('   :undoc-members:\n')

        if basename == 'qcdb' and basename2 == 'interface_dftd3':
            fdriver.write('   :exclude-members: run_dftd3\n')

        fdriver.write('\n')
    fdriver.write('\n')
fdriver.close()
