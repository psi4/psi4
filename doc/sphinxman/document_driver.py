#!/usr/bin/python

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
    print 'Auto-documenting %s file %s' % (category, pyfile)


# Main driver modules in psi4/share/python
fdriver = open('source/autodoc_driver.rst', 'w')
fdriver.write('\n.. include:: /autodoc_abbr_options_c.rst\n\n')
fdriver.write('.. _`sec:driver`:\n\n')
fdriver.write('=============\n')
fdriver.write('Python Driver\n')
fdriver.write('=============\n\n')

for pyfile in glob.glob(DriverPath + '../../share/python/*.py'):
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
            fdriver.write('   :exclude-members: _nbody_gufunc\n')
        elif basename == 'driver_cbs':
            fdriver.write('   :exclude-members: cbs, complete_basis_set, xtpl_highest_1,\n')
            fdriver.write('       scf_xtpl_helgaker_3, scf_xtpl_helgaker_2, corl_xtpl_helgaker_2, n_body\n')
#        elif basename == 'physconst':
#            fdriver.write('\n.. literalinclude:: %sshare/python/%s\n' % (IncludePath, filename))
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


# Python-only plugin modules in psi4/share/python
for basename in os.walk(DriverPath + '../../share/python').next()[1]:
    div = '=' * len(basename)

    if basename not in ['grendel']:

        pts('driver', basename)
        fdriver.write(basename + '\n')
        fdriver.write(div + '\n\n')

        fdriver.write('.. automodule:: %s\n' % (basename))
        fdriver.write('   :members:\n')
        fdriver.write('   :undoc-members:\n')

        for pyfile in glob.glob(DriverPath + '../../share/python/' + basename + '/*py'):
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
