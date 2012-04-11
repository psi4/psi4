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
    import apply_relpath
    IncludePath = apply_relpath.get_topsrcdir_asrelativepathto_objdirsfnxsource()[1]


def pts(category, pyfile):
    print 'Auto-documenting %s file %s' % (category, pyfile)


# Main driver modules in psi4/lib/python
fdriver = open('source/autodoc_driver.rst', 'w')
fdriver.write('\n.. include:: /autodoc_abbr_options_c.rst\n\n')
fdriver.write('.. _`sec:driver`:\n\n')
fdriver.write('=============\n')
fdriver.write('Python Driver\n')
fdriver.write('=============\n\n')

for pyfile in glob.glob(DriverPath + '../../lib/python/*.py'):
    filename = os.path.split(pyfile)[1]
    basename = os.path.splitext(filename)[0]
    div = '=' * len(basename)

    if basename not in ['inpsight', 'pep8', 'sampmod']:

        pts('driver', basename)
        fdriver.write(basename + '\n')
        fdriver.write(div + '\n\n')
    
        fdriver.write('.. automodule:: %s\n' % (basename))
        fdriver.write('   :members:\n')
        fdriver.write('   :undoc-members:\n')
    
        if basename == 'driver':
            fdriver.write('   :exclude-members: energy, optimize, opt, response, frequency, frequencies, freq, property, prop\n')
        elif basename == 'wrappers':
            fdriver.write('   :exclude-members: nbody, cp, counterpoise_correct, counterpoise_correction,\n')
            fdriver.write('       db, database, cbs, complete_basis_set, highest_1, scf_xtpl_helgaker_3,\n')
            fdriver.write('       scf_xtpl_helgaker_2, corl_xtpl_helgaker_2\n')
        elif basename == 'physconst':
            fdriver.write('\n.. literalinclude:: %slib/python/%s\n' % (IncludePath, filename))

    fdriver.write('\n')
fdriver.close()

