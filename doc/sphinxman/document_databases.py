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


# Available databases in psi4/lib/databases
fdriver = open('source/autodoc_available_databases.rst', 'w')
fdriver.write('\n\n')

for pyfile in glob.glob(DriverPath + '../../lib/databases/*.py'):
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

