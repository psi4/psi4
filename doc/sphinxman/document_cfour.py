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
    print 'Auto-documenting %s module %s' % (category, pyfile)


# Available psi variables in psi4/lib/python/qcdb/cfour.py
fdriver = open('source/autodir_psivariables/module__cfour.rst', 'w')
fdriver.write('\n\n')

psivars = []
for pyfile in glob.glob(DriverPath + '../../lib/python/qcdb/cfour.py'):
    filename = os.path.split(pyfile)[1]
    basename = os.path.splitext(filename)[0]
    div = '=' * len(basename)
    print pyfile, filename, basename

    if basename not in []:
        pts('psi variables', basename)
    
        fdriver.write('.. _`apdx:%s_psivar`:\n\n' % (basename.lower()))
        fdriver.write('\n%s\n%s\n\n' % (basename.upper(), '"' * len(basename)))
        fdriver.write('.. hlist::\n   :columns: 1\n\n')

        f = open(pyfile)
        contents = f.readlines()
        f.close()

        for line in contents:
            mobj = re.search(r"""^\s*psivar\[\'(.*)\'\]\s*=""", line)
            if mobj:
                if mobj.group(1) not in psivars:
                    psivars.append(mobj.group(1))

for pv in sorted(psivars):
    pvsquashed = pv.replace(' ', '')
    fdriver.write('   * :psivar:`%s <%s>`\n\n' % (pv, pvsquashed))

fdriver.write('\n')
fdriver.close()

fdriver = open('source/autodoc_psivariables_bymodule.rst', 'a')
fdriver.write('   autodir_psivariables/module__cfour\n\n')
fdriver.close()

