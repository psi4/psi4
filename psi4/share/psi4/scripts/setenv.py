#!/usr/bin/env python

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

from __future__ import print_function
import os
import re
import sys
import platform
import subprocess
import collections

# [LAB, 6 Mar 2015]
# This script can be run with *prefix* variable set to install directory
#   (a /bin/psi4 should be present) and used to assess the runtime
#   environment, particularly the following. Designed to be run after
#   download of binary distribution where conda will adjust *prefix*
#   variable properly.
#   * set PSI_SCRATCH
#   * set PSIPATH
#   * set psi4 in PATH
#   * check linked libraries all found
#   * for binary distributions
#     * check libc compatibility
#     * check right libpython linked

prefix = '/opt/anaconda1anaconda2anaconda3'
envarray = collections.OrderedDict()
# PYTHONPATH
# PYTHONHOME


def which(program, envvar='PATH', startswith=False):
    """Thanks, http://stackoverflow.com/a/377028"""
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ[envvar].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

if platform.system() == 'Linux':

    # <<<  inspect libc  >>>
    libcdist, builtlibc = platform.libc_ver(executable=prefix + '/bin/psi4')
    builtver = tuple(builtlibc.split('.'))
    ans = subprocess.Popen('ldd --version', shell=True, stdout=subprocess.PIPE).stdout.read()
    sysver = tuple(ans.split('\n')[0].split()[-1].split('.'))
    print('\n  System libc is version %s, and conda psi4 requires libc version %s or higher' %
        ('.'.join(sysver), '.'.join(builtver)), end='')
    if (int(sysver[0]) < int(builtver[0])) or (int(sysver[1]) < int(builtver[1])):
        print('\n  ERROR: incompatible libc (%s < %s), no simple remedy.' % ('.'.join(sysver), '.'.join(builtver)))
        sys.exit()
    else:
        print(' .... GOOD')

    # <<<  inspect psi4 library linkage  >>>
    quitafterstep = False
    ans = subprocess.Popen('ldd %s' % (prefix + '/bin/psi4'), shell=True, stdout=subprocess.PIPE).stdout.read()
    libraries = {}
    for line in ans.splitlines():
        match = re.match(r'\t(.*) => (.*) \(0x', line)
        ack = re.match(r'\t(.*) => not found', line)
        if match:
            libraries[match.group(1)] = match.group(2)
        if ack:
            libraries[ack.group(1)] = None
    print('\n  Library dependencies of conda psi4 in current environment:\n%s' % (ans))

    for name, fl in libraries.items():
        #print('%30s    %s' % (name, fl))
        if fl is None:
            print('  ERROR: needed library %s not found. no advice at this time.' % (name))
            quitafterstep = True
        if name.startswith('libpython'):
            if 'conda' in fl:
                print('  {} in conda distribution .... GOOD'.format(fl))
            else:
                print('  {} in system distribution .... FIXABLE'.format(fl))

    if quitafterstep:
        print('\n  Your conda psi4 installation not expected to work.')
        sys.exit()

    # <<<  inspect ld_library_path  >>>
    print("""
  Psi4 uses Python to drive C++ code. This means that the executable was
    linked against a python library at conda-build time, and psi4 finds a
    python library and uses the python interpreter at runtime. For compatibility,
    the executable must find the *conda* libpython at runtime rather than
    any system libpython. Creating a conda environment containing python and
    psi4 sets up the optimal conditions for library resolving. However, the
    environment variable LD_LIBRARY_PATH has precedence over the conda
    environment settings so that a libpython in a LD_LIBRARY_PATH directory
    will likely cause psi4 to fail. Below scans your LD_LIBRARY_PATH for any
    potentially meddlesome directories and suggests a safe truncated
    LD_LIBRARY_PATH. Only you understand the delicate balance of entries in
    this variable and runtimes requirements of all your programs, so you're
    free to seek another solution.""")
    temp = []
    ans = os.getenv('LD_LIBRARY_PATH', None)
    if ans is None:
        print('  Environment variable LD_LIBRARY_PATH is {} .... GOOD'.format(ans))
        envarray['LD_LIBRARY_PATH'] = ''
    else:
        print('  Environment variable LD_LIBRARY_PATH contains:')
        ans = [os.path.abspath(os.path.expandvars(os.path.expanduser(p))) for p in ans.split(os.pathsep)]
        for dpath in ans:
            print('    {}'.format(dpath), end='')
            try:
                meddlesome = [fl for fl in os.listdir(dpath) if os.path.isfile(os.path.join(dpath, fl)) and fl.startswith('libpython')]
            except OSError:
                print(' .... FIXABLE')
                print('      directory does not exist so discarding.')
            else:
                if meddlesome:
                    print(' .... FIXABLE')
                    print('      directory contains interfering {} so discarding.'.format(meddlesome[0]))
                else:
                    print(' .... GOOD')
                    temp.append(dpath)

        print('  Environment variable LD_LIBRARY_PATH contains:')
        envarray['LD_LIBRARY_PATH'] = ':'.join(temp)
        print('    {} .... GOOD'.format(envarray['LD_LIBRARY_PATH']))

    # <<<  inspect psi_scratch  >>>
    print("""
  Quantum chemistry software writes many temporary files of large size.
    The environment variable PSI_SCRATCH controls where these files go.
    You must have write permissions to this location. If at *all* possible,
    this location should be a local, not NFS-mounted, disk.""")
    user_obedient = False
    while not user_obedient:
        ans = os.getenv('PSI_SCRATCH', None)
        print('  Environment variable PSI_SCRATCH is %s' % (ans), end='')
        if ans is not None and os.access(ans, os.W_OK | os.X_OK):
            print(' .... GOOD')
            user_obedient = True
        else:
            print(' .... FIXABLE')
            if ans is not None and not os.access(ans, os.W_OK | os.X_OK):
                print('  Path %s is not existant/writeable .... FIXABLE' % (ans))
                print('  NOTE: if this is just a fixable permissions problem: cancel script, adjust permissions, re-run script.')
            temp = raw_input('  TASK: specify a scratch directory.\n    PSI_SCRATCH [/tmp] =  ').strip()
            if temp == '':
                temp = '/tmp'
            trial = os.path.abspath(os.path.expandvars(os.path.expanduser(temp)))
            os.environ['PSI_SCRATCH'] = trial
    envarray['PSI_SCRATCH'] = ans

    # <<<  inspect psipath and path  >>>
    print("""
  Psi4 can use certain files outside its library for basis sets, efp
    fragments, interfaced executables, etc. The environment variable
    PSIPATH is a colon-separated list of directories to search for
    these files, akin to PATH for executables.""")
    ans = os.getenv('PSIPATH', None)
    envarray['PSIPATH'] = '' if ans is None else ans
    print('  Environment variable PSIPATH is {} .... GOOD'.format(ans))
    for exe in ['dftd3', 'dmrcc', 'xcfour']:
        if ans is None:
            exeinpath = which(exe)
            if exeinpath is not None:
                temp = os.path.dirname(exeinpath)
                print('  Executable {} found in PATH so adding {} to PSIPATH .... GOOD'.format(exe, temp))
                envarray['PSIPATH'] += ':' + temp
        else:
            exeinpsipath = which(exe, envvar='PSIPATH')
            if exeinpsipath is None:
                exeinpath = which(exe)
                if exeinpath is not None:
                    temp = os.path.dirname(exeinpath)
                    print('  Executable {} found in PATH so adding {} to PSIPATH .... GOOD'.format(exe, temp))
                    envarray['PSIPATH'] += ':' + temp
            else:
                print('  Executable {} found in PSIPATH .... GOOD'.format(exe))
    temp = raw_input('  TASK: specify additional paths for auxiliary items (colon separated).\n    PSIPATH [empty] +=  ').strip()
    envarray['PSIPATH'] += ':' + temp
    print('  Environment variable PSIPATH is %s .... GOOD' % (envarray['PSIPATH']))

    # <<<  inspect path to psi4  >>>
    print()
    psi4inpath = which('psi4')
    if psi4inpath is None:
        print('  Psi4 executable not in PATH .... FIXABLE')
        envarray['PATH'] = prefix + '/bin:' + os.getenv('PATH', '')
    else:
        temp = os.path.dirname(psi4inpath)
        if temp != (prefix + '/bin'):
            print('  Wrong psi4 executable in PATH .... FIXABLE')
            envarray['PATH'] = prefix + '/bin:' + os.getenv('PATH', '')
        else:
            print('  Psi4 executable in PATH .... GOOD')
    envarray['PSIDATADIR'] = prefix + '/share/psi4'

    # <<<  write out environment files  >>>
    print("""
  The runtime environment changes from this script are written to files
    "psi4setup.sh" and "psi4setup.csh" for bash- and csh-type shells,
    respectively. These files are echoed below.""")
    setup = prefix + '/share/psi/scripts/psi4setup.sh'
    with open(setup, 'w') as handle:
        for ev, val in envarray.items():
            handle.write("""export {}={}\n""".format(ev.upper(), val))
    print("""\n  bash-prompt>>> source {}""".format(setup))
    os.system('cat {}'.format(setup, 'sh'))
    setup = prefix + '/share/psi/scripts/psi4setup.csh'
    with open(setup, 'w') as handle:
        for ev, val in envarray.items():
            handle.write("""setenv {} {}\n""".format(ev.upper(), val))
    print("""\n  csh-prompt>>> source {}""".format(setup))
    os.system('cat {}'.format(setup))
    with open('setenvtest.in', 'w') as handle:
        handle.write("""
# Any line starting with the # character is a comment line
#! Sample HF/cc-pVDZ H2O computation

memory 250 mb

molecule h2o {
  O
  H 1 0.96
  H 1 0.96 2 104.5
}

set basis cc-pVDZ
energy('scf')

compare_values(-76.0266327341067125, get_variable('SCF TOTAL ENERGY'), 6, 'SCF energy')  #TEST

""")

    # <<<  test the environment  >>>
    lenv = os.environ.copy()
    for ev, val in envarray.items():
        lenv[ev.upper()] = val
    print("""
  Your shell environment has not been changed by this script. A psi4 test
    job will be run under the environment conditions suggested by this script.""")
    ans = subprocess.Popen('psi4 setenvtest.in'.format(prefix), shell=True, stdout=subprocess.PIPE, env=lenv).stdout.read()
    print(ans)
