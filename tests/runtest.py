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
import sys
import time
import subprocess

if len(sys.argv) not in [5, 6, 7, 8, 9, 10]:
    print("""Usage: %s input_file logfile doperltest top_srcdir doreaptest alt_output_file alt_psi4_exe alt_psi4datadir""" % (sys.argv[0]))
    sys.exit(1)

# extract run condition from arguments
python_exec = sys.argv[0]
infile = sys.argv[1]
logfile = sys.argv[2]
psiautotest = sys.argv[3]
top_srcdir = sys.argv[4]
sowreap = sys.argv[5]

if len(sys.argv) >= 7:
    outfile = sys.argv[6]
else:
    outfile = 'output.dat'

if len(sys.argv) >= 8:
    psi = sys.argv[7]
else:
    psi = '../../bin/psi4'

if len(sys.argv) >= 9:
    psidatadir = sys.argv[8]
else:
    psidatadir = os.path.dirname(os.path.realpath(psi)) + '/../share/psi4'

if len(sys.argv) >= 10:
    psilibdir = sys.argv[9] + os.path.sep
else:
    psilibdir = os.path.abspath('/../')

# open logfile and print test case header
try:
    loghandle = open(logfile, 'a')
except IOError as e:
    print("""I can't write to %s: %s""" % (logfile, e))
loghandle.write("""\n%s\n%s\n""" % (os.path.dirname(infile).split(os.sep)[-1], time.strftime("%Y-%m-%d %H:%M")))


def backtick(exelist):
    """Executes the command-argument list in *exelist*, directing the
    standard output to screen and file logfile and string p4out. Returns
    the system status of the call.

    """
    try:
        retcode = subprocess.Popen(exelist, bufsize=0, stdout=subprocess.PIPE, universal_newlines=True)
    except OSError as e:
        sys.stderr.write('Command %s execution failed: %s\n' % (exelist, e.strerror))
        sys.exit(1)

    p4out = ''
    while True:
        data = retcode.stdout.readline()
        if not data:
            break
        sys.stdout.write(data)  # screen
        loghandle.write(data)  # file
        loghandle.flush()
        p4out += data  # string
    while True:
        retcode.poll()
        exstat = retcode.returncode
        if exstat is not None:
            return exstat
        time.sleep(0.1)
    loghandle.close()
    # not sure why 2nd while loop needed, as 1st while loop has always
    #   been adequate for driver interfaces. nevertheless, to collect
    #   the proper exit code, 2nd while loop very necessary.

# run psi4 and collect testing status from any compare_* in input file
if os.path.isfile(infile):
    pyexitcode = backtick([psi, infile, outfile, '-l', psidatadir])
elif os.path.isfile(infile.replace(".dat", ".py")):
    infile = infile.replace(".dat", ".py")
    os.environ["PYTHONPATH"] = psilibdir
    outfile = os.path.dirname(infile) + os.path.sep + outfile
    pyexitcode = backtick(["python", infile, " > ", outfile])
else:
    raise Exception("\n\nError: Input file %s not found\n" % infile)

if sowreap == 'true':
    try:
        retcode = subprocess.Popen([sys.executable, '%s/tests/reap.py' %
                  (top_srcdir), infile, outfile, logfile, psi, psidatadir])
    except OSError as e:
        print("""Can't find reap script: %s """ % (e))
    while True:
        retcode.poll()
        exstat = retcode.returncode
        if exstat is not None:
            reapexitcode = exstat
            break
        time.sleep(0.1)
else:
    reapexitcode = None

# additionally invoke autotest script comparing output.dat to output.ref
if psiautotest == 'true':
    os.environ['SRCDIR'] = os.path.dirname(infile)
    try:
        retcode = subprocess.Popen(['perl', '%s/tests/psitest.pl' % (top_srcdir), infile, logfile])
    except IOError as e:
        print("""Can't find psitest script: %s""" % (e))
    while True:
        retcode.poll()
        exstat = retcode.returncode
        if exstat is not None:
            plexitcode = exstat
            break
        time.sleep(0.1)
else:
    plexitcode = None

# combine, print, and return (0/1) testing status
exitcode = 0 if (pyexitcode == 0 and (plexitcode is None or plexitcode == 0) and (reapexitcode is None or reapexitcode == 0)) else 1
print('Exit Status: infile (', pyexitcode, '); autotest (', plexitcode, '); sowreap (', reapexitcode, '); overall (', exitcode, ')')
sys.exit(exitcode)
