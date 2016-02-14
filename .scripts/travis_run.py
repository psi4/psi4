#!/usr/bin/env python
import re
import sys
import time
import subprocess


# <<<  run ctest  >>>
retcode = subprocess.Popen(['ctest', '-j2', '-L', 'quick'], bufsize=0,
                            stdout=subprocess.PIPE, universal_newlines=True)
ctestout = ''
while True:
    data = retcode.stdout.readline()
    if not data:
        break
    sys.stdout.write(data)  # screen
    #tciout.write(data)  # file
    #tciout.flush()
    ctestout += data  # string
while True:
    retcode.poll()
    exstat = retcode.returncode
    if exstat is not None:
        ctest_exit_status = exstat
        break
    time.sleep(0.1)

# <<<  identify failed tests and cat their output  >>>
sys.stdout.write("""\n  <<<  CTest complete with status %d. Failing outputs follow.  >>>\n\n""" %
                 (ctest_exit_status))
badtests = []
testfail = re.compile(r'^\s*(?P<num>\d+) - (?P<name>\w+(?:-\w+)*) \(Failed\)\s*$')

for line in ctestout.split('\n'):
    linematch = testfail.match(line)
    if linematch:
        bad = linematch.group('name')
        sys.stdout.write("""\n%s failed. Here is the output:\n""" % (bad))

        badoutfile = bad
        for oddity in ['pcmsolver', 'cfour', 'libefp', 'dmrg', 'dftd3', 'mrcc']:
            if bad.startswith(oddity):
                badoutfile = oddity + '/' + bad
        badoutfile = 'tests/' + badoutfile + '/output.dat'

        with open(badoutfile, 'r') as ofile:
            sys.stdout.write(ofile.read())

# <<<  return ctest error code  >>>
sys.exit(ctest_exit_status)
