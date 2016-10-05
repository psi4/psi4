#!/usr/bin/env python
import re
import sys
import time
import subprocess


badtests = []
testfail = re.compile(r'^\s*(?P<num>\d+) - (?P<name>\w+(?:-\w+)*) \(Failed\)\s*$')

with open('full_ctest_output.dat', 'r') as outfile:
    ctestout = outfile.readlines()

ctest_exit_status = ctestout[0]

for line in ctestout[1:]:
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
