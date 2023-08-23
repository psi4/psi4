#!/usr/bin/env python
import json
import os
import re
import sys

badtests = []
testfail = re.compile(r'^\s*(?P<num>\d+) - (?P<name>\w+(?:-\w+)*) \(Failed\)\s*$')

with open('full_ctest_output.dat', 'r') as outfile:
    ctestout = outfile.readlines()

ctest_exit_status = int(ctestout[0])
if len(ctestout[1:]) == 0:
    sys.stdout.write("""\n  <<<  All test cases have passed!  >>>\n\n""")
else:
    sys.stdout.write("""\n  <<<  Failing outputs follow.  >>>\n\n""")

for line in ctestout[1:]:
    linematch = testfail.match(line)
    if linematch:
        bad = linematch.group('name')
        sys.stdout.write("""\n\n%s failed. Here is the output:\n""" % (bad))

        badoutfile = bad
        for oddity in [
                'pcmsolver', 'cfour', 'libefp', 'chemps2', 'dftd3', 'mrcc', 'psi4numpy', 'python', 'json', 'cookbook',
                'dkh', 'erd', 'gcp', 'gdma', 'simint', 'snsmp2', 'v2rdm_casscf', 'gpu_dfcc'
        ]:
            if bad.startswith(oddity):
                badoutfile = oddity + '/' + bad
        if "json" in badoutfile:
            badoutfile = 'tests/' + badoutfile + '/output.json'
            if os.path.exists(badoutfile):
                with open(badoutfile, 'r') as ofile:
                    data = json.load(ofile)
                sys.stdout.write(json.dumps(data, indent=2))
                sys.stdout.write("\n")
            else:
                sys.stdout.write("No JSON output file found\n")
        else:
            badoutfile = 'tests/' + badoutfile + '/output.dat'

            if os.path.exists(badoutfile):
                with open(badoutfile, 'r') as ofile:
                    sys.stdout.write(ofile.read())
            else:
                sys.stdout.write("No output file found\n")

# <<<  return ctest error code  >>>
sys.exit(ctest_exit_status)
