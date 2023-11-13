#!/usr/bin/env python
import subprocess
import sys
import time

# <<<  run ctest  >>>
retcode = subprocess.Popen(['ctest', '-j2', '-L', "plug|smoke"], bufsize=0,
                            stdout=subprocess.PIPE, universal_newlines=True)
print_all = False
ctestout = ''
while True:
    data = retcode.stdout.readline()
    if not data:
        break

    if '% tests passed,' in data:
        print_all = True

    sdata = data.split()
    test_line = ('Test' in sdata) and ('sec' in sdata)
    start_line = ('Start' in sdata)
    if test_line or start_line or print_all:
        sys.stdout.write(data)  # screen
    #print sys.stdout.write(data)
    ctestout += data  # string

while True:
    retcode.poll()
    exstat = retcode.returncode
    if exstat is not None:
        ctest_exit_status = exstat
        break
    time.sleep(0.1)

# <<<  identify failed tests and cat their output  >>>
sys.stdout.write("""\n  <<<  CTest complete with status %d.  >>>\n\n""" %
                 (ctest_exit_status))

ctestout = str(ctest_exit_status) + "\n" + ctestout

with open('full_ctest_output.dat', 'w') as outfile:
    outfile.write(ctestout)

# if ctest_exit_status:
#     sys.stdout.write("""\n  <<<  CTest failed, printing LastTest.log  >>>\n\n""")
#     with open('Testing/Temporary/LastTest.log', 'r') as ttllog:
#         sys.stdout.write(ttllog.read())
