import os
import sys
import time
import subprocess

if len(sys.argv) not in [5, 6, 7]:
    print """Usage: %s input_file logfile doperltest top_srcdir [alt_output_file alt_psi4_exe]""" % (sys.argv[0])
    sys.exit(1)

# extract run condition from arguments
infile = sys.argv[1]
logfile = sys.argv[2]
psiautotest = sys.argv[3]
top_srcdir = sys.argv[4]
psidatadir = top_srcdir + '/lib'
if len(sys.argv) in [6, 7]:
    outfile = sys.argv[5]
else:
    outfile = 'output.dat'
if len(sys.argv) == 7:
    psi = sys.argv[6]
else:
    psi = '../../bin/psi4'

# open logfile and print test case header
try:
    loghandle = open(logfile, 'a')
except IOError as e:
    print """I can't write to %s: %s""" % (logfile, e)
loghandle.write("""\n%s\n%s\n""" % (os.path.dirname(infile).split(os.sep)[-1], time.strftime("%Y-%m-%d %H:%M")))


def backtick(exelist):
    """Executes the command-argument list in *exelist*, directing the
    standard output to screen and file logfile and string p4out. Returns
    the system status of the call.

    """
    try:
        retcode = subprocess.Popen(exelist, bufsize=0, stdout=subprocess.PIPE)
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
    # not sure why 2nd while loop needed, as 1st while loop has always
    #   been adequate for driver interfaces. nevertheless, to collect
    #   the proper exit code, 2nd while loop very necessary.

# run psi4 and collect testing status from any compare_* in input file
pyexitcode = backtick([psi, infile, outfile, '-l', psidatadir])
loghandle.close()

# additionally invoke autotest script comparing output.dat to output.ref
if psiautotest == 'true':
    os.environ['SRCDIR'] = os.path.dirname(infile)
    try:
        retcode = subprocess.Popen(['perl', '%s/tests/psitest.pl' % (top_srcdir), infile, logfile])
    except IOError as e:
        print """Can't find psitest script: %s""" % (e)
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
exitcode = 0 if (pyexitcode == 0 and (plexitcode is None or plexitcode == 0)) else 1
print 'Exit Status: infile (', pyexitcode, '); autotest (', plexitcode, '); overall (', exitcode, ')'
sys.exit(exitcode)
