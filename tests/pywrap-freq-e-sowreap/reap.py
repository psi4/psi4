# /usr/bin/env python
# vim:ft=python

# reap-sow helper script
# executed after the sow step
# AJ
import os
import re
import sys
import subprocess
import time

root_directory = os.path.dirname(os.path.realpath(__file__))


def sowList(outfile):
    """read the output from the 'sow' step to the list of files
    to run before the 'reap' step """
    first_out = os.path.join(root_directory,"output.dat")
    find_cmd=re.compile("^#\s+ psi4\s+-i(?P<infile> (?P<tag>[a-zA-Z]+)-(?P<index>[0-9]+)\.in)")
    the_list=[]
    master = ""
    with open(first_out,'r') as sow_out:
        for line in sow_out:
            #print line
            matchobj=find_cmd.match(line)
            if (matchobj):
                the_list.append(matchobj.group('infile'))
                master = "{0}-master.in".format(matchobj.group('tag'))


    return the_list,master


def storeTests():
    testFile=os.path.join(root_directory,"tests")
    retTests=[]
    with open(testFile) as F:
        for line in F:
            retTests.append(line.strip())

    return retTests


def addTests(theMaster,theTests):
    with open(theMaster,'a') as F:
        for test in theTests:
            F.write("{}\n".format(test))
    pass

def runFiles(psi4Exe,theFileList):

    for infile in theFileList:
        thisFile=os.path.join(root_directory,infile.strip())
        cmd = [psi4Exe,"-i",thisFile]
        subprocess.call(cmd)

def runMaster(psi4Exe,theMasterFile,logfile):
    cmd = [psi4Exe,"-i",theMasterFile]
    try:
        loghandle = open(logfile, 'a')
    except IOError as e:
        print ("""I can't write to %s: %s""" %(logfile, e))
        sys.exit(1)

    try:
        retcode = subprocess.Popen(cmd, bufsize=0, stdout=subprocess.PIPE, universal_newlines=True)
    except OSError as e:
        sys.stderr.write('Command %s execution failed: %s\n' % cmd, e.strerror)
        sys.exit(1)
    p4out = ''
    while True:
        data = retcode.stdout.readline()
        if not data:
            break
        sys.stdout.write(data) # screen
        loghandle.write(data) # file
        loghandle.flush()
        p4out+= data # string

    while True:
        retcode.poll()
        exstat = retcode.returncode
        if exstat is not None:
            return exstat
        time.sleep(0.1)





def main(argv):
    # check the psi4 path is there and is executable
    if os.path.isfile(argv[0]) and os.access(argv[0],os.X_OK):
        # get the list of intermediate input files
        filesList,masterFile= sowList(argv[1])
        # get the commands from tests
        tests=storeTests()
        # set the "reapmode" master file
        reapMaster=os.path.join(root_directory,masterFile)
        # append the tests
        addTests(reapMaster,tests)
        # run intermediates
        runFiles(argv[0],filesList)
        # run the master
        runMaster(argv[0],reapMaster,argv[2])
        # sucess
        sys.exit(0)
    else:
        print("""psi4Exe incorrect check args path %s is not executable""" % (argv[0]))
        sys.exit(1)



if __name__=='__main__':
    print sys.argv[1:]
    main(sys.argv[1:])







