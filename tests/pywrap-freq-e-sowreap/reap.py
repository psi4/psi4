# /usr/bin/env python
# vim:ft=python

# reap-sow helper script
# executed after the sow step
# AJ
import os
import re
import sys
import subprocess

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

def runMaster(psi4Exe,theMasterFile):
    cmd = [psi4Exe,"-i",theMasterFile]
    subprocess.call(cmd)


def main(argv):
    if os.path.isfile(argv[0]):
        filesList,masterFile= sowList(argv[1])
        tests=storeTests()
        reapMaster=os.path.join(root_directory,masterFile)
        addTests(reapMaster,tests)

        runFiles(argv[0],filesList)
        runMaster(argv[0],reapMaster)

        sys.exit(0)
    else:
        sys.exit(1)



if __name__=='__main__':
    print sys.argv[1:]
    main(sys.argv[1:])







