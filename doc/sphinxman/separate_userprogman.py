#!/usr/bin/python

import sys
import os
import glob

findex = open('source/template_index.rst', 'r')

fhtml = open('source/index_html.rst', 'w')
fuser = open('source/index_latexuser.rst', 'w')
fprog = open('source/index_latexprog.rst', 'w')

bhtml = False
buser = False
bprog = False
for line in findex.readlines():
    if line.startswith('.. #####'):
        sline = line.split()
        if (sline[2] == 'HTML') and (sline[3] == 'START'): bhtml = True
        if (sline[2] == 'HTML') and (sline[3] == 'STOP'): bhtml = False
        if (sline[2] == 'LATEX-USER') and (sline[3] == 'START'): buser = True
        if (sline[2] == 'LATEX-USER') and (sline[3] == 'STOP'): buser = False
        if (sline[2] == 'LATEX-PROG') and (sline[3] == 'START'): bprog = True
        if (sline[2] == 'LATEX-PROG') and (sline[3] == 'STOP'): bprog = False
    else:
        if bhtml:
            fhtml.write(line)
        if buser:
            fuser.write(line)
        if bprog:
            fprog.write(line)

findex.close()
fhtml.close()
fuser.close()
fprog.close()

