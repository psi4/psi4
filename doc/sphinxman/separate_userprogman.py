#!/usr/bin/python

import sys
import os
import glob


def separate_documents():

    bhtmlall = False
    bhtmluser = False
    bhtmlprog = False
    blatexuser = False
    blatexprog = False
    for line in ftemplate.readlines():
        if line.startswith('.. #####'):
            sline = line.split()
            if (sline[2] == 'HTML-ALL') and (sline[3] == 'START'): bhtmlall = True
            if (sline[2] == 'HTML-ALL') and (sline[3] == 'STOP'): bhtmlall = False
            if (sline[2] == 'HTML-USER') and (sline[3] == 'START'): bhtmluser = True
            if (sline[2] == 'HTML-USER') and (sline[3] == 'STOP'): bhtmluser = False
            if (sline[2] == 'HTML-PROG') and (sline[3] == 'START'): bhtmlprog = True
            if (sline[2] == 'HTML-PROG') and (sline[3] == 'STOP'): bhtmlprog = False
            if (sline[2] == 'LATEX-USER') and (sline[3] == 'START'): blatexuser = True
            if (sline[2] == 'LATEX-USER') and (sline[3] == 'STOP'): blatexuser = False
            if (sline[2] == 'LATEX-PROG') and (sline[3] == 'START'): blatexprog = True
            if (sline[2] == 'LATEX-PROG') and (sline[3] == 'STOP'): blatexprog = False
        else:
            if bhtmlall: fhtmlall.write(line)
            if bhtmluser: fhtmluser.write(line)
            if bhtmlprog: fhtmlprog.write(line)
            if blatexuser: flatexuser.write(line)
            if blatexprog: flatexprog.write(line)
    
    ftemplate.close()
    fhtmlall.close()
    fhtmluser.close()
    fhtmlprog.close()
    flatexuser.close()
    flatexprog.close()


# Partition template_index.rst to form index.rst
ftemplate = open('source/template_index.rst', 'r')
fhtmlall = open('source/autodoc_index_html.rst', 'w')
fhtmluser = open('source/autodoc_index_htmluser.rst', 'w')
fhtmlprog = open('source/autodoc_index_htmlprog.rst', 'w')
flatexuser = open('source/autodoc_index_latexuser.rst', 'w')
flatexprog = open('source/autodoc_index_latexprog.rst', 'w')

separate_documents()

# Partition template_appendices.rst to form appendices.rst
ftemplate = open('source/template_appendices.rst', 'r')
fhtmlall = open('source/autodoc_appendices_html.rst', 'w')
fhtmluser = open('source/autodoc_appendices_htmluser.rst', 'w')
fhtmlprog = open('source/autodoc_appendices_htmlprog.rst', 'w')
flatexuser = open('source/autodoc_appendices_latexuser.rst', 'w')
flatexprog = open('source/autodoc_appendices_latexprog.rst', 'w')

separate_documents()

