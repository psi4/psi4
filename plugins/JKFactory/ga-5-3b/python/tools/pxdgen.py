#!/usr/bin/env python
"""

Generate the gah.pxd source from the ga.h header.
We basically remove the "extern" and the ";" since our ga.h header is
one-function-per-line.  We then add other cython things to the top.

Usage:
    pxdgen.py path/to/ga.h > gah.pxd

"""
import sys

if len(sys.argv) != 2:
    print 'incorrect number of arguments'
    print 'usage: pxdgen.py <ga.h> > <gah.pxd>'
    sys.exit(len(sys.argv))

# print headers
print '''from libc.stdio  cimport FILE
from libc.stdint cimport int64_t

cdef extern from "typesf2c.h":
    ctypedef int Integer
    ctypedef float Real
    ctypedef double DoublePrecision
    ctypedef struct DoubleComplex:
        DoublePrecision real
        DoublePrecision imag
    ctypedef struct SingleComplex:
        Real real
        Real imag

cdef extern from "ga.h":
    ctypedef Integer ga_nbhdl_t
'''
for line in open(sys.argv[1]):
    line = line.strip()
    if 'extern' in line and '{' in line:
        continue
    elif 'extern' in line:
        line = line.replace('extern','')
        line = line.replace(';','')
        if '(void)' in line:
            line = line.replace('(void)', '()')
        line = line.strip()
        print "    %s" % line
