#!/usr/bin/env python
import os, sys, re

filename = sys.argv[1]

fh = open(filename, 'r')
lines = fh.readlines()
fh.close()

indices = []
spaces = []
var1 = []
var2 = []
argv = []

re1 = re.compile(r'^(\s+)double (\S+) = exp\((\S+)\);\s*$')
re2 = re.compile(r'^(\s+)double (\S+) = Ei\(\-(\S+)\);\s*$')

for index in range(len(lines)):

    line1 = lines[index]
    if (index < (len(lines) - 1)):
        line2 = lines[index + 1]
    else:
        line2 = ''

    mobj1 = re.match(re1, line1)
    mobj2 = re.match(re2, line2)
    
    if (mobj1 and mobj2 and mobj1.group(3) == mobj2.group(3)):
        indices.append(index)    
        spaces.append(mobj1.group(1))
        var1.append(mobj1.group(2))
        var2.append(mobj2.group(2))
        argv.append(mobj1.group(3))

if (len(indices)):
    
    fh = open(filename, 'w')
    offset = 0;
    instance = 0;
    for index in indices:
        fh.writelines(lines[offset:index])
        fh.write('%sdouble %s;\n' % (spaces[instance], var1[instance])) 
        fh.write('%sdouble %s;\n' % (spaces[instance], var2[instance])) 
        fh.write('%sif (%s > expei_cutoff) {\n' % (spaces[instance], argv[instance]))
        fh.write('%s    %s = 1.0;\n' % (spaces[instance], var1[instance])) 
        fh.write('%s    %s = expei(%s);\n' % (spaces[instance], var2[instance], argv[instance])) 
        fh.write('%s} else {\n' % (spaces[instance]))
        fh.write('%s    %s = exp(%s);\n' % (spaces[instance], var1[instance], argv[instance])) 
        fh.write('%s    %s = Ei(-%s);\n' % (spaces[instance], var2[instance], argv[instance])) 
        fh.write('%s}\n' % ( spaces[instance]))
        offset = index + 2 
        instance = instance + 1;
    fh.writelines(lines[offset:])
    fh.close()
