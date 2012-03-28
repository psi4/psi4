#!/usr/bin/env python
import os, re, sys

name = sys.argv[1]
path = sys.argv[2]

os.system('mv %sfunctional.* %s/' %(name, path))

fh = open('%s/factory.cc' %(path), 'r')
lines = fh.readlines();
fh.close()

header_found = False;
re_header = re.compile(r'^#include "%sfunctional\.h"' %(name))

for line in lines:
    if (re.match(re_header, line)):
        header_found = True 

if not header_found:

    re_inc = re.compile(r'^#include')
    inc = True;
    re_el = re.compile(r'^\s+\} else \{')
    
    fh = open('%s/factory.cc' %(path), 'w')
    
    for line in lines:
        if inc:
            if not re.match(re_inc, line):
                fh.write('#include "%sfunctional.h"\n' %(name))
                inc = False
    
        if re.match(re_el, line):
            fh.write('    } else if (alias == "%s") {\n' %(name))
            fh.write('        fun = new %sFunctional();\n' %(name))
            el = False
    
        fh.write(line)    
    
    fh.close()

