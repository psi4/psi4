#!/usr/bin/env python
import os, re, sys

name = sys.argv[1]
path = sys.argv[2]
python = sys.argv[3]

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


fh = open('%s/functional.py' %(python),'r')
lines = fh.readlines()
fh.close()

hook_found = False;
key = name.lower();
re_hook = re.compile(r"^(\s+)'%s'(\s+): build_primitive_functional" %(key))

for line in lines:
    if (re.match(re_hook, line)):
        hook_found = True 

if not hook_found:

    functional_hook = "        '%s'        : build_primitive_functional,\n" %(key)
    superfunctional_hook = "        '%s'   : build_primitive_superfunctional,\n" %(key)

    fh = open('%s/functional.py' %(python),'w')
    
    re_fun2 = re.compile(r"^\s+\S+\s+: build_\S+_functional");
    re_super2 = re.compile(r"^\s+\S+\s+: build_\S+_superfunctional");

    in_fun = False;
    in_super = False;

    for line in lines:

        if (in_fun and (not re.match(re_fun2, line))):
            fh.write(functional_hook)
            in_fun = False            
        if (in_super and (not re.match(re_super2, line))):
            fh.write(superfunctional_hook)
            in_super = False            

        if (re.match(re_fun2, line)):
            in_fun = True;
        if (re.match(re_super2, line)):
            in_super = True;

        fh.write(line)    

    fh.close()
