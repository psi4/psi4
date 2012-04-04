#!/usr/bin/env python
import os, re, sys

name = sys.argv[1]

os.system("sed -i 's/NAME/%s/g' %sfunctional.cc" %(name,name))
os.system("sed -i 's/NAME/%s/g' %sfunctional.h" %(name,name))

fh = open('preamble', 'r')
preamble = fh.readlines()
fh.close()

fh = open('parameters', 'r')
parameters = fh.readlines()
fh.close()

fh = open('functional', 'r')
functional = fh.readlines()
fh.close()

fh = open('functional_rho_a0', 'r')
functional_rho_a0 = fh.readlines()
fh.close()

fh = open('functional_rho_b0', 'r')
functional_rho_b0 = fh.readlines()
fh.close()

fh = open('%sfunctional.cc' % (name), 'r')
template = fh.readlines()
fh.close()

fh = open('%sfunctional.cc' % (name), 'w')

pre_re    = re.compile(r'^(\s+)PREAMBLE$')
par_re    = re.compile(r'^(\s+)PARAMETERS$')
fun_re    = re.compile(r'^(\s+)FUNCTIONAL$')
fun_a0_re = re.compile(r'^(\s+)FUNCTIONAL_RHO_A0$')
fun_b0_re = re.compile(r'^(\s+)FUNCTIONAL_RHO_B0$')

for line in template:
    mobj = re.match(pre_re, line);
    if mobj:
        spaces = mobj.group(1)
        for l in preamble:
            fh.write(spaces + l)
        continue
    mobj = re.match(par_re, line);
    if mobj:
        spaces = mobj.group(1)
        for l in parameters:
            fh.write(spaces + l)
        continue
    mobj = re.match(fun_re, line);
    if mobj:
        spaces = mobj.group(1)
        for l in functional:
            fh.write(spaces + l)
        continue
    mobj = re.match(fun_a0_re, line);
    if mobj:
        spaces = mobj.group(1)
        for l in functional_rho_a0:
            fh.write(spaces + l)
        continue
    mobj = re.match(fun_b0_re, line);
    if mobj:
        spaces = mobj.group(1)
        for l in functional_rho_b0:
            fh.write(spaces + l)
        continue
    fh.write(line)

os.system('rm preamble')
os.system('rm parameters')
os.system('rm functional')
os.system('rm functional_rho_a0')
os.system('rm functional_rho_b0')

fh.close()
