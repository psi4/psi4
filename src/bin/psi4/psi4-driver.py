#!/usr/bin/python

# 
# PSI4 driver
# 
# August 2010
#
# The idea here is to parse a simple input file and turn it into a 
# PsiMod python script.  We need to parse the whole thing at once,
# not line-by-line, because we might want to incorporate multi-line
# python code blocks (which might include PSI directives/options inside
# of them which need to be turned into executable PsiMod python code)

import sys;
import re;               # regular expressions

#import psi;
#from psi import *;
#from PsiMod import *;


if (len(sys.argv) < 2):
  ifname = "input.dat"
else:
  ifname = sys.argv[1]

if (len(sys.argv) == 1):
  pfname = "input.py"
elif (len(sys.argv) == 2):
  pfname = ifname + ".py"
else:
  pfname = sys.argv[2]

print "Running with input file ", ifname
print "and python output file ", pfname

infile = open(ifname, 'r')
lines = infile.readlines()
infile.close()
pfile = open(pfname, 'w')

for line in lines:
  matchstr = re.match("(?i)\s*option\s+|\s*!\s+", line);
  if (matchstr) :
    newline = line[0:matchstr.start()]
    tmp = line[matchstr.end():]
    tmp2 = re.split("\s+", tmp)
    key = tmp2[0]
    tmp2.remove(tmp2[0])
    vals = ' '.join(tmp2)
    newline += 'Psimod.set_option("%s", %s);\n' % (key, vals)
  else:
    newline = line  

  pfile.write(newline)

pfile.close()

