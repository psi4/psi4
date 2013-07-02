#!/usr/bin/python
import sys
if(len(sys.argv) != 2):
    sys.exit("Requires nMAX as an argument");
nMAX = sys.argv[1]
for l in range(0, int(nMAX)):
    for n in range(l+1,int(sys.argv[1])+1):
        print n, l, "0"
