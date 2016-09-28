#! /Users/daniel/anaconda/bin/python
import sys
import psi4

args = sys.argv
if len(args) == 1:
    filename = 'input.dat'
elif len(args) == 2:
    filename = args[-1]
elif len(args) == 5:
    _, filename, outname, _, lib = args
elif len(args) > 2:
    raise Exception("Incomplete number of arguements only psi4 and filename allowed currently.")

with open(filename) as f:
    content = f.read()

processed_content = psi4.process_input(content)
exec(processed_content)


