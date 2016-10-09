#!/usr/bin/env python
import sys, os

root = os.path.dirname(os.path.realpath(__file__))

# => Driver Code <= #

if __name__ == '__main__':

    # > Working Dirname < #

    if len(sys.argv) == 1:
        dirname = '.'
    elif len(sys.argv) == 2:
        dirname = sys.argv[1]
    else:
        raise Exception('Usage: fsapt.py [dirname]')

    # > Copy Files < #

    os.system('cp %s/pymol2/*pymol %s' % (root, dirname))
