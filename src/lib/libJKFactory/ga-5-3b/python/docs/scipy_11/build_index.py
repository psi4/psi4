#!/usr/bin/env python

import glob
import os
import sys

sys.path.insert(0, 'publisher')

import options

output_dir = 'output'
dirs = [d for d in glob.glob('%s/*' % output_dir) if os.path.isdir(d)]

pages = []
cum_pages = [1]

for d in sorted(dirs):
    try:
        stats = options.cfg2dict(os.path.join(d, 'paper_stats.cfg'))
        pages.append(int(stats['pages']))
        cum_pages.append(cum_pages[-1] + pages[-1])

        print '"%s" from p. %s to %s' % (os.path.basename(d), cum_pages[-2],
                                         cum_pages[-1] - 1)

        f = open(os.path.join(d, 'page_numbers.tex'), 'w')
        f.write('\setcounter{page}{%s}' % cum_pages[-2])
        f.close()
    except IOError, e:
        continue
