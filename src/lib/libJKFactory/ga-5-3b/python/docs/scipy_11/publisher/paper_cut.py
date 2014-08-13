"""
Parse pdflatex output for paper count, and store in a .ini file.
"""

import sys

import re

regexp = re.compile('Output written on paper.pdf \((\d+) pages')

f = open('paper_stats.cfg', 'w')

line = ''
for line in sys.stdin:
    m = regexp.match(line)
    if m:
        pages = m.groups()[0]
        f.writelines(['[default]\n', 'pages = %s\n' % pages])

f.close()

