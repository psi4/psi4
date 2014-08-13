#!/usr/bin/env python

# --------------------------------------------------------------------

import sys

sys.path.append('..')

import mpi4py.MPI
import ga

try:
    from signal import signal, SIGPIPE, SIG_IGN
    signal(SIGPIPE, SIG_IGN)
except ImportError:
    pass

# --------------------------------------------------------------------

try:
    from docutils.nodes import NodeVisitor
    NodeVisitor.unknown_visit = lambda self, node: None
    NodeVisitor.unknown_departure =  lambda self, node: None
except ImportError:
    pass

try: # epydoc 3.0.1 + docutils 0.6
    from docutils.nodes import Text
    from UserString import UserString
    if not isinstance(Text, UserString):
        def Text_get_data(s):
            try:
                return s._data
            except AttributeError:
                return s.astext()
        def Text_set_data(s, d):
            s.astext = lambda: d
            s._data = d
        Text.data = property(Text_get_data, Text_set_data)
except ImportError:
    pass

# --------------------------------------------------------------------

from epydoc.docwriter import dotgraph

import re
dotgraph._DOT_VERSION_RE = \
    re.compile(r'dot (?:- Graphviz )version ([\d\.]+)')

try:

    dotgraph.DotGraph.DEFAULT_HTML_IMAGE_FORMAT
    dotgraph.DotGraph.DEFAULT_HTML_IMAGE_FORMAT = 'png'

except AttributeError:

    DotGraph_to_html = dotgraph.DotGraph.to_html
    DotGraph_run_dot = dotgraph.DotGraph._run_dot

    def to_html(self, image_file, image_url, center=True):
        if image_file[-4:] == '.gif':
            image_file = image_file[:-4] + '.png'
        if image_url[-4:] == '.gif':
            image_url = image_url[:-4] +  '.png'
        return DotGraph_to_html(self, image_file, image_url)

    def _run_dot(self, *options):
        if '-Tgif' in options:
            opts = list(options)
            for i, o in enumerate(opts):
                if o == '-Tgif': opts[i] = '-Tpng'
            options = type(options)(opts)
        return DotGraph_run_dot(self, *options)

    dotgraph.DotGraph.to_html = to_html
    dotgraph.DotGraph._run_dot = _run_dot

# --------------------------------------------------------------------

import re

_SIGNATURE_RE = re.compile(
    # Class name (for builtin methods)
    r'^\s*((?P<class>\w+)\.)?' +
    # The function name
    r'(?P<func>\w+)' +
    # The parameters
    r'\(((?P<self>(?:self|cls|mcs)),?)?(?P<params>.*)\)' +
    # The return value (optional)
    r'(\s*(->)\s*(?P<return>\S.*?))?'+
    # The end marker
    r'\s*(\n|\s+(--|<=+>)\s+|$|\.\s+|\.\n)')

from epydoc import docstringparser as dsp
dsp._SIGNATURE_RE = _SIGNATURE_RE

# --------------------------------------------------------------------

import sys, os
import epydoc.cli

def epydocify():
    dirname = os.path.dirname(__file__)
    config = os.path.join(dirname, 'epydoc.cfg')
    sys.argv.append('--config=' + config)
    epydoc.cli.cli()

if __name__ == '__main__':
    epydocify()

# --------------------------------------------------------------------
