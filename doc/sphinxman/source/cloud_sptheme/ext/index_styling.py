"""cloud_sptheme.ext.index_styling - improves css styling for genindex"""
import logging; log = logging.getLogger(__name__)
import re
from jinja2 import Markup as literal, escape

prefix = r"^(?P<name>.*)\("
suffix = r"\)$"
_attr_re = re.compile(prefix + r"(?P<left>)(?P<sub>.*)(?P<right> attribute)" + suffix)
_meth_re = re.compile(prefix + r"(?P<left>)(?P<sub>.*)(?P<right> method)" + suffix)
_fc_re = re.compile(prefix + r"(?P<left>class in |in module )(?P<sub>.*)(?P<right>)" + suffix)
_mod_re = re.compile(prefix + r"module" + suffix)

def format_index_name(name):
    while True:
        m = _attr_re.match(name)
        if m:
            name, left, sub, right = m.group("name","left", "sub", "right")
            type = "attribute"
            break
        m = _meth_re.match(name)
        if m:
            name, left, sub, right = m.group("name","left", "sub", "right")
            type = "method"
            break
        m = _fc_re.match(name)
        if m:
            name, left, sub, right = m.group("name","left", "sub", "right")
            if left.startswith("class"):
                type = "class"
            else:
                type = "function"
            break
        m = _mod_re.match(name)
        if m:
            name = m.group("name")
            left = "module"
            sub = right = ''
            type = "module"
            break
        return name
    if sub:
        sub = literal('<span class="subject">') + escape(sub) + literal("</span>")
    cat = left + sub + right
    return escape(name) + literal('<span class="category ' + type + '">') + escape(cat) + literal("</span>")

def mangle_index(app, pagename, templatename, ctx, event_arg):
    if pagename != "genindex":
        return
    fmt = format_index_name
    for key, entries in ctx['genindexentries']:
        for idx, entry in enumerate(entries):
            name, (links, subitems) = entry
            entries[idx] = fmt(name), (links, subitems)
            for idx, entry in enumerate(subitems):
                name, links = entry
                subitems[idx] = fmt(name), links

def setup(app):
    app.connect('html-page-context', mangle_index)
