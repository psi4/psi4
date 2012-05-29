"""cloud_sptheme.ext.relbar_toc - addes TOC entry to relbar"""
import os.path
import re
import logging; log = logging.getLogger(__name__)

def insert_toc(app, pagename, templatename, ctx, event_arg):
    links = ctx['rellinks']

    #remove any existing toc (present on some pages)
    for idx,  elem in enumerate(links):
        if elem[3] == "toc":
            del links[idx]
            break

    #place toc right after "next" / "previous"
    idx = -1
    for idx, entry in enumerate(links):
        if entry[3] in ("next","previous"):
            break
    else:
        idx += 1

    #insert our toc entry
    #FIXME: there's probably a MUCH better / less broken way to do this
    path = os.path.split(os.path.splitext(ctx['pathto']("contents"))[0])[1]
    ##path = os.path.splitext(ctx['pathto']("contents"))[0]
    ##if path == '':
    ##    path = pagename
    links.insert(idx, (path, "Table Of Contents", "C", "toc"))

def setup(app):
    app.connect('html-page-context', insert_toc)
