"""psi4doc_sptheme.ext.relbar_toc - adds TOC entry and icons to relbar.
Modeled from cloud_sptheme.ext.relbar_toc

"""
import os.path
import logging; log = logging.getLogger(__name__)

def insert_toc(app, pagename, templatename, ctx, event_arg):
    links = ctx['rellinks']

    # remove any existing toc (present on some pages)
    for idx, entry in enumerate(links):
        #if entry[3].lower() == "toc":
        if entry[2] == 'C':
            del links[idx]
            break

    # add arrows and caps to existing links
    for idx, entry in enumerate(links):
        if entry[3] == 'next':
            #newlink = u'next \u2192'
            #newlink = u'<i class="fa fa-arrow-circle-right fa-lg"></i>'
            newlink = u'<i class="fa fa-long-arrow-right fa-lg"></i>'
        elif entry[3] == 'previous':
            #newlink = u'\u2190 previous'
            #newlink = u'<i class="fa fa-arrow-circle-left fa-lg"></i>'
            newlink = u'<i class="fa fa-long-arrow-left fa-lg"></i>'
        elif entry[3] == 'index':
            newlink = 'Index'
        else:
            continue
        del links[idx]
        links.insert(idx, (entry[0], entry[1], entry[2], newlink))

    # insert our toc entry
    #FIXME: there's probably a MUCH better / less broken way to do this
    path = os.path.split(os.path.splitext(ctx['pathto']('index'))[0])[1]
    #links.insert(len(links), (path, 'Table Of Contents', 'C', 'TOC'))
    links.insert(len(links), (path, 'Table Of Contents', 'C', '<i class="fa fa-book fa-lg"></i>'))


def setup(app):
    app.connect('html-page-context', insert_toc)
