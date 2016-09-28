"""cloud_sptheme.ext.perpage -- override sphinx config per-page

* perpage_html_logo -- glob map (ala html_sidebars), used to change sidebar logo per-page

.. todo::
    document this extension
"""
#=============================================================================
# imports
#=============================================================================
# core
import os.path
import re
import logging; log = logging.getLogger(__name__)
# site
from sphinx.util.matching import patmatch
# local
#=============================================================================
# helpers
#=============================================================================
def _rank_pattern(pattern):
    """return sorting key for prioritizing which glob pattern should match"""
    # TODO: add more ways to distinguish patterns if both have wildcards
    return not any(char in pattern for char in '*?[')

def bestmatch(patmap, source, default=None, param="source"):
    """return best match given a dictionary mapping glob pattersn -> values"""
    best = None
    best_rank = None
    for pattern in patmap:
        if not patmatch(pattern, source):
            continue
        cur_rank = _rank_pattern(pattern)
        if best is None or cur_rank < best_rank:
            best = pattern
            best_rank = cur_rank
        elif cur_rank == best_rank:
            raise KeyError("%s %r matches too many patterns: %r and %r" %
                           (param, source, best, pattern))
    if best is None:
        return default
    else:
        return patmap[best]

#=============================================================================
# sphinx hooks
#=============================================================================
def perpage_html_logo(app, pagename, templatename, ctx, event_arg):
    """helper to override sidebar logo per-page"""
    patmap = getattr(app.config, "perpage_html_logo", {})
    logo = bestmatch(patmap, pagename, ctx.get("logo"), param="pagename")
    if logo is None:
        ctx.pop("logo", None)
    else:
        ctx['logo'] = logo

# NOTE: this works, just don't have a use for it yet
##def perpage_html_theme_options(app, pagename, templatename, ctx, event_arg):
##    patmap = getattr(app.config, "perpage_html_theme_options")
##    if not patmap:
##        return
##    values = bestmatch(patmap, pagename, None, param="pagename")
##    if values:
##        for key, value in values.items():
##            # TODO: validate 'key' is valid theme option
##            ctx['theme_' + key] = value

#=============================================================================
# sphinx init
#=============================================================================
def setup(app):
    app.add_config_value('perpage_html_logo', None, 'env')
    ##app.add_config_value('perpage_html_theme_options', None, 'env')
    app.connect('html-page-context', perpage_html_logo)
    ##app.connect('html-page-context', perpage_html_theme_options)

#=============================================================================
# eof
#=============================================================================
