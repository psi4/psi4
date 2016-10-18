"""
cloud_sptheme.ext.docfield_markup -- monkeypatches sphinx to allow ``~`` in docfields.
"""
#=============================================================================
# imports
#=============================================================================
# core
import logging; log = logging.getLogger(__name__)
# site
# pkg
from cloud_sptheme.utils import patchapplier, monkeypatch
# local
__all__ = [
    "setup",
]

#=============================================================================
# patch
#=============================================================================
@patchapplier
def _patch_docfield():
    from sphinx.util.docfields import Field, nodes, addnodes

    # NOTE: would like to just wrap make_xref(), but have to override so
    #       many parts that just end up replicating all the code :(
    #       hence this ignored _wrapped() entirely

    @monkeypatch(Field)
    def make_xref(_wrapped, self, rolename, domain, target, innernode=nodes.emphasis, contnode=None):
        # 'contnode' argument added in sphinx 1.3 
        # not sure what it does, aborting if feature is present
        if contnode:
            return _wrapped(self, rolename, domain, target, innernode=innernode, contnode=contnode)
        # NOTE: this tries to replicate the convention used in PyXRefRole.process_link()
        rawtext = title = target
        if rawtext.startswith("~"):
            # if the first character is a tilde, don't display the module/class
            # parts of the contents
            target = target.lstrip("~")
            title = target.rpartition(".")[2]
        if issubclass(innernode, nodes.Text):
            # Text classes want rawtext second
            node = innernode(title, rawtext)
        else:
            # Element classes want rawtext first
            node = innernode(rawtext, title)
        if not rolename:
            return node
        refnode = addnodes.pending_xref(title, refdomain=domain, refexplicit=True,
                                        reftype=rolename, reftarget=target)
        refnode += node
        return refnode

#=============================================================================
# sphinx entrypoint
#=============================================================================
def setup(app):
    # don't apply our patch unless actually loaded by sphinx
    _patch_docfield()

#=============================================================================
# eoc
#=============================================================================
