"""
cloud_sptheme.ext.autoattribute_search_bases -- monkeypatches autodoc so ``autoattribute`` searches base classes for attr doc.
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
def _patch_autoattribute():
    from sphinx.ext.autodoc import AttributeDocumenter, ModuleAnalyzer, PycodeError

    @monkeypatch(AttributeDocumenter)
    def add_content(_wrapped, self, *args, **kwds):
        if not self._datadescriptor and self.analyzer and self.objpath:
            attr_docs = self.analyzer.find_attr_docs()
            key = ('.'.join(self.objpath[:-1]), self.objpath[-1])
            if key not in attr_docs:
                # look for parent class w/ correct attr
                if hasattr(self.parent, "__mro__"):
                    for basecls in self.parent.__mro__[1:]:
                        try:
                            analyzer = ModuleAnalyzer.for_module(basecls.__module__)
                            base_attr_docs = analyzer.find_attr_docs()
                        except PycodeError as err:
                            continue
                        # FIXME: need qualname or equivalent for basecls
                        base_key = (basecls.__name__, self.objpath[-1])
                        if base_key in base_attr_docs:
                            # insert data into existing analyzer
                            # XXX: might be prettier way to handle this,
                            #      (esp so actual source file was reported)
                            #      but would have to duplicate much of add_content()
                            attr_docs[key] = base_attr_docs[base_key]
                            break
        return _wrapped(self, *args, **kwds)

#=============================================================================
# sphinx entry point
#=============================================================================
def setup(app):
    # don't apply our patch unless actually loaded by sphinx
    _patch_autoattribute()

#=============================================================================
# eoc
#=============================================================================
