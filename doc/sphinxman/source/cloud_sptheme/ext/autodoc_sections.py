"""cloud_sptheme.ext.autodoc_sections - support ReST sections in docstrings"""
#=============================================================================
# imports
#=============================================================================
# core
import inspect
import re
import logging; log = logging.getLogger(__name__)
# site
# pkg
from cloud_sptheme.utils import u, patchapplier, monkeypatch
# local
__all__ = [
    "setup",
]

#=============================================================================
# internal helpers used for hacking up sphinx. almost ashamed to be doing this.
#=============================================================================
def get_caller_value(module, var, rtype=None, code=None):
    """
    helper which looks for nearest ancestor in call stack
    which occurred in module, and return specified local variable.

    :param module:
        name of module to look for on stack

    :param var:
        name of local variable to return

    :param rtype:
        do optional type-check to ensure **name** has expected type.

    :param code:
        optionally match based on function name (``frame.f_code.co_name``)
        as well as module.
    """
    frame = None
    try:
        frame = inspect.currentframe().f_back
        while True:
            if ((frame.f_globals.get("__name__") == module) and
                (not code or frame.f_code.co_name == code)):
                    break
            frame = frame.f_back
            if not frame:
                raise RuntimeError("couldn't find module=%r, code=%r in call stack" %
                                   (module, code or '<any>'))
        value = frame.f_locals[var]
        if rtype and not isinstance(value, rtype):
            raise TypeError("%s: expected a %r instance: %r" % (var, rtype, value))
        return value
    finally:
        frame = None

#=============================================================================
# autodoc monkeypatches / hacks
#=============================================================================
@patchapplier
def _patch_sphinx():
    """
    helper which monkeypatches sphinx to install some of our hooks.
    """

    #----------------------------------------------------------------------
    # patch document.note_implicit_target() to look for _modify_new_desc_section
    # attribute, as signal that it should munge up node that's passed to it,
    # to represent new description-level section, rather than document-level section.
    # this flag is then set by RSTState.new_subsection() patch, below.
    # NOTE: ideally, all this action could be done by a hook w/in
    #       RSTState.new_subsection(), but would have to modify source.
    #----------------------------------------------------------------------
    from docutils.nodes import document, make_id

    @monkeypatch(document)
    def note_implicit_target(_wrapped, self, target, *args, **kwds):
        # use default behavior unless signal flag is set
        entry = self._modify_new_desc_section
        if not entry:
            return _wrapped(self, target, *args, **kwds)
        self._modify_new_desc_section = None

        # NOTE: target should be section() node we're modifying,
        #       as we just got called from RSTState.new_subsection(),
        #       which is what sets modify_new_desc_section flag.
        #       'entry' should be last item in memo.desc_stack
        #       (see below).

        # add our custom css classes for styling
        # NOTE: adding 'section-header' class to H<N> node,
        #       so that our css rules don't have to be duplicated for every H<N> value.
        target['classes'].append("desc-section")
        target['classes'].append("desc-section-%d" % entry['level'])
        target[0]['classes'].append("section-header")

        # for duration of call, modify settings.id_prefix to include
        # decription prefix in any auto-generated ids. this helpers
        # section names remaining unique even if used between classes
        # (e.g. common names such as 'Constructor Options')
        settings = self.settings
        orig = settings.id_prefix
        try:
            if entry['prefix']:
                if orig:
                    settings.id_prefix = u("%s-%s") % (settings.id_prefix, entry['prefix'])
                else:
                    settings.id_prefix = entry['prefix']
            return _wrapped(self, target, *args, **kwds)
        finally:
            settings.id_prefix = orig

    document._modify_new_desc_section = None

    #----------------------------------------------------------------------
    # patch RSTState.new_subsection() to generate sections nested within
    # a description. It reads ``memo.desc_stack`` to determine if it's within
    # a description. If set, this attr should be a list of dicts,
    # each entry representing a nested description (e.g. an ObjectDescription)
    # whose content is being parsed, most recent should be last.
    # Each entry should be a dict containing:
    #
    # If there are no description entries active, the normal behavior is used.
    #
    # Each dict should contain the following keys:
    #   * prefix -- None, or string to use as prefix for section identifiers.
    #               helps keep links unique w/in document.
    #   * signode -- signature node used to generate prefix (for debugging)
    #   * owner -- arbitary object (ObjectDescription in our case)
    #              which added this entry to list. intended as sanity check
    #              when popping entries back off stack.
    #   * level -- section level w/in declaration (autoset by code below)
    #----------------------------------------------------------------------
    from docutils.parsers.rst.states import RSTState

    @monkeypatch(RSTState)
    def new_subsection(_wrapped, self, *args, **kwds):
        desc_stack = getattr(self.memo, "desc_stack", None)
        if desc_stack:
            # after new_subsection() creates section node,
            # it will invoke document.note_implicit_target().
            # setting this attr signals our monkeypatch of that method (above)
            # to make changes to that node based on desc_stack entry.
            # NOTE: ideally, the note_implicit_target() monkeypatch,
            #       as well as this code, would be placed inside RSTState.new_subsection(),
            #       but that would require modifying sphinx's source :(
            entry = desc_stack[-1]
            entry['level'] = self.memo.section_level+1 # set level w/in description
            self.document._modify_new_desc_section = entry # enable note hack

        # hand off to real method
        return _wrapped(self, *args, **kwds)

    #----------------------------------------------------------------------
    # monkeypatch ObjectDescription.run() so that:
    # 1. before calling state.nested_parse(), we push a desc context
    #    onto state_machine.desc_context_stack (see above).
    # 2. when it does call state.nested_parse() the first time,
    #    ``match_titles=True`` gets set.
    #    FIXME: using a really awkward way to accomplish this :|
    # 3. pop our context off desc_context_stack when done.
    #----------------------------------------------------------------------
    from sphinx.directives import ObjectDescription
    from sphinx.addnodes import desc as DescNodeType

    @monkeypatch(ObjectDescription)
    def run(_wrapped_run, self):
        # ObjectDescription.before_content() will be invoked right before
        # run calls ``self.state.nested_parse()``. We take advantage of that
        # by wrapping the instance's before_content() call so that the last
        # thing is does is set up ``self.state`` the way we want.

        # NOTE: have to patch per-instance, since subclasses that override this
        #       don't tend to invoke super()
        @monkeypatch(self)
        def before_content(_wrapped_before):
            # let real method do all the setup it wants.
            _wrapped_before()

            #----------------------------------------------------------------------
            # need to figure out prefix to prepend to our description sections,
            # so their IDs are unique. for now, using signature node
            # ObjectDescription.run() has finished generating right before before_content()
            # was called. Unfortunately, it's not available via self,
            # so we have to reach into call stack to grab it...
            # NOTE: ideally we would do this in ObjectDescription.run().
            #----------------------------------------------------------------------
            node = get_caller_value("sphinx.directives", "node",
                                    rtype=DescNodeType, code="run")
            # FIXME: would like a more bullet-proof way of deriving our id prefix...
            signode = node.children[0]
            if signode.get("ids"):
                base = signode['ids'][0]
            elif signode.get("names"):
                base = signode['names'][0]
            else:
                base = signode.astext()
            prefix = re.sub("[^a-zA-Z0-9_.]+", "-", base).strip("-") + "-"

            # now that we've got that info, add our description context entry
            # to the stack
            memo = self.state.memo
            if not hasattr(memo, "desc_stack"):
                memo.desc_stack = []
            memo.desc_stack.append(dict( # see new_subsection() above for dict format
                prefix=prefix,
                owner=self,
                signode=signode,
                level=0,
            ))

            #----------------------------------------------------------------------
            # hack up ``state.nested_parse()`` so that the next time it's called,
            # 'match_titles=True' is set. that call should happen as soon as this
            # function returns back to ObjectDescription.run()
            # NOTE: ideally, we would just set match_titles=True within ObjectDescription.run()
            #----------------------------------------------------------------------
            state = self.state
            if not hasattr(state, "_set_next_match_titles_flag"):
                # state is persistent object, only want to patch it once.
                @monkeypatch(state)
                def nested_parse(_wrapped_parse, *args, **kwds):
                    if state._set_next_match_titles_flag:
                        kwds['match_titles'] = True
                        state._set_next_match_titles_flag = False
                    return _wrapped_parse(*args, **kwds)

            # signal our hack to set match_titles
            state._set_next_match_titles_flag = True

        @monkeypatch(self)
        def after_content(_wrapped_after):
            # remove our description context entry
            # NOTE: ideally would do this in ObjectDescription.run()
            desc = self.state.memo.desc_stack.pop()
            assert desc['owner'] is self, "sanity check failed"

            # let real method do it's work
            return _wrapped_after()

        # now invoke the real run() method.
        # as soon as it calls before_content(), our hack above will patch self.state.
        # after before_content() returns, real run method will call self.state.nested_parse(),
        # and invoke our patched version instead.
        return _wrapped_run(self)

    #----------------------------------------------------------------------
    # make autodoc invoke parse_nested_section_with_titles() for ALL objects
    # if this isn't done, autodoc generates paragraphs instead of sections.
    # this causes all nested content to be omitted
    # FIXME: why is the lack of this causing a problem? should track it down.
    #----------------------------------------------------------------------
    from sphinx.ext.autodoc import Documenter
    Documenter.titles_allowed = True

    #----------------------------------------------------------------------
    # finally, monkeypatch DocFieldTransformer.transform_all()
    # so that it transforms doc fields  within one of our nested sections
    # (default code only looks at top-level nodes)
    #
    # FIXME: find a cleaner way to do this :|
    #----------------------------------------------------------------------
    from sphinx.util.docfields import DocFieldTransformer
    from docutils.nodes import section

    @monkeypatch(DocFieldTransformer)
    def transform_all(_wrapped, self, node):
        # transform immediate node contents like normal
        _wrapped(self, node)

        # our nested sections show up as definition lists,
        # so make sure transform_all is also invoked for the contents
        # of any definition list
        for child in node:
            if isinstance(child, section):
                _wrapped(self, child.children)

    #----------------------------------------------------------------------
    # sigh. done monkeypatching.
    #----------------------------------------------------------------------

#=============================================================================
# docstring mangling
#=============================================================================
def trim_module_header(app, what, name, obj, options, lines):
    """
    helper to remove one-line description from top of module (if preset).
    """
    if what != "module":
        return
    _title_re = re.compile(r"""
        ^ \s*
        ( {0} \s* -- \s* )?
        [a-z0-9 _."']*
        $
    """.format(re.escape(name)), re.X|re.I)
    if len(lines) > 1 and _title_re.match(lines[0]) and lines[1].strip() == '':
        del lines[:2]

#=============================================================================
# sphinx extension entrypoint
#=============================================================================
def setup(app):
    # don't patch sphinx unless this extension is actually in use
    _patch_sphinx()

    # clean up leading bit of module docstring
    app.connect('autodoc-process-docstring', trim_module_header)

#=============================================================================
# documentation helper
#
# NOTE: this function doesn't actually do anything,
#       it exists to test this extension's behavior as part of docs/cloud_theme_test
#=============================================================================
def _doctestfunc():
    """
    The :mod:`~cloud_sptheme.ext.autodoc_sections` extension should generate
    nested sections as found within object docstrings.

    Nested Section
    ==============

    :param arg: xxx

    .. attribute:: foo

        bar

    These sections can in turn contain others:

    Child Section
    -------------

    Which allows breaking long class docstrings up in meaningful ways.

    Child Section 2
    ---------------

    And more content

    Nested Section 2
    ================

    end of class
    """
    pass

#=============================================================================
# eof
#=============================================================================
