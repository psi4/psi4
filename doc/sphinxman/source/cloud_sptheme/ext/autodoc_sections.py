"""cloud_sptheme.ext.autodoc_sections - support ReST sections in docstrings"""
import re
import logging; log = logging.getLogger(__name__)

def indent_sections(lines, reference_prefix=''):
    "replaces any section headers with indented paragraphs"
    end = len(lines)-1
    out = []

    sections = []
    indent_char = ' ' * 4
    indent_level = 0
    SCHARS = '#*=-^"'
    section_chars = [""] #map of section char -> level (built up as scan progresses)
    def get_level(c):
        sc = section_chars[0]
        if c not in sc:
            sc = section_chars[0] = (sc+c)
        return sc.index(c)
    #FIXME: this doesn't detect double-barred sections
    def detect_section(idx):
        if idx == end:
            return None
        line = lines[idx].rstrip()
        if not line or line.lstrip() != line:
            return None
        next = lines[idx+1].rstrip()
        if next.lstrip() != next:
            return None
        for c in SCHARS:
            if next.startswith(c * len(line)):
                return c
        return None
    idx = 0
    lss = False #set to true last non-empty line was a section heading
    while idx <= end:
        line = lines[idx].rstrip()
        if not line:
            if not lss:
                out.append("")
            idx += 1
            continue
        new_char = detect_section(idx)
        if new_char:
            new_level = get_level(new_char)
            while sections and sections[-1] > new_level:
                sections.pop()
            if not sections or sections[-1] < new_level:
                sections.append(new_level)
            name = line.lower().strip().replace(" ", "-").replace("--", "-")
            indent = indent_char * (indent_level-1)
            out.extend([
                indent + ".. _%s:" % (reference_prefix + name),
                indent + ".. rst-class:: nested-section nested-section-%d" % (new_level+1,),
                "",
                indent + "%s\n" % line.rstrip(),
                ])
            idx += 2 #skip section header
            indent_level = max(0, len(sections))
            lss = True
            continue
        lss = False
        indent = indent_char * indent_level
        out.append(indent + line)
        idx += 1
    return out

def _remove_oneline(name, lines):
    #remove one-line description from top of module, if present,
    #cause we don't want it being duplicated (should already be listed in module's header)
    _title_re = re.compile(r"""
        ^ \s*
        ( {0} \s* -- \s* )?
        [a-z0-9 _."']*
        $
    """.format(re.escape(name)), re.X|re.I)
    if len(lines) > 1 and _title_re.match(lines[0]) and lines[1].strip() == '':
        del lines[:2]

def mangle_docstrings(app, what, name, obj, options, lines):
    if what == 'module':
        _remove_oneline(name, lines)
    elif what in ('class', 'exception', 'function', 'method'):
        name = "%s.%s" % (obj.__module__, obj.__name__)
        name = name.replace(".", "-").lower()
        lines[:] = indent_sections(lines, reference_prefix=name + "-")
    elif what in ('attribute',):
        pass
    else:
        #FIXME: handle other cases
        raise NotImplementedError, "unknown node: %r %r" % (what, obj)

def setup(app):
    app.connect('autodoc-process-docstring', mangle_docstrings)
