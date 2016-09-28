"""cloud_sptheme.ext.table_styling -- add directives for styling tables"""
#=============================================================================
# imports
#=============================================================================
# core
from cloud_sptheme.utils import PY3
if PY3:
    from itertools import zip_longest as izip_longest
else:
    # FIXME: not present in py25
    from itertools import izip_longest
import os
from shutil import copyfile
# site
from docutils import nodes
from docutils.parsers.rst import directives
from docutils.parsers.rst.directives.tables import RSTTable
from sphinx.builders.html import StandaloneHTMLBuilder
# pkg
from cloud_sptheme import _root_dir, is_cloud_theme
# local
__all__ = [
    "setup",
]

#=============================================================================
# constants
#=============================================================================

# name of key controlling whether css file is included
EMBED_KEY = "table_styling_embed_css"

# name of key controlling css class name
CLASS_KEY = "table_styling_class"

#=============================================================================
# field option parsers
#=============================================================================
def _split_argument_list(argument):
    if "," in argument:
        return argument.split(",")
    else:
        return argument.split()

def _parse_argument_map(argument, argmap, param):
    args = _split_argument_list(argument)
    if len(args) == 1 and all(c in argmap for c in args[0]):
        args = args[0]
    def norm(arg):
        try:
            return argmap[arg]
        except KeyError:
            raise ValueError("invalid %s: %r" % (param, arg))
    return [norm(arg) for arg in args]

_alignment_map = dict(
    l="left",
    r="right",
    c="center",
    j="justify",
    left="left",
    right="right",
    center="center",
    centered="center", # compat alias
    justify="justify",
    justified="justify", # compat alias
)

def alignment_list(argument):
    """convert into list of alignment options.
    raise ``ValueError`` if no args found, or invalid strings.
    """
    return _parse_argument_map(argument, _alignment_map, "alignment")

_bool_map = {"true": True, "t": True,
             "yes": True, "y": True,
             "false": False, "f": False,
             "no": False, "n": False,
             }

def bool_list(argument):
    """convert to list of true/false values"""
    return _parse_argument_map(argument, _bool_map, "boolean value")

def class_option_list(argument):
    """convert to list of list of classes"""
    args = _split_argument_list(argument)
    return [directives.class_option(arg) for arg in args]

_divider_map = {
    "0": "no",
    "1": "single",
    "2": "double",
    "none": "no",
    "single": "single",
    "double": "double",
    }

def divider_list(argument):
    return _parse_argument_map(argument, _divider_map, "divider style")

#=============================================================================
# replacement for table directive
#=============================================================================
class ExtendedRSTTable(RSTTable):
    # TODO: could have this auto-generate tabularcolumns directive for latex,
    #       based on alignment, dividers, and wrapping settings:
    #           |,|| for dividers
    #           l,c,r - align, nowrap
    #           L,C,R,J - align, wrap
    #       XXX: how to handle widths for tabularcolumns?

    option_spec = RSTTable.option_spec.copy()
    option_spec.update({
        # class, name already present
        ##'header-rows': directives.nonnegative_int,
        'header-columns': directives.nonnegative_int,
        # TODO: column-widths: support limited set of units (em/in/%)
        #       expressable under both css & latex
        'widths': directives.positive_int_list,
        'column-alignment': alignment_list,
        'column-wrapping': bool_list,
        'column-classes': class_option_list,
        'column-dividers': divider_list,
    })

    def run(self):
        result = RSTTable.run(self)
        if result and isinstance(result[0], nodes.table):
            self._update_table_classes(result[0])
        return result

    def _update_table_classes(self, table):
        assert isinstance(table, nodes.table)
        config = self.state.document.settings.env.config
        classes = getattr(config, CLASS_KEY)
        if classes:
            for cls in classes.split():
                table['classes'].append(cls)
        header_cols = self.options.get("header-columns") or 0
##        header_rows = self.options.get("header-rows")
        widths = self.options.get("widths")
        dividers = self.options.get("column-dividers")
        if dividers is None:
            get_divider = None
        else:
            def get_divider(idx):
                try:
                    return dividers[idx]
                except IndexError:
                    return "no"
        EMPTY = ()
        opts = (
            self.options.get("column-alignment", EMPTY),
            self.options.get("column-wrapping", EMPTY),
            self.options.get("column-classes", EMPTY),
        )
        def locate(cls):
            for child in table.children:
                if isinstance(child, cls):
                    return child
            return None
        tgroup = locate(nodes.tgroup)
        if not tgroup:
            return
        col = 0
        for child in tgroup:
            if isinstance(child, nodes.colspec):
                if widths and col < len(widths):
                    child['colwidth'] = widths[col]
                if col < header_cols:
                    child['stub'] = 1
                col += 1
                continue
            assert isinstance(child, (nodes.thead, nodes.tbody))
            for row in child.children:
                # add alignment and wrap classes to each entry (would add to
                # colspec, but html doesn't inherit much from colgroup)
                assert isinstance(row, nodes.row)
                for idx, (entry, align, wrap, clist) in \
                  enumerate(izip_longest(row, *opts)):
                    if entry is None:
                        # FIXME: make into propert rst error
                        raise ValueError("not enough columns for field options")
                    assert isinstance(entry, nodes.entry)
                    classes = entry['classes']
                    if align:
                        classes.append(align + "-align")
                    if wrap is False:
                        classes.append("nowrap")
                    if clist:
                        classes.extend(clist)
                    if get_divider:
                        classes.append(get_divider(idx) + "-left-divider")
                        classes.append(get_divider(idx+1) + "-right-divider")
        # untested - might be missing some docutils node framework bits
        ##if header_rows > 1:
        ##    thead = locate(nodes.thead)
        ##    if not thead:
        ##        thead = nodes.thead()
        ##        idx = table.children.index(tgroup)
        ##        table.children.insert(idx, thead)
        ##    thead.extend(tgroup.children[:header_rows])
        ##    del tgroup.children[:header_rows]

#=============================================================================
# patch builder to copy css file (if needed)
#=============================================================================
def prepare_builder(app):
    # make sure needed css styling gets included when building html
    builder = app.builder
    if not isinstance(builder, StandaloneHTMLBuilder):
        return
    value = getattr(app.config, EMBED_KEY)
    if value is None:
        value = not is_cloud_theme(app.config.html_theme)
    if not value:
        return

    # add custom css stylesheet
    name = "table_styling.css"
    app.add_stylesheet(name)

    # monkeypatch builder to copy over css file
    orig = builder.copy_static_files
    def wrapper():
        orig()
        source = os.path.join(_root_dir, "ext", name)
        target = os.path.join(builder.outdir, "_static", name)
        copyfile(source, target)
    builder.copy_static_files = wrapper

#=============================================================================
# register extension
#=============================================================================
def setup(app):
    # whether it will embed raw css
    app.add_config_value(EMBED_KEY, None, "html")

    # default class to add to all styled tables
    app.add_config_value(CLASS_KEY, "styled-table", "html")

    # replace existing table directive with custom one
    app.add_directive("table", ExtendedRSTTable)

    # add extra resources
    app.connect("builder-inited", prepare_builder)

#=============================================================================
# eof
#=============================================================================
