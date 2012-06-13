"""cloud_sptheme.ext.escaped_samp_literals - allow escaping { and } in samp."""
from docutils import nodes, utils

def emph_literal_role(typ, rawtext, text, lineno, inliner,
                      options={}, content=[]):
    """replacement for sphinx's ``:samp:`` role handler.
    this is a bit stricter in it's parsing, and allows escaping of literal
    ``{`` and ``}`` characters.
    """
    def make_error(pos, value):
        value = "%s at char %d of %s" % (value, pos, rawtext)
        msg = inliner.reporter.error(value, line=lineno)
        prb = inliner.problematic(rawtext, rawtext, msg)
        return [prb], [msg]
    text = utils.unescape(text)
    retnode = nodes.literal(role=typ.lower(), classes=[typ])
    buffer = u"" # contains text being accumulated for next node
    in_escape = False # True if next char is part of escape sequence
    in_var = False # True if parsing variable section instead of plain text
    var_start = None # marks start of var section if in_var is True
    i = 0
    for c in text:
        i += 1
        if in_escape:
            # parse escape sequence
            if c in u"{}\\":
                buffer += c
                in_escape = False
            else:
                return make_error(i-2, "unknown samp-escape '\\\\%s'" % (c,))
        elif c == u"\\":
            # begin escape sequence
            in_escape = True
            i += 1 # account for extra escape char in rawtext
        elif in_var:
            # parsing variable section
            if c == u"{":
                return make_error(i, "unescaped '{'")
            elif c == u"}":
                # finalize variable section, return to plaintext
                if not buffer:
                    return make_error(i-1, "empty variable section")
                retnode += nodes.emphasis(buffer, buffer)
                buffer = u""
                in_var = False
            else:
                buffer += c
        else:
            # parsing plaintext section
            if c == u"{":
                # finalize plaintext section, start variable section
                if buffer:
                    retnode += nodes.Text(buffer, buffer)
                buffer = u""
                in_var = True
                var_start = i
            elif c == u"}":
                return make_error(i, "unescaped '}'")
            else:
                buffer += c
    if in_escape:
        return make_error(i, "unterminated samp-escape sequence")
    elif in_var:
        return make_error(var_start, "unterminated variable section")
    elif buffer:
        retnode += nodes.Text(buffer, buffer)
    return [retnode], []

def setup(app):
    # register our handler to overrride sphinx.roles.emph_literal_role
    from docutils.parsers.rst import roles
    import sphinx.roles as mod
    names = [
        key for key,value in mod.specific_docroles.items()
        if value is mod.emph_literal_role
    ]
    for name in names:
        roles.register_local_role(name, emph_literal_role)
