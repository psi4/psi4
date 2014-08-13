# This code is from: http://pypi.python.org/pypi/rstex/

#!/usr/bin/python2
from docutils import utils, nodes
from docutils.core import publish_cmdline
from docutils.writers.latex2e import Writer, LaTeXTranslator
from docutils.parsers.rst import roles, Directive, directives


class InlineMath(nodes.Inline, nodes.TextElement):
    pass

class PartMath(nodes.Part, nodes.Element):
    pass

class PartLaTeX(nodes.Part, nodes.Element):
    pass

def mathEnv(math, label, type):
    if label:
        eqn_star = ''
    else:
        eqn_star = '*'

    if type in ("split", "gathered"):
        begin = "\\begin{equation%s}\n\\begin{%s}\n" % (type, eqn_star)
        end = "\\end{%s}\n\\end{equation%s}\n" % (type, eqn_star)
    else:
        begin = "\\begin{%s%s}\n" % (type, eqn_star)
        end = "\\end{%s%s}" % (type, eqn_star)
    if label:
        begin += "\\label{%s}\n" % label
    return begin + math + '\n' + end

def mathRole(role, rawtext, text, lineno, inliner, options={}, content=[]):
    latex = utils.unescape(text, restore_backslashes=True)
    return [InlineMath(latex=latex)], []

class MathDirective(Directive):
    has_content = True
    required_arguments = 0
    optional_arguments = 2
    final_argument_whitespace = True
    option_spec = {
        'type': directives.unchanged,
        'label': directives.unchanged,
    }
    def run(self):
        latex = '\n'.join(self.content)
        if self.arguments and self.arguments[0]:
            latex = self.arguments[0] + '\n\n' + latex
        node = PartMath()
        node['latex'] = latex
        node['label'] = self.options.get('label', None)
        node['type'] = self.options.get('type', "equation")
        ret = [node]
        return ret

class LaTeXDirective(Directive):
    has_content = True
    required_arguments = 0
    optional_arguments = 1
    final_argument_whitespace = True
    option_spec = {
        'usepackage': directives.unchanged
    }
    def run(self):
        latex = '\n'.join(self.content)
        if self.arguments and self.arguments[0]:
            latex = self.arguments[0] + '\n\n' + latex
        node = PartLaTeX()
        node['latex'] = latex
        node['usepackage'] = self.options.get("usepackage", "").split(",")
        ret = [node]
        return ret


roles.register_local_role("math", mathRole)
directives.register_directive("math", MathDirective)
directives.register_directive("latex", LaTeXDirective)

