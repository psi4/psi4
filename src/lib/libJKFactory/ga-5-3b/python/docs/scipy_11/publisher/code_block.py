# --- Code-block directive from Sphinx ---

from docutils import nodes
from docutils.parsers.rst import Directive, directives

class CodeBlock(Directive):
    """
    Directive for a code block with special highlighting or line numbering
    settings.
    """

    has_content = True
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec = {
        'linenos': directives.flag,
        'linenostart': directives.nonnegative_int,
    }

    def run(self):
        code = u'\n'.join(self.content)
        literal = nodes.literal_block(code, code)
        literal['language'] = self.arguments[0]
        literal['linenos'] = 'linenos' in self.options
        literal['linenostart'] = self.options.get('linenostart', 1)
        return [literal]

directives.register_directive('code-block', CodeBlock)

# --- End code-block directive from Sphinx ---
