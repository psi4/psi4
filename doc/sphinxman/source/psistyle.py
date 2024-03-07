from sphinx.pygments_styles import SphinxStyle
from pygments.style import Style
from pygments.token import Error

class PsiStyle(Style):
    """
    Turn off error highlighting for the sphinx code style, as there is psithon code that
    is not valid python.

    To change the code highlighting style, simply import something other than SphinxStyle.
    """
    background_color = SphinxStyle.background_color
    styles = SphinxStyle.styles
    styles.update({Error: ''})
