from pygments.style import Style
from pygments.styles.friendly import FriendlyStyle
from pygments.token import Generic, Comment, Number

class SphinxStyle(Style):
    """
    Like friendly, but a bit darker to enhance contrast on the green
    background.
    """

    background_color = '#eeffcc'
    default_style = ''

    styles = FriendlyStyle.styles
    styles.update({
        Generic.Output: '#333',
        Comment: 'italic #408090',
        Number: '#208050',
    })
