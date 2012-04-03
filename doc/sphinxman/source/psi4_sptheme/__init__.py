"""
This module contains a few small sphinx extensions.
They are mainly used to help with the generation
of BPS's own documentation, but some other projects
use them as well, so they are kept here.
"""
import re
import os.path

__version__ = "1.2"

def get_theme_dir():
    "return path to directory containing sphinx themes in this package"
    return os.path.abspath(os.path.join(__file__,os.path.pardir, "themes"))

def get_version(release):
    "derive short version string from longer release"
    return re.match("(\d+\.\d+)", release).group(1)
