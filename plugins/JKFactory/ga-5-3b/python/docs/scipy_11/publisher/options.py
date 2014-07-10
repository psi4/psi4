"""
Expose options from scipy_conf.cfg as a dictionary.
"""

__all__ = ['options']

from ConfigParser import ConfigParser
import os.path

default_filename = os.path.join(os.path.dirname(__file__),
                                '../scipy_proc.cfg')

def cfg2dict(filename=default_filename):
    """Return the content of a .ini file as a dictionary.

    """
    options = {}

    if not os.path.isfile(filename):
        print "*** Warning: Could not load config file '%s'." % filename
    else:
        cp = ConfigParser()
        cp.read(filename)
        for key in cp.options('default'):
            options[key] = cp.get('default', key)

    return options

options = cfg2dict()

