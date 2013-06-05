""" A variety of utility functions which create, modify, or are in some way related to strings and text.
"""
import inspect
import string
import re
from grendel.util.aliasing import function_alias

def classname(class_or_str_or_obj):
    """ Convenience method for getting the typename as a string from the <type 'Something'> string sent back by str(type(obj))
    """
    tystr = None
    if inspect.isclass(class_or_str_or_obj):
        tystr = str(class_or_str_or_obj)
    elif isinstance(class_or_str_or_obj, basestring):
        tystr = class_or_str_or_obj
    else:
        tystr = str(type(class_or_str_or_obj))
    return re.sub(r'^.*\.([^\.]+)$', r'\1', re.sub(r'^<.*\'(.+)\'>$', r'\1', tystr))


# TODO Support for letters and/or parenthesis?
def superscript(num_or_str):
    """ Returns the unicode string superscript equivalent of `num` as a unicode string
    """
    supers = {
        '0' : u'\u2070',
        '1' : u'\u00B9',
        '2' : u'\u00B2',
        '3' : u'\u00B3',
        '4' : u'\u2074',
        '5' : u'\u2075',
        '6' : u'\u2076',
        '7' : u'\u2077',
        '8' : u'\u2078',
        '9' : u'\u2079',
    }
    minus = u'\u207B'
    i = u'\u2071'
    plus = u'\u207A'
    if (raises_error(int, num_or_str)) or (not raises_error(float, num_or_str) and not round(float(num_or_str)) == float(num_or_str)):
        raise TypeError("Only integers can be superscripts")
    ret_val = u''
    for digit in str(num_or_str):
        if digit == '-':
            ret_val += minus
        elif digit == '+':
            ret_val += plus
        elif digit == 'i':
            raise NotImplementedError("Complex numbers as superscripts are not implemented yet.")
        elif digit in supers:
            ret_val += supers[digit]
    return ret_val


def subscript(num_or_str):
    """ Returns the unicode string subscript equivalent of `num` as a unicode string
    """
    if (raises_error(int, num_or_str)) or (not raises_error(float, num_or_str) and not round(float(num_or_str)) == float(num_or_str)):
        raise TypeError("Only integers can be subscripts")
    plus = u'\u208A'
    minus = u'\u208B'
    ret_val = u''
    for digit in str(num_or_str):
        if digit == '-':
            ret_val += minus
        elif digit == '+':
            ret_val += plus
        elif digit in '0123456789':
            ret_val += unichr(0x2080 + int(digit))
    return ret_val

def make_safe_identifier(instr, replace_with='_'):
    return re.sub(r'[^A-Za-z0-9_]', replace_with, instr)

def camel_to_lower(instr, divider = "_"):
    """ Converts a CamalCase string to a lower_case_string_with_underscores
    If the `divider` argument is given, it is used instead of underscores to join words.

    Examples
    --------

    """
    new_str = instr[0].lower()
    for letter in instr[1:]:
        if letter in string.lowercase:
            new_str += letter
        elif letter in string.uppercase:
            # word boundary, add divider
            new_str += divider + letter.lower()
        elif letter in string.digits:
            new_str += letter
        else:
            raise ValueError("String {0} is not in proper CamelCase because it contains the character {1}".format(instr, letter))
    return new_str

def indented(instr, nspaces=4):
    ret_val = ""
    for line in instr.splitlines(True):
        ret_val += " "*nspaces + line
    return ret_val

def short_str(obj):
    if hasattr(obj, '__short_str__'):
        return obj.__short_str__()
    else:
        return str(obj)

def andjoin(iterable, fxn=str, oxford_comma=True):
    l = list(iterable)
    length = len(l)
    if length == 0:
        return '<empty list>'
    elif length == 1:
        return fxn(l[0])
    elif length == 2:
        return fxn(l[0]) + ' and ' + fxn(l[1])
    else:
        return ', '.join(fxn(i) for i in l[:-1]) + (',' if oxford_comma else '') + ' and ' + fxn(l[-1])

def orjoin(iterable, fxn=str, oxford_comma=True):
    l = list(iterable)
    length = len(l)
    if length == 0:
        return '<empty list>'
    elif length == 1:
        return fxn(l[0])
    elif length == 2:
        return fxn(l[0]) + ' or ' + fxn(l[1])
    else:
        return ', '.join(fxn(i) for i in l[:-1]) + (',' if oxford_comma else '') + ' or ' + fxn(l[-1])


#####################
# Dependent Imports #
#####################

from grendel.util.exceptions import raises_error

###########
# Aliases #
###########

indent = function_alias('indent', indented)
shortstr = function_alias('shortstr', short_str)

