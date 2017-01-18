#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import fcntl
import os
import sys
import struct
import termios
import textwrap

if sys.version_info[0] == 2:
    from StringIO import StringIO
elif sys.version_info[0] > 2:
    from io import StringIO

from .color import *

_debug = False
_verbose = False
indent = "  "


def is_verbose():
    return _verbose


def is_debug():
    return _debug


def set_verbose(flag):
    global _verbose
    _verbose = flag


def set_debug(flag):
    global _debug
    _debug = False


def msg(message, *args):
    cprint("@*b{==>} %s" % cescape(message))
    for arg in args:
        print(indent + str(arg))


def info(message, *args, **kwargs):
    fmt = kwargs.get('format', '*b')
    stream = kwargs.get('stream', sys.stdout)
    wrap = kwargs.get('wrap', False)

    cprint("@%s{==>} %s" % (fmt, cescape(str(message))), stream=stream)
    for arg in args:
        if wrap:
            lines = textwrap.wrap(
                str(arg), initial_indent=indent, subsequent_indent=indent
            )
            for line in lines:
                stream.write(line + '\n')
        else:
            stream.write(indent + str(arg) + '\n')


def verbose(message, *args, **kwargs):
    if _verbose:
        kwargs.setdefault('format', 'c')
        info(message, *args, **kwargs)


def debug(message, *args, **kwargs):
    if _debug:
        kwargs.setdefault('format', 'g')
        kwargs.setdefault('stream', sys.stderr)
        info(message, *args, **kwargs)


def error(message, *args, **kwargs):
    kwargs.setdefault('format', '*r')
    kwargs.setdefault('stream', sys.stderr)
    info("Error: " + str(message), *args, **kwargs)


def warn(message, *args, **kwargs):
    kwargs.setdefault('format', '*Y')
    kwargs.setdefault('stream', sys.stderr)
    info("Warning: " + str(message), *args, **kwargs)


def die(message, *args, **kwargs):
    error(message, *args, **kwargs)
    sys.exit(1)


def hline(label=None, **kwargs):
    """Draw a labeled horizontal line.
       Options:
           char       Char to draw the line with. Default '-'
           max_width  Maximum width of the line. Default is 64 chars.
    """
    char = kwargs.pop('char', '-')
    max_width = kwargs.pop('max_width', 64)
    if kwargs:
        raise TypeError("'%s' is an invalid keyword argument for this function."
                        % next(kwargs.iterkeys()))

    rows, cols = terminal_size()
    if not cols:
        cols = max_width
    else:
        cols -= 2
    cols = min(max_width, cols)

    label = str(label)
    prefix = char * 2 + " "
    suffix = " " + (cols - len(prefix) - clen(label)) * char

    out = StringIO()
    out.write(prefix)
    out.write(label)
    out.write(suffix)

    print(out.getvalue())


def terminal_size():
    """Gets the dimensions of the console: (rows, cols)."""

    def ioctl_GWINSZ(fd):
        try:
            rc = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ, '1234'))
        except:
            return
        return rc

    rc = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not rc:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            rc = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not rc:
        rc = (os.environ.get('LINES', 25), os.environ.get('COLUMNS', 80))

    return int(rc[0]), int(rc[1])
