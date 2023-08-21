#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2023 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
"""Module with utility classes and functions related
to data tables and text.

"""

__all__ = [
    "banner",
    "find_approximate_string_matches",
    "levenshtein",
    "message_box",
]

import sys
import warnings
from typing import List, Optional

from psi4 import core

from .. import constants


def banner(text: str, type: int = 1, width: int = 35, strNotOutfile: bool = False) -> Optional[str]:
    """Format `text` into a banner style and print or return it.

    Parameters
    ----------
    text
        String to be emphasized.
    type
        Style 1 has minimum three-line height. Style 2 has minimum one-light
        height.
    width
        Minimum length of banner string.
    strNotOutfile
        Controls mode of return.

    Returns
    -------
    str
        If *strNotOutfile* is True, return string.
    None
        If *strNotOutfile* is False, print it to output file.

    """
    lines = text.split('\n')
    max_length = 0
    for line in lines:
        max_length = max(len(line), max_length)

    max_length = max(width, max_length)

    null = ''
    if type == 1:
        banner = '  //' + null.center(max_length, '>') + '//\n'
        for line in lines:
            banner += '  //' + line.center(max_length) + '//\n'
        banner += '  //' + null.center(max_length, '<') + '//\n'

    if type == 2:
        banner = ''
        for line in lines:
            banner += (' ' + line + ' ').center(max_length, '=')

    if strNotOutfile:
        return banner
    else:
        core.print_out(banner)


def levenshtein(seq1: str, seq2:str) -> int:
    """Compute the Levenshtein distance between two strings.

    Parameters
    ----------
    seq1
        First string.
    seq2
        Second string.

    """
    oneago = None
    thisrow = list(range(1, len(seq2) + 1)) + [0]
    for x in range(len(seq1)):
        twoago, oneago, thisrow = oneago, thisrow, [0] * len(seq2) + [x + 1]
        for y in range(len(seq2)):
            delcost = oneago[y] + 1
            addcost = thisrow[y - 1] + 1
            subcost = oneago[y - 1] + (seq1[x] != seq2[y])
            thisrow[y] = min(delcost, addcost, subcost)
    return thisrow[len(seq2) - 1]


def find_approximate_string_matches(seq1: str, options: List[str], max_distance: int) -> List[str]:
    """Find list of approximate (within `max_distance`) matches to string `seq1` among `options`.

    Parameters
    ----------
    seq1
        Target string to look for near matches to.
    options
        Alternatives among which to look for near matches to `seq1`.
    max_distance
        Maximum Levenshtein distance from `seq1` to return.

    """
    return [seq2 for seq2 in options if (levenshtein(seq1, seq2) <= max_distance)]


def message_box(message: str, max_width: int = 80, min_width: int = 30) -> str:
    """Put a message string into a box for extra attention.

    Parameters
    ----------
    message
       Message string to be boxed.
    max_width
        Maximal character width of the box.
    min_width
        Minimal character width of the box.

    Returns
    -------
    str
       Box containing the message as a multiline string.
    """
    from textwrap import wrap

    # ensure box is within min/max boundaries
    msg = message.splitlines()
    max_line = len(max(msg, key=len))
    box_width = max(min(max_width, max_line), min_width)

    error_str = []
    error_str.append('\n!' + '-' * box_width + '--!\n')
    error_str.append('!' + ' ' * box_width + '  !\n')

    fmt = "! {:" + str(box_width) + "} !\n"
    for line in msg[:]:
        error_str.extend([fmt.format(x) for x in wrap(line, box_width, subsequent_indent="    ")])

    error_str.append('!' + ' ' * box_width + '  !\n')
    error_str.append('!' + '-' * box_width + '--!\n')
    error_str = ''.join(error_str)

    return error_str
