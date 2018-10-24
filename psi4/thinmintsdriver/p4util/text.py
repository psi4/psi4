#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2018 The Psi4 Developers.
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

from __future__ import print_function
"""Module with utility classes and functions related
to data tables and text.

"""
import sys

from thinmints.driver import constants
from .exceptions import *


def print_stdout(stuff):
    """Function to print *stuff* to standard output stream."""
    print(stuff, file=sys.stdout)


def print_stderr(stuff):
    """Function to print *stuff* to standard error stream."""
    print(stuff, file=sys.stderr)

def levenshtein(seq1, seq2):
    """Function to compute the Levenshtein distance between two strings."""
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

def find_approximate_string_matches(seq1,options,max_distance):
    """Function to compute approximate string matches from a list of options."""
    matches = []
    for seq2 in options:
        distance = levenshtein(seq1,seq2)
        if distance <= max_distance:
            matches.append(seq2)
    return matches
