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

"""
This module contains a few small sphinx extensions.
They are mainly used to help with the generation
of BPS's own documentation, but some other projects
use them as well, so they are kept here.
"""
import re
import os.path

__version__ = "1.4"

def get_theme_dir():
    "return path to directory containing sphinx themes in this package"
    return os.path.abspath(os.path.join(__file__,os.path.pardir, "themes"))

def get_version(release):
    "derive short version string from longer release"
    return re.match("(\d+\.\d+)", release).group(1)
