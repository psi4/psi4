#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
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

import os
import atexit

from psi4 import core

# Numpy place holder for files and cleanup
numpy_files = []
def register_numpy_file(filename):
    if filename not in numpy_files:
        numpy_files.append(filename)

def clean_numpy_files():
    for nfile in numpy_files:
        os.unlink(nfile)

atexit.register(clean_numpy_files)

# Exit printing
def exit_printing():
    if _success_flag_:
        core.print_out( "\n*** Psi4 exiting successfully. Buy a developer a beer!\n")
    else:
        core.print_out( "\n*** Psi4 encountered an error. Buy a developer more coffee!\n")

_success_flag_ = True


# Working directory
_input_dir_ = os.getcwd()

def get_input_directory():
    return _input_dir_
