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

"""Dummy setup.py file solely for the purposes of getting an on-the-fly
computed version number into the conda recipe.

"""
import sys
from distutils.core import setup


def version_func():
    import subprocess

    command = 'python psi4/versioner.py --formatonly --format={versionlong}'
    process = subprocess.Popen(command.split(), shell=False, stdout=subprocess.PIPE)
    (out, err) = process.communicate()
    if sys.version_info >= (3, 0):
        return out.decode('utf-8').strip()
    else:
        return out.strip()

setup(
    version=version_func(),
)
