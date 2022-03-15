#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2021 The Psi4 Developers.
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

import os

from psi4 import core

from psi4.driver.p4util.exceptions import ValidationError


def cubeprop_compute_properties(self):
    """Filesystem wrapper for CubeProperties::raw_compute_properties."""

    filepath = core.get_global_option("CUBEPROP_FILEPATH")

    # Is filepath a valid directory?
    if not os.path.isdir(os.path.abspath(os.path.expandvars(filepath))):
        raise ValidationError("""Filepath "{}" is not valid.  Please create this directory.""".format(filepath))

    geomfile = filepath + os.sep + 'geom.xyz'
    xyz = self.basisset().molecule().to_string(dtype='xyz', units='Angstrom')
    with open(geomfile, 'w') as fh:
        fh.write(xyz)

    self.raw_compute_properties()


core.CubeProperties.compute_properties = cubeprop_compute_properties
