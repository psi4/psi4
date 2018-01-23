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

"""Module with utility functions that act on molecule objects."""
from __future__ import absolute_import

from psi4 import core
from psi4.driver.p4util import constants, filter_comments
from psi4.driver.inputparser import process_pubchem_command, pubchemre
from psi4.driver import qcdb


def molecule_set_attr(self, name, value):
    """Function to redefine __setattr__ method of molecule class."""
    fxn = object.__getattribute__(self, "is_variable")
    isvar = fxn(name)
    if isvar:
        fxn = object.__getattribute__(self, "set_variable")
        fxn(name, value)
        return

    object.__setattr__(self, name, value)


def molecule_get_attr(self, name):
    """Function to redefine __getattr__ method of molecule class."""
    fxn = object.__getattribute__(self, "is_variable")
    isvar = fxn(name)

    if isvar:
        fxn = object.__getattribute__(self, "get_variable")
        return fxn(name)

    return object.__getattribute__(self, name)




def dynamic_variable_bind(cls):
    """Function to dynamically add extra members to
    the core.Molecule class.

    """
    cls.__setattr__ = molecule_set_attr
    cls.__getattr__ = molecule_get_attr

    cls.BFS = BFS


dynamic_variable_bind(core.Molecule)  # pass class type, not class instance


#
# Define geometry to be used by PSI4.
# The molecule created by this will be set in options.
#
# geometry("
#   O  1.0 0.0 0.0
#   H  0.0 1.0 0.0
#   H  0.0 0.0 0.0
#
def geometry(geom, name="default"):
    """Function to create a molecule object of name *name* from the
    geometry in string *geom*. Permitted for user use but deprecated
    in driver in favor of explicit molecule-passing. Comments within
    the string are filtered.

    """
    core.efp_init()
    geom = pubchemre.sub(process_pubchem_command, geom)
    geom = filter_comments(geom)
    molecule = core.Molecule.create_molecule_from_string(geom)
    molecule.set_name(name)

    # Attempt to go ahead and construct the molecule
    try:
        molecule.update_geometry()
    except:
        core.print_out("Molecule: geometry: Molecule is not complete, please use 'update_geometry'\n"
                       "                    once all variables are set.\n")

    activate(molecule)

    return molecule


def activate(mol):
    """Function to set molecule object *mol* as the current active molecule.
    Permitted for user use but deprecated in driver in favor of explicit
    molecule-passing.

    """
    core.set_active_molecule(mol)
