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
"""Module with utility functions that act on molecule objects."""

from typing import Dict, Tuple, Union

import numpy as np
import qcelemental as qcel

from psi4 import core
from psi4.driver.p4util import temp_circular_import_blocker
from psi4.driver import qcdb
from psi4.driver.p4util.exceptions import *


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


@classmethod
def _molecule_from_string(cls,
                         molstr,
                         dtype=None,
                         name=None,
                         fix_com=None,
                         fix_orientation=None,
                         fix_symmetry=None,
                         return_dict=False,
                         enable_qm=True,
                         enable_efp=True,
                         missing_enabled_return_qm='none',
                         missing_enabled_return_efp='none',
                         verbose=1):
    molrec = qcel.molparse.from_string(
        molstr=molstr,
        dtype=dtype,
        name=name,
        fix_com=fix_com,
        fix_orientation=fix_orientation,
        fix_symmetry=fix_symmetry,
        return_processed=False,
        enable_qm=enable_qm,
        enable_efp=enable_efp,
        missing_enabled_return_qm=missing_enabled_return_qm,
        missing_enabled_return_efp=missing_enabled_return_efp,
        verbose=verbose)
    if return_dict:
        return core.Molecule.from_dict(molrec['qm']), molrec
    else:
        return core.Molecule.from_dict(molrec['qm'])


@classmethod
def _molecule_from_arrays(cls,
                         geom=None,
                         elea=None,
                         elez=None,
                         elem=None,
                         mass=None,
                         real=None,
                         elbl=None,
                         name=None,
                         units='Angstrom',
                         input_units_to_au=None,
                         fix_com=None,
                         fix_orientation=None,
                         fix_symmetry=None,
                         fragment_separators=None,
                         fragment_charges=None,
                         fragment_multiplicities=None,
                         molecular_charge=None,
                         molecular_multiplicity=None,
                         comment=None,
                         provenance=None,
                         connectivity=None,
                         missing_enabled_return='error',
                         tooclose=0.1,
                         zero_ghost_fragments=False,
                         nonphysical=False,
                         mtol=1.e-3,
                         verbose=1,
                         return_dict=False):
    """Construct Molecule from unvalidated arrays and variables.

    Light wrapper around :py:func:`~qcelemental.molparse.from_arrays`
    that is a full-featured constructor to dictionary representa-
    tion of Molecule. This follows one step further to return
    Molecule instance.

    Parameters
    ----------
    See :py:func:`~qcelemental.molparse.from_arrays`.

    Returns
    -------
    :py:class:`psi4.core.Molecule`

    """
    molrec = qcel.molparse.from_arrays(
        geom=geom,
        elea=elea,
        elez=elez,
        elem=elem,
        mass=mass,
        real=real,
        elbl=elbl,
        name=name,
        units=units,
        input_units_to_au=input_units_to_au,
        fix_com=fix_com,
        fix_orientation=fix_orientation,
        fix_symmetry=fix_symmetry,
        fragment_separators=fragment_separators,
        fragment_charges=fragment_charges,
        fragment_multiplicities=fragment_multiplicities,
        molecular_charge=molecular_charge,
        molecular_multiplicity=molecular_multiplicity,
        comment=comment,
        provenance=provenance,
        connectivity=connectivity,
        domain='qm',
        missing_enabled_return=missing_enabled_return,
        tooclose=tooclose,
        zero_ghost_fragments=zero_ghost_fragments,
        nonphysical=nonphysical,
        mtol=mtol,
        verbose=verbose)
    if return_dict:
        return core.Molecule.from_dict(molrec), molrec
    else:
        return core.Molecule.from_dict(molrec)


@classmethod
def _molecule_from_schema(cls, molschema: Dict, return_dict: bool = False, nonphysical: bool = False, verbose: int = 1) -> Union[core.Molecule, Tuple[core.Molecule, Dict]]:
    """Construct Molecule from non-Psi4 schema.

    Light wrapper around :py:func:`~psi4.core.Molecule.from_arrays`.

    Parameters
    ----------
    molschema
        Dictionary form of Molecule following known schema.
    return_dict
        Additionally return Molecule dictionary intermediate.
    nonphysical
        Do allow masses outside an element's natural range to pass validation?
    verbose
        Amount of printing.

    Returns
    -------
    mol : :py:class:`psi4.core.Molecule`
    molrec : dict
        Dictionary representation of instance.
        Only provided if `return_dict` is True.

    """
    molrec = qcel.molparse.from_schema(molschema, nonphysical=nonphysical, verbose=verbose)

    qmol = core.Molecule.from_dict(molrec)
    geom = np.array(molrec["geom"]).reshape((-1, 3))
    qmol._initial_cartesian = core.Matrix.from_array(geom)

    if return_dict:
        return qmol, molrec
    else:
        return qmol


def dynamic_variable_bind(cls):
    """Function to dynamically add extra members to
    the core.Molecule class.

    """
    cls.__setattr__ = molecule_set_attr
    cls.__getattr__ = molecule_get_attr

    cls.to_arrays = qcdb.Molecule.to_arrays
    cls.to_dict = qcdb.Molecule.to_dict
    cls.BFS = qcdb.Molecule.BFS
    cls.B787 = qcdb.Molecule.B787
    cls.scramble = qcdb.Molecule.scramble
    cls.from_arrays = _molecule_from_arrays
    cls.from_string = _molecule_from_string
    cls.to_string = qcdb.Molecule.to_string
    cls.from_schema = _molecule_from_schema
    cls.to_schema = qcdb.Molecule.to_schema
    cls.run_dftd3 = qcdb.Molecule.run_dftd3
    cls.run_dftd4 = qcdb.Molecule.run_dftd4
    cls.run_gcp = qcdb.Molecule.run_gcp
    cls.format_molecule_for_mol = qcdb.Molecule.format_molecule_for_mol


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
def geometry(geom: str, name: str = "default") -> core.Molecule:
    """Function to create a molecule object of name *name* from the
    geometry in string *geom*. Permitted for user use but deprecated
    in driver in favor of explicit molecule-passing. Comments within
    the string are filtered.

    """
    molrec = qcel.molparse.from_string(
        geom, enable_qm=True, missing_enabled_return_qm='minimal', enable_efp=True, missing_enabled_return_efp='none')

    molecule = core.Molecule.from_dict(molrec['qm'])
    if "geom" in molrec["qm"]:
        geom = np.array(molrec["qm"]["geom"]).reshape((-1, 3))
        if molrec["qm"]["units"] == "Angstrom":
            geom = geom / qcel.constants.bohr2angstroms
        molecule._initial_cartesian = core.Matrix.from_array(geom)
    molecule.set_name(name)

    if 'efp' in molrec:
        try:
            import pylibefp
        except ImportError as e:  # py36 ModuleNotFoundError
            raise ImportError("""Install pylibefp to use EFP functionality. `conda install pylibefp -c psi4` Or build with `-DENABLE_libefp=ON`""") from e
        #print('Using pylibefp: {} (version {})'.format(pylibefp.__file__, pylibefp.__version__))
        efpobj = pylibefp.from_dict(molrec['efp'])
        # pylibefp.core.efp rides along on molecule
        molecule.EFP = efpobj

    # Attempt to go ahead and construct the molecule
    try:
        molecule.update_geometry()
    except Exception:
        core.print_out("Molecule: geometry: Molecule is not complete, please use 'update_geometry'\n"
                       "                    once all variables are set.\n")

    activate(molecule)

    return molecule


def activate(mol: core.Molecule):
    """Function to set molecule object *mol* as the current active molecule.
    Permitted for user use but deprecated in driver in favor of explicit
    molecule-passing.

    """
    core.set_active_molecule(mol)
