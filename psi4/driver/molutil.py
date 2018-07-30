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

"""Module with utility functions that act on molecule objects."""
from __future__ import absolute_import

import numpy as np

from psi4 import core
from psi4.driver.p4util import constants, filter_comments
from psi4.driver.inputparser import process_pubchem_command, pubchemre
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
def molecule_from_string(cls,
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
    molrec = qcdb.molparse.from_string(molstr=molstr,
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
def molecule_from_arrays(cls,
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

                         missing_enabled_return='error',
                         tooclose=0.1,
                         zero_ghost_fragments=False,
                         nonphysical=False,
                         mtol=1.e-3,
                         verbose=1,

                         return_dict=False):
        """Construct Molecule from unvalidated arrays and variables.

        Light wrapper around :py:func:`~qcdb.molparse.from_arrays`
        that is a full-featured constructor to dictionary representa-
        tion of Molecule. This follows one step further to return
        Molecule instance.

        Parameters
        ----------
        See :py:func:`~qcdb.molparse.from_arrays`.

        Returns
        -------
        :py:class:`psi4.core.Molecule`

        """
        molrec = qcdb.molparse.from_arrays(geom=geom,
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
def molecule_from_schema(cls, molschema, return_dict=False, verbose=1):
        """Construct Molecule from non-Psi4 schema.

        Light wrapper around :py:func:`~psi4.core.Molecule.from_arrays`.

        Parameters
        ----------
        molschema : dict
            Dictionary form of Molecule following known schema.
        return_dict : bool, optional
            Additionally return Molecule dictionary intermediate.
        verbose : int, optional
            Amount of printing.

        Returns
        -------
        mol : :py:class:`psi4.core.Molecule`
        molrec : dict, optional
            Dictionary representation of instance.
            Only provided if `return_dict` is True.

        """

        if (molschema.get('schema_name', '').startswith('qc_schema') and
            (molschema.get('schema_version', '') == 1)):
            # Lost Fields
            # -----------
            # * 'comment'
            # * 'provenance'
            ms = molschema['molecule']

            if 'fragments' in ms:
                frag_pattern = ms['fragments']
            else:
                frag_pattern = [np.arange(len(ms['symbols']))]

            dcontig = qcdb.Molecule.contiguize_from_fragment_pattern(frag_pattern,
                                                                     geom=ms['geometry'],
                                                                     elea=None,
                                                                     elez=None,
                                                                     elem=ms['symbols'],
                                                                     mass=ms.get('masses', None),
                                                                     real=ms.get('real', None),
                                                                     elbl=None,
                                                                     throw_reorder=True)

            molrec = qcdb.molparse.from_arrays(geom=dcontig['geom'],
                                               elea=None,
                                               elez=None,
                                               elem=dcontig['elem'],
                                               mass=dcontig['mass'],
                                               real=dcontig['real'],
                                               elbl=None,
                                               name=ms.get('name', None),
                                               units='Bohr',
                                               input_units_to_au=None,
                                               fix_com=ms.get('fix_com', None),
                                               fix_orientation=ms.get('fix_orientation', None),
                                               fix_symmetry=None,
                                               fragment_separators=dcontig['fragment_separators'],
                                               fragment_charges=ms.get('fragment_charges', None),
                                               fragment_multiplicities=ms.get('fragment_multiplicities', None),
                                               molecular_charge=ms.get('molecular_charge', None),
                                               molecular_multiplicity=ms.get('molecular_multiplicity', None),
                                               domain='qm',
                                               #missing_enabled_return=missing_enabled_return,
                                               #tooclose=tooclose,
                                               #zero_ghost_fragments=zero_ghost_fragments,
                                               #nonphysical=nonphysical,
                                               #mtol=mtol,
                                               verbose=verbose)

        else:
            raise ValidationError("""Schema not recognized, schema_name/schema_version: {}/{} """.format(molschema.get('schema_name', 'NA'), molschema.get('schema_version', 'NA')))

        if return_dict:
            return core.Molecule.from_dict(molrec), molrec
        else:
            return core.Molecule.from_dict(molrec)


def dynamic_variable_bind(cls):
    """Function to dynamically add extra members to
    the core.Molecule class.

    """
    cls.__setattr__ = molecule_set_attr
    cls.__getattr__ = molecule_get_attr

    # for py3-only, all these `cls.fn = qcdb.Molecule._raw_fn` can be
    #   replaced by `cls.fn = qcdb.Molecule.fn` and in qcdb/molecule.py
    #   itself, the raw, staticmethod fns can use their official names again
    #   and not need the append line at bottom of file. what I do for you,
    #   py2 ...
    cls.to_arrays = qcdb.Molecule._raw_to_arrays
    cls.to_dict = qcdb.Molecule._raw_to_dict
    cls.BFS = qcdb.Molecule._raw_BFS
    cls.B787 = qcdb.Molecule._raw_B787
    cls.scramble = qcdb.Molecule._raw_scramble
    cls.from_arrays = molecule_from_arrays
    cls.from_string = molecule_from_string
    cls.to_string = qcdb.Molecule._raw_to_string
    cls.from_schema = molecule_from_schema
    cls.to_schema = qcdb.Molecule._raw_to_schema


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
