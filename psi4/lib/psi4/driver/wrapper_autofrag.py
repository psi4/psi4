#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
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

__all__ = ["auto_fragments"]

from typing import List, Optional

from psi4 import core


def auto_fragments(
    molecule: Optional[core.Molecule] = None,
    seed_atoms: Optional[List[List[int]]] = None,
) -> core.Molecule:
    r"""Detects fragments in unfragmented molecule using BFS algorithm.
    Currently only used for the WebMO implementation of SAPT.

    Parameters
    ----------
    molecule : :ref:`molecule <op_py_molecule>`, optional
        The target molecule, if not the last molecule defined.
    seed_atoms
        List of lists of atoms (0-indexed) belonging to independent fragments.
        Useful to prompt algorithm or to define intramolecular fragments through
        border atoms. Example: `[[1, 0], [2]]`

    Returns
    -------
    :py:class:`~psi4.core.Molecule` |w--w| fragmented molecule in
    Cartesian, fixed-geom (no variable values), no dummy-atom format.

    Examples
    --------
    >>> # [1] prepare unfragmented (and non-adjacent-atom) HHFF into (HF)_2 molecule ready for SAPT
    >>> molecule mol {\nH 0.0 0.0 0.0\nH 2.0 0.0 0.0\nF 0.0 1.0 0.0\nF 2.0 1.0 0.0\n}
    >>> print mol.nfragments()  # 1
    >>> fragmol = auto_fragments()
    >>> print fragmol.nfragments()  # 2

    """
    # Make sure the molecule the user provided is the active one
    if molecule is None:
        molecule = core.get_active_molecule()
    molecule.update_geometry()
    molname = molecule.name()

    frag, bmol = molecule.BFS(seed_atoms=seed_atoms, return_molecule=True)

    bmol.set_name(molname)
    bmol.print_cluster()
    core.print_out("""  Exiting auto_fragments\n""")

    return bmol
