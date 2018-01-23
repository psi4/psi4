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

from __future__ import print_function
from __future__ import absolute_import
from psi4 import core


def auto_fragments(**kwargs):
    r"""Detects fragments if the user does not supply them.
    Currently only used for the WebMO implementation of SAPT.

    :returns: :py:class:`~psi4.core.Molecule`) |w--w| fragmented molecule.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :examples:

    >>> # [1] replicates with cbs() the simple model chemistry scf/cc-pVDZ: set basis cc-pVDZ energy('scf')
    >>> molecule mol {\nH 0.0 0.0 0.0\nH 2.0 0.0 0.0\nF 0.0 1.0 0.0\nF 2.0 1.0 0.0\n}
    >>> print mol.nfragments()  # 1
    >>> fragmol = auto_fragments()
    >>> print fragmol.nfragments()  # 2

    """
    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()
    molname = molecule.name()


    core.print_out("""  Exiting auto_fragments\n""")

