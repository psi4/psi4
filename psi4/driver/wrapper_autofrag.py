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

from __future__ import print_function
from __future__ import absolute_import
import math
from psi4 import core

def _autofragment_convert(p, symbol):
    # Finding radii for auto-fragmenter
    if symbol[p] == 'H':
        d = 1.001
    if symbol[p] == 'He':
        d = 1.012
    if symbol[p] == 'Li':
        d = 0.825
    if symbol[p] == 'Be':
        d = 1.408
    if symbol[p] == 'B':
        d = 1.485
    if symbol[p] == 'C':
        d = 1.452
    if symbol[p] == 'N':
        d = 1.397
    if symbol[p] == 'O':
        d = 1.342
    if symbol[p] == 'F':
        d = 1.287
    if symbol[p] == 'Ne':
        d = 1.243
    if symbol[p] == 'Na':
        d = 1.144
    if symbol[p] == 'Mg':
        d = 1.364
    if symbol[p] == 'Al':
        d = 1.639
    if symbol[p] == 'Si':
        d = 1.716
    if symbol[p] == 'P':
        d = 1.705
    if symbol[p] == 'S':
        d = 1.683
    if symbol[p] == 'Cl':
        d = 1.639
    if symbol[p] == 'Ar':
        d = 1.595

    return d / 1.5


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

    geom = molecule.save_string_xyz()

    numatoms = molecule.natom()
    VdW = [1.2, 1.7, 1.5, 1.55, 1.52, 1.9, 1.85, 1.8]

    symbol = list(range(numatoms))
    X = [0.0] * numatoms
    Y = [0.0] * numatoms
    Z = [0.0] * numatoms

    Queue = []
    White = []
    Black = []
    F = geom.split('\n')
    for f in range(numatoms):
        A = F[f + 1].split()
        symbol[f] = A[0]
        X[f] = float(A[1])
        Y[f] = float(A[2])
        Z[f] = float(A[3])
        White.append(f)
    Fragment = [[] for i in range(numatoms)]  # stores fragments

    start = 0  # starts with the first atom in the list
    Queue.append(start)
    White.remove(start)

    frag = 0

    while((len(White) > 0) or (len(Queue) > 0)):  # Iterates to the next fragment
        while(len(Queue) > 0):  # BFS within a fragment
            for u in Queue:  # find all nearest Neighbors
                             #   (still coloured white) to vertex u
                for i in White:
                    Distance = math.sqrt((X[i] - X[u]) * (X[i] - X[u]) +
                                         (Y[i] - Y[u]) * (Y[i] - Y[u]) +
                                         (Z[i] - Z[u]) * (Z[i] - Z[u]))
                    if Distance < _autofragment_convert(u, symbol) + _autofragment_convert(i, symbol):
                        Queue.append(i)  # if you find you, put it in the que
                        White.remove(i)  # and remove it from the untouched list
            Queue.remove(u)  # remove focus from Queue
            Black.append(u)
            Fragment[frag].append(int(u + 1))  # add to group (adding 1 to start
                                           #   list at one instead of zero)

        if(len(White) != 0):  # cant move White->Queue if no more exist
            Queue.append(White[0])
            White.remove(White[0])
        frag += 1

    new_geom = """\n"""
    for i in Fragment[0]:
        new_geom = new_geom + F[i].lstrip() + """\n"""
    new_geom = new_geom + """--\n"""
    for j in Fragment[1]:
        new_geom = new_geom + F[j].lstrip() + """\n"""
    new_geom = new_geom + """units angstrom\n"""

    moleculenew = core.Molecule.create_molecule_from_string(new_geom)
    moleculenew.set_name(molname)
    moleculenew.update_geometry()
    moleculenew.print_cluster()
    core.print_out("""  Exiting auto_fragments\n""")

    return moleculenew

