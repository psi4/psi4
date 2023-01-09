#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2022 The Psi4 Developers.
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

import math
import collections

import numpy as np

import qcelemental as qcel


def BFS(geom, elem, seed_atoms=None, bond_threshold=1.20):
    """Detect fragments among real atoms through a breadth-first search (BFS) algorithm.

    Parameters
    ----------
    geom : ndarray of float
        (nat x 3) Cartesian coordinates [a0] of real atoms.
    elem : ndarray of str or int
        (nat) Either element symbols or atomic numbers corresponding to `geom`.
        Used for selecting van der Waals radius.
    seed_atoms : list (optional)
        List of lists of atoms (0-indexed) belonging to independent fragments.
        Useful to prompt algorithm or to define intramolecular fragments through
        border atoms. Example: `[[1, 0], [2]]`.
    bond_threshold : float (optional)
        Factor beyond average of covalent radii to determine bond cutoff.

    Returns
    -------
    list of lists
        Array of atom indices (0-indexed) of detected fragments. See example
        below for how to transform inputs.

    Notes
    -----
    Relies upon van der Waals radii and so faulty for close (especially
    hydrogen-bonded) fragments. `seed_atoms` can help.

    Authors
    -------
    Original code from Michael S. Marshall, linear-scaling algorithm from
    Trent M. Parker, revamped by Lori A. Burns

    Usage
    -----
    >>> # [1] BFS on large array of jumbled coordinates `geom` and element
    >>> #     symbols `elem`. Use the output `fragments` to form list of small
    >>> #     per-fragment arrays.
    >>> fragments = BFS(geom, elem)
    >>> frag_geoms = [geom[fr] for fr in fragments]
    >>> frag_elems = [elem[fr] for fr in fragments]

    """
    radii = _get_covalent_radii(elem)
    max_covalent_radius = np.max(radii)
    blocksize = int(math.ceil(2.0 * bond_threshold * max_covalent_radius))
    allblocks = _get_blocks(geom, blocksize)

    bond_tree = _get_bond_tree(radii, geom, allblocks, blocksize, bond_threshold)
    if seed_atoms is None:
        seed_atoms = []
    allfragments = seed_atoms

    # bare queues
    new_list = []
    break_list = []
    unfound_list = list(range(geom.shape[0]))

    # seed queues from intrafrag atom hints
    for ifr, fr in enumerate(allfragments):
        new_list.append([])
        for at in fr:
            new_list[ifr].append(at)
            break_list.append(at)
            unfound_list.remove(at)

    # perform BFS
    while len(unfound_list) > 0:
        for ifr, fr in enumerate(new_list):
            while len(fr) > 0:
                for at1 in reversed(fr):
                    for at2 in bond_tree[at1]:
                        if at2 in unfound_list and at2 not in break_list:
                            allfragments[ifr].append(at2)
                            new_list[ifr].append(at2)
                            unfound_list.remove(at2)
                    new_list[ifr].remove(at1)
        if len(unfound_list) > 0:
            at_new = unfound_list[0]
            allfragments.append([at_new])
            new_list.append([at_new])
            unfound_list.remove(at_new)

    for fr in range(len(allfragments)):
        allfragments[fr] = sorted(allfragments[fr])

    return allfragments


def _get_covalent_radii(elem):
    """Return covalent radii [a0] for all atoms

    Look-up values for covalent (or ionic) radii by atomic element [A] from
    "Inorganic Chemistry" 3rd ed, Housecroft, Appendix 6, pgs 1013-1014

    """
    covalent_radii_lookup = {
        'H' : 0.37,                                                                                     'He': 0.30,
        'Li': 1.02, 'Be': 0.27,             'B' : 0.88, 'C' : 0.77, 'O' : 0.73, 'N' : 0.75, 'F' : 0.71, 'Ne': 0.84,
        'Na': 1.02, 'Mg': 0.72,             'Al': 1.30, 'Si': 1.18, 'P' : 1.10, 'S' : 1.03, 'Cl': 0.99, 'Ar': 1.00,
        'K' : 1.38, 'Ca': 1.00,
                                'Sc': 0.75, 'Ti': 0.86, 'V' : 0.79, 'Cr': 0.73, 'Mn': 0.67,
                                'Fe': 0.61, 'Co': 0.64, 'Ni': 0.55, 'Cu': 0.46, 'Zn': 0.60,
                                            'Ga': 1.22, 'Ge': 1.22, 'As': 1.22, 'Se': 1.17, 'Br': 1.14, 'Kr': 1.03,
                                                                                            'I' : 1.33,
                                                                                                        'X' : 0.00}  # yapf: disable
    #'RN': 2.40 / 1.5,  # extrapolation
    #'H': 1.06 / 1.5,  # Bondi JPC 68 441 (1964)
    #'SN': 2.16 / 1.5,  # Bondi JPC 68 441 (1964)
    #'SB': 2.12 / 1.5,  # Bondi JPC 68 441 (1964)
    #'TE': 2.08 / 1.5,  # Bondi JPC 68 441 (1964)
    #'XE': 2.05 / 1.5}  # Bondi JPC 68 441 (1964)
    nat = elem.shape[0]
    try:
        caps = [el.capitalize() for el in elem]
    except AttributeError:
        caps = [qcel.periodictable.to_E(z) for z in elem]

    covrad = np.fromiter((covalent_radii_lookup[caps[at]] for at in range(nat)), dtype=float, count=nat)
    return np.divide(covrad, qcel.constants.bohr2angstroms)


def _get_key(x, y, z, b):
    """Return key string from point values and block resolution"""

    return """{},{},{}""".format(x - x % b, y - y % b, z - z % b)


def _distance2(v, u):
    """Compute the square distance between points defined by vectors *v* and *u*."""

    return sum(((v[i] - u[i]) * (v[i] - u[i]) for i in range(len(v))))


def _get_blocks(geom, blocksize):
    """Parition atoms into spatial blocks"""

    allblocks = collections.defaultdict(list)
    for at in range(geom.shape[0]):
        x, y, z = (int(math.floor(geom[at][j])) for j in range(3))
        xyz_key = _get_key(x, y, z, blocksize)
        allblocks[xyz_key].append(at)
    return allblocks


def _get_bond_tree(radii, geom, allblocks, blocksize, bond_threshold):
    """Create bond tree from atomic coordinates"""

    bond_tree = [[] for at in range(geom.shape[0])]
    for blk in allblocks:
        atom_list = _get_atoms_from_blocks(_get_neighbor_blocks(blk, blocksize, allblocks), allblocks)
        for at1 in allblocks[blk]:
            for at2 in atom_list:
                r2_ij = _distance2(geom[at1], geom[at2])
                r2_thresh = bond_threshold * (radii[at1] + radii[at2])**2
                if at1 != at2 and r2_ij <= r2_thresh:
                    if at2 not in bond_tree[at1]:
                        bond_tree[at1].append(at2)
                    if at1 not in bond_tree[at2]:
                        bond_tree[at2].append(at1)
    return bond_tree


def _get_neighbor_blocks(block, blocksize, allblocks):
    """Find occupied blocks which neighbor `block`, including self"""

    x, y, z = (int(block.split(',')[j]) for j in range(3))
    neighbor_blocks = [_get_key(x + blocksize * (i - 1),
                                y + blocksize * (j - 1),
                                z + blocksize * (k - 1),
                                blocksize)
                       for i in range(3)
                       for j in range(3)
                       for k in range(3)]  # yapf: disable
    active_blocks = list(set(neighbor_blocks) & set(allblocks))
    return active_blocks


def _get_atoms_from_blocks(blocks, master_blocks):
    """Get list of atoms in a set of blocks"""

    atoms_nested = [master_blocks[blk] for blk in blocks]
    atoms = [at for sublist in atoms_nested for at in sublist]
    return atoms
