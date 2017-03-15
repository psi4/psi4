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

"""Module with utility functions that act on molecule objects."""
from __future__ import absolute_import
import math
import re

from psi4 import core
from psi4.driver.p4util import constants, filter_comments
from psi4.driver.inputparser import process_pubchem_command, pubchemre


def extract_clusters(mol, ghost=True, cluster_size=0):
    """Function to return all subclusters of the molecule *mol* of
    real size *cluster_size* and all other atoms ghosted if *ghost*
    equals true, all other atoms discarded if *ghost* is false. If
    *cluster_size* = 0, returns all possible combinations of cluster size.

    """
    # How many levels of clusters are possible?
    nfrag = mol.nfragments()

    # Initialize the cluster array
    clusters = []

    # scope the arrays
    reals = []
    ghosts = []

    # counter
    counter = 0

    # loop over all possible cluster sizes
    for nreal in range(nfrag, 0, -1):

        # if a specific cluster size size is requested, only do that
        if (nreal != cluster_size and cluster_size > 0):
            continue

        # initialize the reals list
        reals = []

        # setup first combination [3,2,1] lexical ordering
        # fragments indexing is 1's based, bloody hell
        for index in range(nreal, 0, -1):
            reals.append(index)

        # start loop through lexical promotion
        while True:

            counter = counter + 1

            # Generate cluster from last iteration
            if (ghost):
                ghosts = []
                for g in range(nfrag, 0, -1):
                    if (g not in reals):
                        ghosts.append(g)
                clusters.append(mol.extract_subsets(reals, ghosts))
            else:
                clusters.append(mol.extract_subsets(reals))

            # reset rank
            rank = 0

            # look for lexical promotion opportunity
            # i.e.: [4 2 1] has a promotion opportunity at
            #   index 1 to produce [4 3 1]
            for k in range(nreal - 2, -1, -1):
                if (reals[k] != reals[k + 1] + 1):
                    rank = k + 1
                    break

            # do the promotion
            reals[rank] = reals[rank] + 1

            # demote the right portion of the register
            val = 1
            for k in range(nreal - 1, rank, -1):
                reals[k] = val
                val = val + 1

            # boundary condition is promotion into
            # [nfrag+1 nfrag-1 ...]
            if (reals[0] > nfrag):
                break

    return clusters


def extract_cluster_indexing(mol, cluster_size=0):
    """Function to returns a LIST of all subclusters of the molecule *mol* of
    real size *cluster_size*. If *cluster_size* = 0, returns all possible
    combinations of cluster size.

    """
    import copy

    # How many levels of clusters are possible?
    nfrag = mol.nfragments()

    # Initialize the cluster array
    clusters = []

    # scope the arrays
    reals = []

    # counter
    counter = 0

    # loop over all possible cluster sizes
    for nreal in range(nfrag, 0, -1):

        # if a specific cluster size size is requested, only do that
        if (nreal != cluster_size and cluster_size > 0):
            continue

        # initialize the reals list
        reals = []

        # setup first combination [3,2,1] lexical ordering
        # fragments indexing is 1's based, bloody hell
        for index in range(nreal, 0, -1):
            reals.append(index)

        # start loop through lexical promotion
        while True:

            counter = counter + 1

            # Generate cluster from last iteration
            clusters.append(copy.deepcopy(reals))

            # reset rank
            rank = 0

            # look for lexical promotion opportunity
            # i.e.: [4 2 1] has a promotion opportunity at
            #   index 1 to produce [4 3 1]
            for k in range(nreal - 2, -1, -1):
                if (reals[k] != reals[k + 1] + 1):
                    rank = k + 1
                    break

            # do the promotion
            reals[rank] = reals[rank] + 1

            # demote the right portion of the register
            val = 1
            for k in range(nreal - 1, rank, -1):
                reals[k] = val
                val = val + 1

            # boundary condition is promotion into
            # [nfrag+1 nfrag-1 ...]
            if (reals[0] > nfrag):
                break

    return clusters


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


def BFS(self):
    """Perform a breadth-first search (BFS) on the real atoms
    in molecule, returning an array of atom indices of fragments.
    Relies upon van der Waals radii and so faulty for close
    (esp. hydrogen-bonded) fragments. Original code from
    Michael S. Marshall.

    """
    vdW_diameter = {
        'H':  1.001 / 1.5,
        'HE': 1.012 / 1.5,
        'LI': 0.825 / 1.5,
        'BE': 1.408 / 1.5,
        'B':  1.485 / 1.5,
        'C':  1.452 / 1.5,
        'N':  1.397 / 1.5,
        'O':  1.342 / 1.5,
        'F':  1.287 / 1.5,
        'NE': 1.243 / 1.5,
        'NA': 1.144 / 1.5,
        'MG': 1.364 / 1.5,
        'AL': 1.639 / 1.5,
        'SI': 1.716 / 1.5,
        'P':  1.705 / 1.5,
        'S':  1.683 / 1.5,
        'CL': 1.639 / 1.5,
        'AR': 1.595 / 1.5}

    Queue = []
    White = range(self.natom())  # untouched
    Black = []  # touched and all edges discovered
    Fragment = []  # stores fragments

    start = 0  # starts with the first atom in the list
    Queue.append(start)
    White.remove(start)

    # Simply start with the first atom, do a BFS when done, go to any
    #   untouched atom and start again iterate until all atoms belong
    #   to a fragment group
    while White or Queue:                # Iterates to the next fragment
        Fragment.append([])

        while Queue:                     # BFS within a fragment
            for u in Queue:              # find all white neighbors to vertex u
                for i in White:
                    dist = constants.bohr2angstroms * math.sqrt(
                           (self.x(i) - self.x(u)) ** 2 +
                           (self.y(i) - self.y(u)) ** 2 +
                           (self.z(i) - self.z(u)) ** 2)
                    if dist < vdW_diameter[self.symbol(u)] + \
                       vdW_diameter[self.symbol(i)]:
                        Queue.append(i)  # if you find you, put in the queue
                        White.remove(i)  # & remove it from the untouched list
            Queue.remove(u)              # remove focus from Queue
            Black.append(u)
            Fragment[-1].append(int(u))  # add to group (0-indexed)
            Fragment[-1].sort()          # preserve original atom ordering

        if White:                        # can't move White -> Queue if empty
            Queue.append(White[0])
            White.remove(White[0])

    return Fragment


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

    activate(molecule)

    return molecule


def activate(mol):
    """Function to set molecule object *mol* as the current active molecule.
    Permitted for user use but deprecated in driver in favor of explicit
    molecule-passing.

    """
    core.set_active_molecule(mol)
