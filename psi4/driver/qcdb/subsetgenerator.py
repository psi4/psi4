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

"""Module containing functions that, when passed a qcdb.WrappedDatabase instance
*dbinstance*, return an array of reaction names that are a subset of
dbinstance.hrxn.keys(). Since the full database is passed in, reactions
can be filtered by molecule characteristics, reaction names, existing
subsets, etc. The official name of the subset is specified by the function
docstring. Second line of docstring becomes tagl.

"""


def genset_MXuDD(dbinstance):
    """mxdd
    near-equilibrium systems also in mxdd

    """
    try:
        ssA = set(dbinstance.sset['mx'].keys())
    except KeyError:
        ssA = set()
    try:
        ssB = set(dbinstance.sset['dd'].keys())
    except KeyError:
        ssB = set()
    return ssA.union(ssB)


def genset_HBn5min(dbinstance):
    """HB-5min
    near-equilibrium systems also in hb

    """
    try:
        ssA = set(dbinstance.sset['hb'].keys())
    except KeyError:
        ssA = set()
    try:
        ssB = set(dbinstance.sset['5min'].keys())
    except KeyError:
        ssB = set()
    return ssA.intersection(ssB)


def genset_MXn5min(dbinstance):
    """MX-5min
    near-equilibrium systems also in mx

    """
    try:
        ssA = set(dbinstance.sset['mx'].keys())
    except KeyError:
        ssA = set()
    try:
        ssB = set(dbinstance.sset['5min'].keys())
    except KeyError:
        ssB = set()
    return ssA.intersection(ssB)


def genset_DDn5min(dbinstance):
    """DD-5min
    near-equilibrium systems also in dd

    """
    try:
        ssA = set(dbinstance.sset['dd'].keys())
    except KeyError:
        ssA = set()
    try:
        ssB = set(dbinstance.sset['5min'].keys())
    except KeyError:
        ssB = set()
    return ssA.intersection(ssB)


def genset_MXDDPPn5min(dbinstance):
    """MXDDPP-5min
    near-equilibrium systems also in mxddpp

    """
    try:
        ssA = set(dbinstance.sset['mxddpp'].keys())
    except KeyError:
        ssA = set()
    try:
        ssB = set(dbinstance.sset['5min'].keys())
    except KeyError:
        ssB = set()
    return ssA.intersection(ssB)


def genset_MXDDNPn5min(dbinstance):
    """MXDDNP-5min
    near-equilibrium systems also in mxddnp

    """
    try:
        ssA = set(dbinstance.sset['mxddnp'].keys())
    except KeyError:
        ssA = set()
    try:
        ssB = set(dbinstance.sset['5min'].keys())
    except KeyError:
        ssB = set()
    return ssA.intersection(ssB)


def genset_allneutral(dbinstance):
    """neutral
    systems where all components are neutral

    """
    eligible = []
    for rxn, orxn in dbinstance.hrxn.items():
        if all([True if rgt.charge == 0 else False for rgt in orxn.rxnm['default'].keys()]):
            eligible.append(rxn)
    return eligible


def genset_anyanion(dbinstance):
    """anion
    systems where any component is an anion

    """
    eligible = []
    for rxn, orxn in dbinstance.hrxn.items():
        for rgt in orxn.rxnm['default'].keys():
            if rgt.charge < 0:
                eligible.append(rxn)
                break
    return eligible


def genset_anycation(dbinstance):
    """cation
    systems where any component is a cation

    """
    eligible = []
    for rxn, orxn in dbinstance.hrxn.items():
        for rgt in orxn.rxnm['default'].keys():
            if rgt.charge > 0:
                eligible.append(rxn)
                break
    return eligible
