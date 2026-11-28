#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2026 The Psi4 Developers.
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
"""
AO operator matrices of the one-electron properties of SAPT property calculations.

The interaction-induced SAPT properties are linear in the property operator, so every
property that is linear in the one-particle density is evaluated the same way: one
AO matrix per component, whose trace with a density matrix is the electronic part of
that component. The components are named and shaped as :py:class:`psi4.core.OEProp`
names and shapes them.
"""

import re
from dataclasses import dataclass
from typing import Callable, List, Tuple, Union

import numpy as np

from psi4 import core

from ...constants import constants
from ...p4util import text
from ...p4util.exceptions import ValidationError
from ...p4util.python_helpers import _multipole_order, _multipole_plumper

__all__ = ["Component", "Quantity", "validate_properties", "build_operators"]

# Properties linear in the one-particle density, with few components
_OPERATOR_PROPERTIES = ["DIPOLE", "QUADRUPOLE", "ESP_AT_NUCLEI", "MULLIKEN_CHARGES", "LOWDIN_CHARGES"]

# Linear in the density, but with one component per grid point
_GRID_PROPERTIES = ["GRID_ESP", "GRID_FIELD"]

_MULTIPOLE_RE = re.compile(r"^MULTIPOLES?\((\d+)\)$")

# OEProp names multipoles up to this order
_MAX_MULTIPOLE_ORDER = len(_multipole_order) - 1


@dataclass
class Component:
    """One property component: the trace of `matrix` with a density matrix is its
    electronic part, and `nuclear` gives its nuclear part from the nuclear charges
    of the atoms."""

    matrix: np.ndarray
    nuclear: Callable[[np.ndarray], float]


@dataclass
class Quantity:
    """A property as OEProp reports it, e.g. "DIPOLE" or "ESP AT CENTER 2": the
    `components` it is made of, how to `assemble` their values into its value, and
    `labels` of the components for printing."""

    name: str
    components: List[int]
    assemble: Callable[[np.ndarray], Union[float, np.ndarray]]
    labels: List[str]


def validate_properties(properties: List[str]) -> List[str]:
    """Upper-cased `properties`, or a ValidationError if any of them is not a property
    linear in the one-particle density that SAPT properties support."""

    if not len(properties):
        raise ValidationError("SAPT properties: No properties specified!")

    # imported here, as the procrouting modules import each other
    from ..response.scf_response import property_dicts
    response_properties = list(property_dicts)

    supported = []
    for prop in properties:
        prop = prop.upper()

        if prop in _OPERATOR_PROPERTIES:
            supported.append(prop)
            continue

        match = _MULTIPOLE_RE.match(prop)
        if match:
            order = int(match.group(1))
            if not (1 <= order <= _MAX_MULTIPOLE_ORDER):
                raise ValidationError(
                    f"SAPT properties: {prop} must have an order between 1 and {_MAX_MULTIPOLE_ORDER}.")
            supported.append(prop)
            continue

        if prop in _GRID_PROPERTIES:
            raise ValidationError(f"SAPT properties: {prop} is not available yet. It needs one operator per grid "
                                  "point, so it will come with the interaction-induced density matrices.")

        if prop in core.OEProp.valid_methods or prop in response_properties:
            raise ValidationError(f"SAPT properties: {prop} is not linear in the one-particle density, so it has "
                                  "no interaction-induced SAPT correction. Supported properties: "
                                  f"{', '.join(_OPERATOR_PROPERTIES)} and MULTIPOLE(n).")

        alternatives = text.find_approximate_string_matches(prop, _OPERATOR_PROPERTIES, 2)
        alternatives = f" Did you mean? {' '.join(alternatives)}" if alternatives else ""
        raise ValidationError(f"SAPT properties: Feature '{prop}' is not recognized.{alternatives}")

    return supported


def properties_origin(molecule: core.Molecule) -> List[float]:
    """The origin [a0] of the multipoles, from the PROPERTIES_ORIGIN option, as
    :py:class:`psi4.core.OEProp` takes it."""

    if not core.has_global_option_changed("PROPERTIES_ORIGIN"):
        return [0.0, 0.0, 0.0]

    spec = core.get_global_option("PROPERTIES_ORIGIN")
    if len(spec) == 1:
        center = str(spec[0]).upper()
        if center == "COM":
            weights = np.array([molecule.mass(atom) for atom in range(molecule.natom())])
        elif center == "NUCLEAR_CHARGE":
            weights = np.array([molecule.charge(atom) for atom in range(molecule.natom())])
        else:
            raise ValidationError("Invalid specification of PROPERTIES_ORIGIN.  Please consult the manual.")
        geometry = np.asarray(molecule.geometry())
        return list(weights @ geometry / weights.sum())

    if len(spec) == 3:
        origin = np.array(spec, dtype=float)
        if molecule.units() == "Angstrom":
            origin /= constants.bohr2angstroms
        return list(origin)

    raise ValidationError("Invalid specification of PROPERTIES_ORIGIN.  Please consult the manual.")


def _cartesian_powers(order: int) -> List[Tuple[int, int, int]]:
    """(lx, ly, lz) of the components of a multipole of `order`, in the order of
    :py:meth:`psi4.core.MintsHelper.ao_multipoles`."""

    powers = []
    for ii in range(order + 1):
        lx = order - ii
        for lz in range(ii + 1):
            powers.append((lx, ii - lz, lz))
    return powers


def _multipole_name(order: int) -> str:
    """OEProp's name of a multipole of `order`."""

    return "DIPOLE" if order == 1 else _multipole_order[order]


def _multipole_orders(properties: List[str]) -> List[int]:
    """The multipole orders the properties need, as OEProp computes them: MULTIPOLE(n)
    gives every order from 1 to n."""

    orders = set()
    for prop in properties:
        if prop == "DIPOLE":
            orders.add(1)
        elif prop == "QUADRUPOLE":
            orders.add(2)
        else:
            match = _MULTIPOLE_RE.match(prop)
            if match:
                orders.update(range(1, int(match.group(1)) + 1))
    return sorted(orders)


def build_operators(properties: List[str], basisset: core.BasisSet,
                    molecule: core.Molecule) -> Tuple[List[Component], List[Quantity]]:
    """
    The AO operators of the properties.

    Parameters
    ----------
    properties
        Property names, as returned by :py:func:`validate_properties`.
    basisset
        The dimer-centered basis set the operators are built in.
    molecule
        The dimer, whose atoms the per-atom properties and the multipole origin refer to.
        Its atoms are in the order of the atoms of `basisset`.

    Returns
    -------
    components, quantities
        Every operator component once, and the properties made of them.
    """

    mints = core.MintsHelper(basisset)
    geometry = np.asarray(molecule.geometry())
    natom = molecule.natom()

    components = []
    quantities = []

    # Multipoles about the origin of PROPERTIES_ORIGIN
    orders = _multipole_orders(properties)
    if orders:
        origin = properties_origin(molecule)
        relative = geometry - np.array(origin)
        matrices = mints.ao_multipoles(max(orders), origin)

        offset = 0
        for order in range(1, max(orders) + 1):
            powers = _cartesian_powers(order)
            if order in orders:
                indices = []
                for num, (lx, ly, lz) in enumerate(powers):
                    nuclear_moments = relative[:, 0]**lx * relative[:, 1]**ly * relative[:, 2]**lz
                    components.append(
                        Component(np.asarray(matrices[offset + num]),
                                  lambda charges, moments=nuclear_moments: float(charges @ moments)))
                    indices.append(len(components) - 1)

                if order == 1:
                    assemble = np.array
                else:
                    assemble = lambda values, order=order: _multipole_plumper(np.asarray(values), order)
                labels = ["X" * lx + "Y" * ly + "Z" * lz for lx, ly, lz in powers]
                quantities.append(Quantity(_multipole_name(order), indices, assemble, labels))
            offset += len(powers)

    # Electrostatic potential at the nuclei. The potential integrals are those of a
    # positive charge, so the electronic operator takes a minus sign. The nuclear part
    # leaves out the nucleus the potential is taken at.
    if "ESP_AT_NUCLEI" in properties:
        for atom in range(natom):
            potential = np.asarray(mints.ao_multipole_potential(0, list(geometry[atom]))[0])
            distances = np.linalg.norm(geometry - geometry[atom], axis=1)
            inverse = np.zeros(natom)
            others = np.arange(natom) != atom
            inverse[others] = 1.0 / distances[others]

            components.append(Component(-potential, lambda charges, inverse=inverse: float(charges @ inverse)))
            quantities.append(
                Quantity(f"ESP AT CENTER {atom + 1}", [len(components) - 1], lambda values: float(values[0]), [""]))

    # Population charges: the nuclear charge minus the population of the basis
    # functions of the atom
    charge_properties = [prop for prop in ("MULLIKEN_CHARGES", "LOWDIN_CHARGES") if prop in properties]
    if charge_properties:
        overlap = np.asarray(mints.ao_overlap())
        centers = np.array([basisset.function_to_center(func) for func in range(basisset.nbf())])

        if "LOWDIN_CHARGES" in charge_properties:
            evals, evecs = np.linalg.eigh(overlap)
            overlap_half = (evecs * evals**0.5) @ evecs.T

        for prop in charge_properties:
            indices = []
            for atom in range(natom):
                projector = np.diag((centers == atom).astype(float))
                if prop == "MULLIKEN_CHARGES":
                    matrix = -0.5 * (overlap @ projector + projector @ overlap)
                else:
                    matrix = -overlap_half @ projector @ overlap_half

                components.append(Component(matrix, lambda charges, atom=atom: float(charges[atom])))
                indices.append(len(components) - 1)

            labels = [f"{molecule.symbol(atom)}{atom + 1}" for atom in range(natom)]
            quantities.append(Quantity(prop.replace("_", " "), indices, np.array, labels))

    return components, quantities
