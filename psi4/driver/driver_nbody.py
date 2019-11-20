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

import copy
import itertools
import math
from typing import Any, Callable, Dict, List, Union
from ast import literal_eval

import pydantic

import pprint
pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)
import logging

import numpy as np
import pydantic
from qcelemental.models import DriverEnum, AtomicResult

from psi4 import core
from psi4.driver import constants, driver_nbody_multilevel, p4util, pp
from psi4.driver.p4util.exceptions import *
from psi4.driver.task_base import AtomicComputer, BaseComputer

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

### Math helper functions


def nCr(n, r):
    f = math.factorial
    return f(n) // f(r) // f(n - r)


def _sum_cluster_ptype_data(ptype,
                            ptype_dict,
                            compute_list,
                            fragment_slice_dict,
                            fragment_size_dict,
                            ret,
                            vmfc=False,
                            n=0,
                            level=1):
    """
    Sum gradients/hessians from n-body computations to obtain the BSSE corrected or uncorrected gradient/hessian

    Parameters
    ----------
    ptype : str
        Either "gradient" or "hessian"

    ptype_dict : dict
        Dictionary containing computed gradient or Hessian obtained from each subsystem computation

    compute_list : tuple
        A tuple of (frag, basis) data containing all the required computations

    fragment_slice_dict : dict
        Dictionary containing slices that index the gradient or Hessian matrix for each of the fragments

    fragment_size_dict : dict
        Dictionary containing the number of atoms of each fragment

    ret : np.ndarray
        An array containing the returned gradient or Hessian data. Modified in place

    vmfc : bool
        Is it a VMFC calculation

    n : int
        MBE level; required for VMFC calculations

    Returns
    -------
    None
    """


    if len(compute_list) == 0:
        return

    sign = 1

    # Do ptype
    if ptype == 'gradient':
        for fragn, basisn in compute_list:
            key = str(level) + "_" + str((fragn, basisn))
            start = 0
            grad = np.asarray(ptype_dict[key])

            if vmfc:
                sign = ((-1)**(n - len(fragn)))

            for bas in basisn:
                end = start + fragment_size_dict[bas]
                ret[fragment_slice_dict[bas]] += sign * grad[start:end]
                start += fragment_size_dict[bas]

    elif ptype == 'hessian':
        for fragn, basisn in compute_list:
            key = str(level) + "_" + str((fragn, basisn))
            #hess = np.asarray(ptype_dict[(fragn, basisn)])
            hess = np.asarray(ptype_dict[key])

            if vmfc:
                sign = ((-1)**(n - len(fragn)))

            # Build up start and end slices
            abs_start, rel_start = 0, 0
            abs_slices, rel_slices = [], []
            for bas in basisn:
                rel_end = rel_start + 3 * fragment_size_dict[bas]
                rel_slices.append(slice(rel_start, rel_end))
                rel_start += 3 * fragment_size_dict[bas]

                tmp_slice = fragment_slice_dict[bas]
                abs_slices.append(slice(tmp_slice.start * 3, tmp_slice.stop * 3))

            for abs_sl1, rel_sl1 in zip(abs_slices, rel_slices):
                for abs_sl2, rel_sl2 in zip(abs_slices, rel_slices):
                    ret[abs_sl1, abs_sl2] += sign * hess[rel_sl1, rel_sl2]

    else:
        raise KeyError("ptype can only be gradient or hessian How did you end up here?")


def _print_nbody_energy(energy_body_dict, header, embedding=False):
    core.print_out("""\n   ==> N-Body: %s energies <==\n\n""" % header)
    core.print_out("""   n-Body     Total Energy [Eh]       I.E. [kcal/mol]      Delta [kcal/mol]\n""")
    previous_e = energy_body_dict[1]
    if previous_e == 0.0:
        tot_e = False
    else:
        tot_e = True
    nbody_range = list(energy_body_dict)
    nbody_range.sort()
    for n in nbody_range:
        delta_e = (energy_body_dict[n] - previous_e)
        delta_e_kcal = delta_e * constants.hartree2kcalmol
        int_e_kcal = (
            energy_body_dict[n] - energy_body_dict[1]) * constants.hartree2kcalmol if not embedding else np.nan
        if tot_e:
            core.print_out("""     %4s  %20.12f  %20.12f  %20.12f\n""" % (n, energy_body_dict[n], int_e_kcal,
                                                                       delta_e_kcal))
        else:
            core.print_out("""     %4s  %20s  %20.12f  %20.12f\n""" % (n, "N/A", int_e_kcal,
                                                                       delta_e_kcal))
        previous_e = energy_body_dict[n]
    core.print_out("\n")


    """
    Computes the nbody interaction energy, gradient, or Hessian depending on input.
    This is a generalized universal function for computing interaction and total quantities.

    :returns: *return type of func* |w--w| The data.

    :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| data and wavefunction with energy/gradient/hessian set appropriately when **return_wfn** specified.

    :type func: Callable
    :param func: ``energy`` || etc.

        Python function that accepts method_string and a molecule. Returns a
        energy, gradient, or Hessian as requested.

    :type method_string: str
    :param method_string: ``'scf'`` || ``'mp2'`` || ``'ci5'`` || etc.

        First argument, lowercase and usually unlabeled. Indicates the computational
        method to be passed to func.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :type return_wfn: :ref:`boolean <op_py_boolean>`
    :param return_wfn: ``'on'`` || |dl| ``'off'`` |dr|

        Indicate to additionally return the :py:class:`~psi4.core.Wavefunction`
        calculation result as the second element of a tuple.

    :type bsse_type: str or list
    :param bsse_type: ``'cp'`` || ``['nocp', 'vmfc']`` || |dl| ``None`` |dr| || etc.

        Type of BSSE correction to compute: CP for counterpoise correction, NoCP
        for plain supramolecular interaction energy, or VMFC for Valiron-Mayer
        Function Counterpoise correction. If a list is provided, the first string in
        the list determines which interaction or total energies/gradients/Hessians are
        returned by this function. By default, this function is not called.

    :type max_nbody: int
    :param max_nbody: ``3`` || etc.

        Maximum n-body to compute, cannot exceed the number of fragments in the molecule.

    :type ptype: str
    :param ptype: ``'energy'`` || ``'gradient'`` || ``'hessian'``

        Type of the procedure passed in.

    :type return_total_data: :ref:`boolean <op_py_boolean>`
    :param return_total_data: ``'on'`` || |dl| ``'off'`` |dr|

        If True returns the total data (energy/gradient/Hessian) of the system,
        otherwise returns interaction data. Default is ``'off'`` for energies,
        ``'on'`` for gradients and Hessians. Note that the calculation of total
        counterpoise corrected energies implies the calculation of the energies of
        monomers in the monomer basis, hence specifying ``return_total_data = True``
        may carry out more computations than ``return_total_data = False``.

    :type levels: dict
    :param levels: ``{1: 'ccsd(t)', 2: 'mp2', 'supersystem': 'scf'}`` || ``{1: 2, 2: 'ccsd(t)', 3: 'mp2'}`` || etc

        Dictionary of different levels of theory for different levels of expansion
        Note that method_string is not used in this case. ``supersystem`` computes
        all higher order n-body effects up to the number of fragments.

    :type embedding_charges: dict
    :param embedding_charges: ``{1: [-0.834, 0.417, 0.417], ..}``

        Dictionary of atom-centered point charges. keys: 1-based index of fragment, values: list of charges for each fragment.

    :type charge_method: str
    :param charge_method: ``scf/6-31g`` || ``b3lyp/6-31g*`` || etc

        Method to compute point charges for monomers. Overridden by ``embedding_charges``
        if both are provided.

    :type charge_type: str
    :param charge_type: ``MULLIKEN_CHARGES`` || ``LOWDIN_CHARGES``

        Default is ``MULLIKEN_CHARGES``
    """


def build_nbody_compute_list(bsse_type_list, nbody, max_frag):
    """Generates the list of N-Body computations to be performed for a given BSSE type.

    Parameters
    ----------
    metadata : dict of str
        Dictionary containing N-body metadata.

        Required ``'key': value`` pairs:
        ``'bsse_type_list'``: list of str
            List of requested BSSE treatments.  Possible values include lowercase ``'cp'``, ``'nocp'``,
            and ``'vmfc'``.
        ``'max_nbody'``: int
            Maximum number of bodies to include in the N-Body treatment.
            Possible: `max_nbody` <= `max_frag`
            Default: `max_nbody` = `max_frag`
        ``'max_frag'``: int
            Number of distinct fragments comprising full molecular supersystem.

    Returns
    -------
    metadata : dict of str
        Dictionary containing N-body metadata.

        New ``'key': value`` pair:
        ``'compute_dict'`` : dict of str: dict
            Dictionary containing subdicts enumerating compute lists for each possible BSSE treatment.

            Contents:
            ``'all'``: dict of int: set
                Set containing full list of computations required
            ``'cp'``: dict of int: set
                Set containing list of computations required for CP procedure
            ``'nocp'``: dict of int: set
                Set containing list of computations required for non-CP procedure
            ``'vmfc_compute'``: dict of int: set
                Set containing list of computations required for VMFC procedure
            ``'vmfc_levels'``: dict of int: set
                Set containing list of levels required for VMFC procedure
    """
    # What levels do we need?
    fragment_range = range(1, max_frag + 1)

    # Need nbody and all lower-body in full basis
    cp_compute_list = {x: set() for x in nbody}
    nocp_compute_list = {x: set() for x in nbody}
    vmfc_compute_list = {x: set() for x in nbody}
    vmfc_level_list = {x: set() for x in nbody}  # Need to sum something slightly different

    # Verify proper passing of bsse_type_list
    bsse_type_remainder = set(bsse_type_list) - {'cp', 'nocp', 'vmfc'}
    if bsse_type_remainder:
        raise ValidationError("""Unrecognized BSSE type(s): %s
Possible values are 'cp', 'nocp', and 'vmfc'.""" % ', '.join(str(i) for i in bsse_type_remainder))

    # Build up compute sets
    if 'cp' in bsse_type_list:
        # Everything is in full n-mer basis
        basis_tuple = tuple(fragment_range)

        for nb in nbody:
            if nb == 1:
                # Add monomers in monomer basis
                for x in fragment_range:
                    cp_compute_list[1].add(((x, ), (x, )))
            else:
                for sublevel in range(1, nb + 1):
                    for x in itertools.combinations(fragment_range, sublevel):
                        if nbody == 1: break
                        cp_compute_list[nb].add((x, basis_tuple))

    if 'nocp' in bsse_type_list or metadata['return_total_data']:
        # Everything in monomer basis
        for nb in nbody:
            for sublevel in range(1, nb + 1):
                for x in itertools.combinations(fragment_range, sublevel):
                    nocp_compute_list[nb].add((x, x))

    if 'vmfc' in bsse_type_list:
        # Like a CP for all combinations of pairs or greater
        for nb in nbody:
            for cp_combos in itertools.combinations(fragment_range, nb):
                basis_tuple = tuple(cp_combos)
                for interior_nbody in range(1, nb + 1):
                    for x in itertools.combinations(cp_combos, interior_nbody):
                        combo_tuple = (x, basis_tuple)
                        vmfc_compute_list[nb].add(combo_tuple)
                        vmfc_level_list[len(basis_tuple)].add(combo_tuple)
    # Build a comprehensive compute_range
    compute_list = {x: set() for x in nbody}
    for nb in nbody:
        compute_list[nb] |= cp_compute_list[nb]
        compute_list[nb] |= nocp_compute_list[nb]
        compute_list[nb] |= vmfc_compute_list[nb]
        core.print_out("        Number of %d-body computations:     %d\n" % (nb, len(compute_list[nb])))

    compute_dict = {
        'all': compute_list,
        'cp': cp_compute_list,
        'nocp': nocp_compute_list,
        'vmfc_compute': vmfc_compute_list,
        'vmfc_levels': vmfc_level_list
    }
    return compute_dict


def assemble_nbody_components(metadata, component_results):
    """Assembles N-body components into interaction quantities according to requested BSSE procedure(s).

    Parameters
    -----------
    metadata : dict of str
        Dictionary of N-body metadata.

        Required ``'key': value`` pairs:
            ``'ptype'``: {'energy', 'gradient', 'hessian'}
                   Procedure which has generated the N-body components to be combined.
            ``'bsse_type_list'``: list of str
                   List of requested BSSE treatments.  Possible values include lowercase ``'cp'``, ``'nocp'``,
                   and ``'vmfc'``.
            ``'max_nbody'``: int
                   Maximum number of bodies to include in the N-Body treatment.
                   Possible: `max_nbody` <= `max_frag`
                   Default: `max_nbody` = `max_frag`
            ``'max_frag'``: int
                   Number of distinct fragments comprising full molecular supersystem.
            ``'energies_dict'``: dict of set: float64
                   Dictionary containing all energy components required for given N-body procedure.
            ``'ptype_dict'``: dict of set: float64 or dict of set: psi4.Matrix
                   Dictionary of returned quantities from calls of function `func` during N-body computations
            ``'compute_dict'``: dict of str: dict
                   Dictionary containing {int: set} subdicts enumerating compute lists for each possible
                   BSSE treatment.
            ``'kwargs'``: dict
                   Arbitrary keyword arguments.
    component_results : dict of str: dict
        Dictionary containing computed N-body components.

        Required ``'key': value`` pairs:
        ``'energies'``: dict of set: float64
               Dictionary containing all energy components required for given N-body procedure.
        ``'ptype'``: dict of set: float64 or dict of set: psi4.Matrix
               Dictionary of returned quantities from calls of function `func` during N-body computations
        ``'intermediates'``: dict of str: float64
               Dictionary of psivars for intermediate N-body computations to be set at the end of the
               N-body procedure.
    Returns
    -------
    results : dict of str
        Dictionary of all N-body results.

        Contents:
        ``'ret_energy'``: float64
            Interaction data requested.  If multiple BSSE types requested in `bsse_type_list`, the interaction data associated with the *first* BSSE
            type in the list is returned.
        ``'nbody_dict'``: dict of str: float64
            Dictionary of relevant N-body psivars to be set
        ``'energy_body_dict'``: dict of int: float64
            Dictionary of total energies at each N-body level, i.e., ``results['energy_body_dict'][2]`` is the sum of all 2-body total energies
            for the supersystem. May be empty if ``return_total_data`` is ``False``.
        ``'ptype_body_dict'``: dict or dict of int: array_like
            Empty dictionary if `ptype is ``'energy'``, or dictionary of total ptype
            arrays at each N-body level; i.e., ``results['ptype_body_dict'][2]``
            for `ptype` ``'gradient'``is the total 2-body gradient.
    """

    quiet = metadata['quiet']

    # which level(s) are we assembling?
    level = 0
    level_idx = {int(i.split('_')[0]) for i in component_results['energies'].keys()}

    if len(level_idx) != 1:
        raise ValidationError("Something's wrong")

    # get the range of nbodies for this level
    for idx in level_idx:
        level = idx
        nbody_range = metadata['nbody_list'][level - 1]
        if nbody_range[0] == 'supersystem':
            # range for supersystem sub-components
            nbody_range = metadata['nbody_list'][level]

    compute_dict = build_nbody_compute_list(metadata['bsse_type'], nbody_range, metadata['max_frag'])
    cp_compute_list = compute_dict['cp']
    cp_ptype_compute_list = {x: set() for x in nbody_range}
    nocp_compute_list = compute_dict['nocp']
    nocp_ptype_compute_list = {x: set() for x in nbody_range}
    vmfc_compute_list = compute_dict['vmfc_compute']
    vmfc_level_list = compute_dict['vmfc_levels']

    # Build size and slices dictionaries
    fragment_size_dict = {
        frag: metadata['molecule'].extract_subsets(frag).natom()
        for frag in range(1, metadata['max_frag'] + 1)
    }
    start = 0
    fragment_slice_dict = {}
    for k, v in fragment_size_dict.items():
        fragment_slice_dict[k] = slice(start, start + v)
        start += v

    molecule_total_atoms = sum(fragment_size_dict.values())

    # Final dictionaries
    cp_energy_by_level = {n: 0.0 for n in range(1, nbody_range[-1] + 1)}
    nocp_energy_by_level = {n: 0.0 for n in range(1, nbody_range[-1] + 1)}

    cp_energy_body_dict = {n: 0.0 for n in range(1, nbody_range[-1] + 1)}
    nocp_energy_body_dict = {n: 0.0 for n in range(1, nbody_range[-1] + 1)}
    vmfc_energy_body_dict = {n: 0.0 for n in range(1, nbody_range[-1] + 1)}

    # Build out ptype dictionaries if needed
    if metadata['driver'] != 'energy':
        if metadata['driver'] == 'gradient':
            arr_shape = (molecule_total_atoms, 3)
        elif metadata['driver'] == 'hessian':
            arr_shape = (molecule_total_atoms * 3, molecule_total_atoms * 3)
        else:
            raise KeyError("N-Body: ptype '%s' not recognized" % metadata['driver'])

        cp_ptype_by_level = {n: np.zeros(arr_shape) for n in range(1, nbody_range[-1] + 1)}
        nocp_ptype_by_level = {n: np.zeros(arr_shape) for n in range(1, nbody_range[-1] + 1)}
        vmfc_ptype_by_level = {n: np.zeros(arr_shape) for n in range(1, nbody_range[-1] + 1)}

        cp_ptype_body_dict = {n: np.zeros(arr_shape) for n in range(1, nbody_range[-1] + 1)}
        nocp_ptype_body_dict = {n: np.zeros(arr_shape) for n in range(1, nbody_range[-1] + 1)}
        vmfc_ptype_body_dict = {n: np.zeros(arr_shape) for n in range(1, nbody_range[-1] + 1)}
    else:
        cp_ptype_by_level, cp_ptype_body_dict = {}, {}
        nocp_ptype_by_level, nocp_ptype_body_dict = {}, {}
        vmfc_ptype_body_dict = {}
        cp_mon_in_mon_basis = 0.0

    # Sum up all of the levels
    nbody_dict = {}

    add_cp_mon = False
    cp_add_dict = set()
    nocp_add_dict = set()

    # build the special list for cp gradients/hessians
    if metadata['driver'] != 'energy':
        for n in nbody_range:
            for v in cp_compute_list[n]:
                if len(v[1]) != 1:
                    cp_ptype_compute_list[len(v[0])].add(v)
            for w in nocp_compute_list[n]:
                nocp_ptype_compute_list[len(w[0])].add(w)

    cp_mon_in_mon_basis = 0.0
    for n in nbody_range:

        # Energy
        # Extract energies for monomers in monomer basis for CP total data
        if n == 1:
            monomers_in_monomer_basis = [v for v in cp_compute_list[1] if len(v[1]) == 1]
            monomer_energies = 0.0
            monomer_energy_list = []
            for i in monomers_in_monomer_basis:
                key = str(level) + "_" + str(i)
                monomer_energy_list.append(component_results['energies'][key])
                monomer_energies += component_results['energies'][key]
                cp_compute_list[1].remove(i)
            add_cp_mon = True

        if 'cp' in metadata['bsse_type']:
            for v in cp_compute_list[n]:
                if v not in cp_add_dict:
                    energy = component_results['energies'][str(level) + "_" + str(v)]
                    cp_energy_by_level[len(v[0])] += energy
                    cp_add_dict.add(v)
        if 'nocp' in metadata['bsse_type']:
            for v in nocp_compute_list[n]:
                if v not in nocp_add_dict:
                    energy = component_results['energies'][str(level) + "_" + str(v)]
                    nocp_energy_by_level[len(v[0])] += energy
                    nocp_add_dict.add(v)

        # Special vmfc case
        if n > 1:
            vmfc_energy_body_dict[n] = vmfc_energy_body_dict[n - 1]
        for tup in vmfc_level_list[n]:
            key = str(level) + "_" + str(tup)
            vmfc_energy_body_dict[n] += ((-1)**(n - len(tup[0]))) * component_results['energies'][key]

        # Do ptype
        if metadata['driver'] != 'energy':
            _sum_cluster_ptype_data(metadata['driver'],
                                    component_results['ptype'],
                                    cp_ptype_compute_list[n],
                                    fragment_slice_dict,
                                    fragment_size_dict,
                                    cp_ptype_by_level[n],
                                    level=level)
            _sum_cluster_ptype_data(metadata['driver'],
                                    component_results['ptype'],
                                    nocp_ptype_compute_list[n],
                                    fragment_slice_dict,
                                    fragment_size_dict,
                                    nocp_ptype_by_level[n],
                                    level=level)
            _sum_cluster_ptype_data(metadata['driver'],
                                    component_results['ptype'],
                                    vmfc_level_list[n],
                                    fragment_slice_dict,
                                    fragment_size_dict,
                                    vmfc_ptype_by_level[n],
                                    vmfc=True,
                                    n=n,
                                    level=level)

    if metadata['driver'] != 'energy':
        # Extract ptype data for monomers in monomer basis for CP total data
        monomer_ptype = np.zeros(arr_shape)
        _sum_cluster_ptype_data(metadata['driver'], component_results['ptype'], monomers_in_monomer_basis,
                                fragment_slice_dict, fragment_size_dict, monomer_ptype)

    # Compute cp energy and ptype
    if 'cp' in metadata['bsse_type']:
        for n in range(1, nbody_range[-1] + 1):
            if n == metadata['max_frag']:
                cp_energy_body_dict[n] = cp_energy_by_level[n] - bsse
                if metadata['driver'] != 'energy':
                    cp_ptype_body_dict[n][:] = cp_ptype_by_level[n] - bsse_ptype
                continue

            for k in range(1, n + 1):
                take_nk = nCr(metadata['max_frag'] - k - 1, n - k)
                sign = ((-1)**(n - k))
                value = cp_energy_by_level[k]
                cp_energy_body_dict[n] += sign * value

                if metadata['driver'] != 'energy':
                    value = cp_ptype_by_level[k]
                    cp_ptype_body_dict[n] += sign * value

            if n == 1:
                bsse = cp_energy_body_dict[n] - monomer_energies
                cp_energy_body_dict[n] = monomer_energies
                if metadata['driver'] != 'energy':
                    bsse_ptype = cp_ptype_body_dict[n] - monomer_ptype
                    cp_ptype_body_dict[n] = monomer_ptype.copy()

            else:
                cp_energy_body_dict[n] -= bsse
                if metadata['driver'] != 'energy':
                    cp_ptype_body_dict[n] -= bsse_ptype

        cp_interaction_energy = cp_energy_body_dict[metadata['max_nbody']] - cp_energy_body_dict[1]
        nbody_dict['Counterpoise Corrected Interaction Energy'] = cp_interaction_energy

        for n in nbody_range[1:]:
            var_key = 'CP-CORRECTED %d-BODY INTERACTION ENERGY' % n
            nbody_dict[var_key] = cp_energy_body_dict[n] - cp_energy_body_dict[1]

        if not quiet:
            _print_nbody_energy(cp_energy_body_dict, "Counterpoise Corrected (CP)", metadata['embedding_charges'])
        cp_interaction_energy = cp_energy_body_dict[metadata['max_nbody']] - cp_energy_body_dict[1]
        if monomer_energies != 0.0:
            nbody_dict['Counterpoise Corrected Total Energy'] = cp_energy_body_dict[metadata['max_nbody']]
        nbody_dict['Counterpoise Corrected Interaction Energy'] = cp_interaction_energy

    # Compute nocp energy and ptype
    if 'nocp' in metadata['bsse_type']:
        for n in range(1, nbody_range[-1] + 1):
            if n == metadata['max_frag']:
                nocp_energy_body_dict[n] = nocp_energy_by_level[n]
                if metadata['driver'] != 'energy':
                    nocp_ptype_body_dict[n][:] = nocp_ptype_by_level[n]
                continue

            for k in range(1, n + 1):
                take_nk = nCr(metadata['max_frag'] - k - 1, n - k)
                sign = ((-1)**(n - k))
                value = nocp_energy_by_level[k]
                nocp_energy_body_dict[n] += take_nk * sign * value

                if metadata['driver'] != 'energy':
                    value = nocp_ptype_by_level[k]
                    nocp_ptype_body_dict[n] += take_nk * sign * value

        if not quiet:
            _print_nbody_energy(nocp_energy_body_dict, "Non-Counterpoise Corrected (NoCP)",
                                metadata['embedding_charges'])
        nocp_interaction_energy = nocp_energy_body_dict[metadata['max_nbody']] - nocp_energy_body_dict[1]
        nbody_dict['Non-Counterpoise Corrected Total Energy'] = nocp_energy_body_dict[metadata['max_nbody']]
        nbody_dict['Non-Counterpoise Corrected Interaction Energy'] = nocp_interaction_energy

        for n in nbody_range[1:]:
            var_key = 'NOCP-CORRECTED %d-BODY INTERACTION ENERGY' % n
            nbody_dict[var_key] = nocp_energy_body_dict[n] - nocp_energy_body_dict[1]

    # Compute vmfc ptype
    if 'vmfc' in metadata['bsse_type']:
        if metadata['driver'] != 'energy':
            for n in nbody_range:
                if n > 1:
                    vmfc_ptype_body_dict[n] = vmfc_ptype_by_level[n - 1]
                vmfc_ptype_body_dict[n] += vmfc_ptype_by_level[n]

        if not quiet:
            _print_nbody_energy(vmfc_energy_body_dict, "Valiron-Mayer Function Couterpoise (VMFC)",
                                metadata['embedding_charges'])
        vmfc_interaction_energy = vmfc_energy_body_dict[metadata['max_nbody']] - vmfc_energy_body_dict[1]
        nbody_dict['Valiron-Mayer Function Couterpoise Total Energy'] = vmfc_energy_body_dict[metadata['max_nbody']]
        nbody_dict['Valiron-Mayer Function Couterpoise Interaction Energy'] = vmfc_interaction_energy

        for n in nbody_range[1:]:
            var_key = 'VMFC-CORRECTED %d-BODY INTERACTION ENERGY' % n
            nbody_dict[var_key] = vmfc_energy_body_dict[n] - vmfc_energy_body_dict[1]

    # Returns
    results = {}
    results['nbody'] = nbody_dict
    for b in ['cp', 'nocp', 'vmfc']:
        if monomer_energies != 0.0:
            results['%s_energy_body_dict' % b] = locals('%s_energy_body_dict' % b)
            results['%s_energy_body_dict' % b] = {str(i) + b: j for i, j in results['%s_energy_body_dict' % b].items()}
        else:
            results['%s_energy_body_dict' % b] = {}

    # Figure out and build return types
    return_method = metadata['bsse_type'][0]

    if return_method == 'cp':
        results['ptype_body_dict'] = cp_ptype_body_dict
        results['energy_body_dict'] = cp_energy_body_dict
    elif return_method == 'nocp':
        results['ptype_body_dict'] = nocp_ptype_body_dict
        results['energy_body_dict'] = nocp_energy_body_dict
    elif return_method == 'vmfc':
        results['ptype_body_dict'] = vmfc_ptype_body_dict
        results['energy_body_dict'] = vmfc_energy_body_dict
    else:
        raise ValidationError(
            "N-Body Wrapper: Invalid return type. Should never be here, please post this error on github.")

    if metadata['return_total_data']:
        results['ret_energy'] = results['energy_body_dict'][metadata['max_nbody']]
    else:
        results['ret_energy'] = results['energy_body_dict'][metadata['max_nbody']]
        results['ret_energy'] -= results['energy_body_dict'][1]

    if metadata['driver'] != 'energy':
        if metadata['return_total_data']:
            np_final_ptype = results['ptype_body_dict'][metadata['max_nbody']].copy()
        else:
            np_final_ptype = results['ptype_body_dict'][metadata['max_nbody']].copy()
            np_final_ptype -= results['ptype_body_dict'][1]

        results['ret_ptype'] = np_final_ptype
    else:
        results['ret_ptype'] = results['ret_energy']

    if monomer_energies == 0.0:
        del results['energy_body_dict']

    return results


def electrostatic_embedding(embedding_charges):
    """
    Add atom-centered point charges
    """
    from psi4.driver import qmmm

    # Add embedding point charges
    Chrgfield = qmmm.QMMMbohr()
    for i in embedding_charges:
        Chrgfield.extern.addCharge(i[0], i[1][0], i[1][1], i[1][2])
    core.set_global_option_python('EXTERN', Chrgfield.extern)


class ManyBodyComputer(BaseComputer):

    molecule: Any
    driver: DriverEnum

    keywords: Dict[str, Any] = {}

    # nbody kwargs
    bsse_type: List[str] = ["cp"]
    max_nbody: int = -1
    max_frag: int = None
    nbody_list: List[int] = []
    nbody_level: int = 1
    return_wfn: bool = False
    return_total_data: bool
    embedding_charges: Dict[int, list] = {}
    quiet: bool = False

    task_list: Dict[str, AtomicComputer] = {}

    def __init__(self, **data):

        if not isinstance(data['bsse_type'], list):
            data['bsse_type'] = [data['bsse_type']]

        BaseComputer.__init__(self, **data)

        self.max_frag = self.molecule.nfragments()
        if self.max_nbody == -1:
            self.max_nbody = self.molecule.nfragments()
        else:
            self.max_nbody = min(self.max_nbody, self.max_frag)

        if not self.return_total_data:
            if self.embedding_charges:
                raise Exception('Cannot return interaction data when using embedding scheme.')

    @pydantic.validator('molecule')
    def set_molecule(cls, mol):
        mol.update_geometry()
        mol.fix_com(True)
        mol.fix_orientation(True)
        return mol

    @pydantic.validator('bsse_type')
    def set_bsse_type(cls, bsse_type):

        bsse_types = []
        for bt in bsse_type:
            if bt.lower() not in ['cp', 'nocp', 'vmfc']:
                raise ValidationError(f"N-Body GUFunc: bsse_type '{bt}' is not recognized")
            bsse_types.append(bt.lower())

        return bsse_types

    @pydantic.validator("return_total_data")
    def set_return_total_data(cls, v, values):
        if "return_total_data" in values:
            return v
        else:
            if values["driver"] in ["gradient", "hessian"]:
                return True
            else:
                return False

    # Build task for a given level
    def build_tasks(self, obj, **kwargs):
        bsse_type = "all"

        nbody_level = kwargs.pop('nlevel', self.max_nbody)
        level = kwargs.pop('level', None)
        template = copy.deepcopy(kwargs)

        # Store tasks by level
        # Get the n-body orders for this level
        nbody_list = self.nbody_list
        nbody = nbody_list[nbody_level]

        # Add supersystem computation if requested
        if level == 'supersystem':
            data = template
            data["molecule"] = self.molecule
            self.task_list[str(level) + '_' + str(self.max_frag)] = obj(**data)
            logger.debug(str(level) + ', ' + str(self.max_frag))
            compute_dict = build_nbody_compute_list(self.bsse_type, [i for i in range(1, self.max_nbody + 1)],
                                                    self.max_frag)
        else:
            # Get compute list
            compute_dict = build_nbody_compute_list(self.bsse_type, nbody, self.max_frag)

        compute_list = compute_dict[bsse_type]

        counter = 0

        # Add current compute list to the master task list
        for count, n in enumerate(compute_list):
            for key, pair in enumerate(compute_list[n]):
                lbl = str(nbody_level + 1) + '_' + str(pair)
                if lbl in self.task_list:
                    continue
                logger.debug(pair)
                data = template
                ghost = list(set(pair[1]) - set(pair[0]))
                data["molecule"] = self.molecule.extract_subsets(list(pair[0]), ghost)
                if self.embedding_charges:
                    embedding_frags = list(set(range(1, self.max_frag + 1)) - set(pair[1]))
                    charges = []
                    for frag in embedding_frags:
                        positions = self.molecule.extract_subsets(frag).geometry().np.tolist()
                        charges.extend([[chg, i] for i, chg in zip(positions, self.embedding_charges[frag])])
                    data['keywords']['function_kwargs'].update({'embedding_charges': charges})

                self.task_list[lbl] = obj(**data)
                counter += 1

        return counter

    def plan(self):
        return [t.plan() for t in self.task_list.values()]

    def compute(self, client=None):
        with p4util.hold_options_state():
            for k, t in self.task_list.items():
                t.compute(client)

        if self.embedding_charges:
            core.set_global_option_python('EXTERN', None)

    def prepare_results(self, results={}, client=None):

        nlevels = len({i.split('_')[0] for i in self.task_list})
        if nlevels > 1 and not results:
            return driver_nbody_multilevel.prepare_results(self, client)

        metadata = self.dict().copy()
        results_list = {k: v.get_results(client=client) for k, v in (results.items() or self.task_list.items())}
        energies = {k: v.properties.return_energy for k, v in results_list.items()}

        ptype = None
        if self.driver.name in ['gradient', 'hessian']:
            ptype = {k: v.return_result for k, v in results_list.items()}
        tmp = {"energies": energies, "ptype": ptype}
        nbody_results = assemble_nbody_components(metadata.copy(), tmp)

        if self.driver == 'hessian':
            # Compute nbody gradient
            gradient = {k: v.extras['qcvars']["CURRENT GRADIENT"] for k, v in results_list.items()}
            tmp['ptype'] = gradient
            metadata = metadata.copy()
            metadata['driver'] = 'gradient'
            grad_result = assemble_nbody_components(metadata, tmp)
            nbody_results.update({'ret_gradient': grad_result['ret_ptype']})

        nbody_results['intermediates'] = {
            "N-BODY (%s)@(%s) TOTAL ENERGY" % (', '.join(str(i) for i in literal_eval(k.split("_")[1])[0]), ', '.join(
                str(i) for i in literal_eval(k.split("_")[1])[1])): v.properties.return_energy
            for k, v in results_list.items()
        }

        nbody_results['intermediates2'] = {str(k): v.properties.return_energy for k, v in results_list.items()}

        if ptype is not None:
            nbody_results['intermediates_ptype'] = {str(k): v for k, v in ptype.items()}

        return nbody_results

    def get_results(self, client=None):
        """Return nbody results as nbody-flavored QCSchema."""

        results = self.prepare_results(client=client)
        ret_energy, ret_ptype, ret_gradient = results.pop('ret_energy'), results.pop('ret_ptype'), results.pop(
            'ret_gradient', None)

        # load QCVariables
        qcvars = {
            'NUCLEAR REPULSION ENERGY': self.molecule.nuclear_repulsion_energy(),
            'NBODY NUMBER': len(self.task_list),
        }

        for k, val in results.items():
            qcvars[k] = val

        qcvars['CURRENT ENERGY'] = ret_energy
        if self.driver == 'gradient':
            qcvars['CURRENT GRADIENT'] = ret_ptype
        elif self.driver == 'hessian':
            qcvars['CURRENT GRADIENT'] = ret_gradient
            qcvars['CURRENT HESSIAN'] = ret_ptype

        component_results = self.dict()['task_list']
        for k, val in component_results.items():
            val['molecule'] = val['molecule'].to_schema(dtype=2)

        nbodyjob = AtomicResult(
            **{
                'driver': self.driver,
                'model': {
                    'method': self.method,
                    'basis': self.basis,
                },
                'molecule': self.molecule.to_schema(dtype=2),
                'properties': {
                    'calcinfo_natom': self.molecule.natom(),
                    'nuclear_repulsion_energy': self.molecule.nuclear_repulsion_energy(),
                    'return_energy': ret_energy,
                },
                'provenance': p4util.provenance_stamp(__name__),
                'extras': {
                    'qcvars': qcvars,
                    'component_results': component_results,
                },
                'return_result': ret_ptype,
                'success': True,
            })

        return nbodyjob

    def get_psi_results(self, return_wfn=False):

        nbody_results = self.get_results()
        ret = nbody_results.return_result

        wfn = core.Wavefunction.build(self.molecule, 'def2-svp')

        # this should be all the setting that's needed (and shouldn't need the isinstance
        #   to avoid dicts). can this be simplified, Asim?
        for qcv, val in nbody_results.extras['qcvars'].items():
            if isinstance(val, (int, float)):
                for obj in [core, wfn]:
                    obj.set_variable(qcv, val)

        dicts = [
            'energies', 'ptype', 'intermediates', 'intermediates2', 'intermediates_ptype', 'energy_body_dict',
            'gradient_body_dict', 'nbody', 'cp_energy_body_dict', 'nocp_energy_body_dict', 'vmfc_energy_body_dict'
        ]
        if self.driver == 'gradient':
            ret = core.Matrix.from_array(ret)
            wfn.set_gradient(ret)
            nbody_results.extras['qcvars']['gradient_body_dict'] = nbody_results.extras['qcvars']['ptype_body_dict']
        elif self.driver == 'hessian':
            ret = core.Matrix.from_array(ret)
            grad = core.Matrix.from_array(nbody_results.extras['qcvars']['CURRENT GRADIENT'])
            wfn.set_hessian(ret)
            wfn.set_gradient(grad)

        for d in nbody_results.extras['qcvars']:
            if d in dicts:
                for var, value in nbody_results.extras['qcvars'][d].items():
                    try:
                        wfn.set_variable(str(var), value)
                        core.set_variable(str(var), value)
                    except:
                        wfn.set_variable(self.driver + ' ' + str(var), value)

        core.set_variable("CURRENT ENERGY", nbody_results.properties.return_energy)
        # wfn.set_energy(nbody_results.properties['return_energy'])  # catches CURRENT ENERGY on Wfn

        if return_wfn:
            return (ret, wfn)
        else:
            return ret
