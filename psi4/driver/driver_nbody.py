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
"""Plan, run, and assemble QC tasks to obtain many-body expansion and basis-set superposition error treatments.

=============
ManyBody Flow
=============
Bullet points are major actions
Lines of dashes denote function calls
e/d/dd=dg/g/h := energy, dipole, dipole derivative = dipole gradient, gradient, Hessian
`mc_(frag, bas)` := a modelchem index, mc; indices of real fragments, frag; set(bas - frag) are indices of ghost fragments. see "intermediates_energy" in big table below for example.
note that there's a lot of natural 1-indexing (1, 2, 3) rather than 0-indexing (0, 1, 2) in manybody. e.g., 2-body energy, Molecule.extract_subsets(1, (1, 2))
note that a "level" can be n-body level (how many real molecular fragments) or a modelchem level (`mc_`; e.g., CC on 1-bodies, MP2 on 2-bodies; "multilevel")

---------------------------
ManyBodyComputer.__init__()
---------------------------
* not an explicit function but pydantic handles some defaults and validation
* fields molecule, nfragments, bsse_type, return_total_data, and initial max_nbody set
* BaseComputer.__init__()

task_planner.py::task_planner()
-------------------------------
* computer gets modified from task_planner outside this file!
* modelchem (method and basis) treatment levels for each n-body level determined from user levels kwarg. fields nbodies_per_mc_level set and max_nbody reset
* for each modelchem treatment level, call build_tasks() below via one of four routes, depending on simple MB or layered MB(FD), MB(CBS), or MB(FD(CBS))

    ------------------------------
    ManyBodyComputer.build_tasks()
    ------------------------------
    * if supersystem requested as a modelchem level, request (frag, bas) indices for full nbody range of nocp treatment from build_nbody_compute_list()
    * otherwise, request (frag, bas) indices for specified nbody range covering specified bsse treatments from build_nbody_compute_list()

        build_nbody_compute_list()
        --------------------------
        * initializes dicts for each of nocp, cp, vmfc (2 for this one) with keys requested n-body levels and values empty sets
        * use combinatorics formulas to fill each key with (frag, bas) indices (what fragments are active and what fragments have basis functions)
          needed to compute the requested bsse treatments at the requested n-body levels.
        * merge by n-body level the sets of indices for each bsse treatment into an "all" dict. return this and all the per-bsse dicts.

    * merge (from different bsse_types) all the requested indices, prepend a modelchem treatment index to form `mc_(frag, bas)`
    * construct a molecule appropriately real/ghosted from active-fragment info in (frag, bas)
    * if embedding_charges active, prepare external_potentials array for atoms not in bas fragments
    * for any new `mc_(frag, bas)` index, append a new computer to self.task_list

--------------------------
ManyBodyComputer.compute()
--------------------------
* compute() for each job in self.task_list

----------------------------------
ManyBodyComputer.get_psi_results()
----------------------------------

    Computer.get_results()
    ----------------------

        Computer._prepare_results()
        ---------------------------
        * if multiple modelchems (multilevel):

            multilevel.prepare_results()
            ----------------------------
            * from the pool of calcs, partition them by modelchem treatment and call _prepare_results on each subpool
            * sums modelchem levels and returns small dict back to get_results()

        * call get_results() for each job in task list
        * assemble all the computed energies, all the computed gradients, and all the computed hessians
        * for each available derivative, call:

            assemble_nbody_components()
            ---------------------------
            * re-call build_nbody_compute_list to get the cp/nocp/vmfc lists again

                build_nbody_compute_list()
                --------------------------

            * slice up the supersystem mol into fragment atom ranges to help g/h arrays build piecemeal
            * prepare empty {bsse_type}_by_level and {bsse_type}_body_dict structs. the former have different contents for vmfc
            * for cp and nocp, resort the build_nbody_compute_list returns into per-body lists suitable for summing
            * note that nb loops often run over more than active nbodies_per_mc_level item due to 1-body for subtraction and multilevel complications
            * for each possibly active n-body level and each active bsse_type, call _sum_cluster_ptype_data to build by_level structs

                _sum_cluster_ptype_data()
                -------------------------
                * sum up ene, grad, or Hess in per-fragment pieces based on list of (frag, bas) subjobs active for that bsse treatment

            * compute special case of monomers in monomer basis
            * for each of cp/nocp/vmfc, apply appropriate formula to build each n-body level of cumulative total energy into body_dict
            * for driver=energy, set several qcvars and call:

                _print_nbody_energy()
                ---------------------
                * prints and logs formatted energy output. called separately for cp, nocp, vmfc

            * collect qcvars and summed levels into a return dictionary with some extra aliases for target bsse_type and target driver

        * merge all the assemble_nbody_components return dictionaries
        * in struct["intermediates"], store dict of `"N-BODY (?)@(?) TOTAL ENERGY" = return_energy` for all in task_list or results kwarg
        * in struct["intermediates_{ptype}"], store dict of `task_list key = return_{ptype}` for all in task_list or results kwarg. ptype=e/g/h
          always for ptype=energy, as available for higher derivatives when driver=g/h

    * form nbody qcvars and properties, inc'l number, current e/g/h as available
    * pull results (incl dicts!) into qcvars
    * form model, including copy of class with mols converted to qcsk at atomicresult.extras["component_results"]

* collect ManyBody-flavored AtomicResult from self.get_results()
* build wfn from nbody mol and basis (always def2-svp)
* push qcvars to P::e and wfn. push various internal dicts to qcvars, too
* convert result to psi4.core.Matrix (non-energy) and set g/h on wfn
* return e/g/h and wfn

"""

__all__ = [
    "BsseEnum",
    "ManyBodyComputer",
    "nbody",
]

import copy
import itertools
import math
from typing import Any, Callable, Dict, List, Literal, Optional, Sequence, Set, Tuple, Union, TYPE_CHECKING
from ast import literal_eval
from enum import Enum
try:
    from pydantic.v1 import Field, validator
except ImportError:
    from pydantic import Field, validator

import logging

import numpy as np
from qcelemental.models import DriverEnum, AtomicResult

from psi4 import core
from .constants import constants, pp
from . import driver_nbody_multilevel, p4util
from .p4util.exceptions import *
from .task_base import BaseComputer, AtomicComputer, EnergyGradientHessianWfnReturn
from .driver_cbs import CompositeComputer
from .driver_findif import FiniteDifferenceComputer

if TYPE_CHECKING:
    import qcportal

logger = logging.getLogger(__name__)

FragBasIndex = Tuple[Tuple[int], Tuple[int]]

SubTaskComputers = Union[AtomicComputer, CompositeComputer, FiniteDifferenceComputer]

def nbody():
    """
    Computes the nbody interaction energy, gradient, or Hessian depending on input.
    This is a generalized universal function for computing interaction and total quantities.

    :returns: *return type of func* |w--w| The data.

    :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| data and wavefunction with energy/gradient/hessian set appropriately when **return_wfn** specified.

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
        returned by this function. By default, many-body treatments are inactive.

    :type max_nbody: int
    :param max_nbody: ``3`` || etc.

        Maximum n-body to compute, cannot exceed the number of fragments in the molecule.

    :type return_total_data: :ref:`boolean <op_py_boolean>`
    :param return_total_data: ``'on'`` || |dl| ``'off'`` |dr|

        If True returns the total data (energy/gradient/Hessian) of the system,
        otherwise returns interaction data. Default is ``'off'`` for energies,
        ``'on'`` for gradients and Hessians. Note that the calculation of total
        counterpoise corrected energies implies the calculation of the energies of
        monomers in the monomer basis, hence specifying ``return_total_data = True``
        may carry out more computations than ``return_total_data = False``.
        For gradients and Hessians, ``return_total_data = False`` is rarely useful.

    :type levels: dict
    :param levels: ``{1: 'ccsd(t)', 2: 'mp2', 'supersystem': 'scf'}`` || ``{1: 2, 2: 'ccsd(t)', 3: 'mp2'}`` || etc

        Dictionary of different levels of theory for different levels of expansion
        Note that method_string is not used in this case. ``supersystem`` computes
        all higher order n-body effects up to the number of fragments.

    :type embedding_charges: dict
    :param embedding_charges: ``{1: [-0.834, 0.417, 0.417], ..}``

        Dictionary of atom-centered point charges. keys: 1-based index of fragment, values: list of charges for each fragment.
        Add atom-centered point charges for fragments whose basis sets are not included in the computation.

    """
    pass


class BsseEnum(str, Enum):
    """Available basis-set superposition error (BSSE) treatments."""

    nocp = "nocp"  # plain supramolecular interaction energy
    cp = "cp"      # counterpoise correction
    vmfc = "vmfc"  # Valiron-Mayer function counterpoise


def _sum_cluster_ptype_data(
    ptype: DriverEnum,
    ptype_dict: Dict,
    compute_list: Set[FragBasIndex],
    fragment_slice_dict: Dict[int, Sequence],
    fragment_size_dict: Dict[int, int],
    mc_level_lbl: int,
    vmfc: bool = False,
    nb: int = 0,
) -> Union[float, np.ndarray]:
    """
    Sum arrays from n-body computations to obtain the BSSE corrected or uncorrected scalar or array.

    Parameters
    ----------
    ptype
        Hint to shape of array data to sum.
    ptype_dict
        Dictionary containing computed energy, gradient, or Hessian obtained from each subsystem computation
    compute_list
        A list of (frag, bas) tuples notating all the required computations.
    fragment_slice_dict
        Dictionary containing slices that index the gradient or Hessian matrix for each of the 1-indexed fragments.
        For He--HOOH--Me cluster, `{1: slice(0, 1, None), 2: slice(1, 5, None), 3: slice(5, 10, None)}`.
    fragment_size_dict
        Dictionary containing the number of atoms of each 1-indexed fragment.
        For He--HOOH--Me cluster, `{1: 1, 2: 4, 3: 5}`.
    vmfc
        Is it a VMFC calculation?
    nb
        n-body level; required for VMFC calculations.
    mc_level_lbl
        User label for what modelchem level results should be pulled out of *ptype_dict*.
        This is the 1-indexed counterpart to 0-indexed mc_level_idx.

    Returns
    -------
    ret
        Scalar or array containing the summed energy, gradient, or Hessian result.
        Formerly, passed in and modified in place and only called for g/h.

    """
    sign = 1
    nat = sum(fragment_size_dict.values())

    def labeler(frag: Tuple, bas:Tuple) -> str:
        return str(mc_level_lbl) + "_" + str((frag, bas))

    if ptype == "energy":
        ret = 0.0

        for frag, bas in compute_list:
            ene = ptype_dict[labeler(frag, bas)]

            if vmfc:
                sign = ((-1)**(nb - len(frag)))

            ret += sign * ene

        return ret

    elif ptype == 'gradient':
        ret = np.zeros((nat, 3))

        for frag, bas in compute_list:
            grad = np.asarray(ptype_dict[labeler(frag, bas)])

            if vmfc:
                sign = ((-1)**(nb - len(frag)))

            start = 0
            for ifr in bas:
                end = start + fragment_size_dict[ifr]
                ret[fragment_slice_dict[ifr]] += sign * grad[start:end]
                start += fragment_size_dict[ifr]

        return ret

    elif ptype == 'hessian':
        ret = np.zeros((nat * 3, nat * 3))

        for frag, bas in compute_list:
            hess = np.asarray(ptype_dict[labeler(frag, bas)])

            if vmfc:
                sign = ((-1)**(nb - len(frag)))

            # Build up start and end slices
            abs_start, rel_start = 0, 0
            abs_slices, rel_slices = [], []
            for ifr in bas:
                rel_end = rel_start + 3 * fragment_size_dict[ifr]
                rel_slices.append(slice(rel_start, rel_end))
                rel_start += 3 * fragment_size_dict[ifr]

                tmp_slice = fragment_slice_dict[ifr]
                abs_slices.append(slice(tmp_slice.start * 3, tmp_slice.stop * 3))

            for abs_sl1, rel_sl1 in zip(abs_slices, rel_slices):
                for abs_sl2, rel_sl2 in zip(abs_slices, rel_slices):
                    ret[abs_sl1, abs_sl2] += sign * hess[rel_sl1, rel_sl2]

        return ret

    else:
        raise KeyError("ptype can only be energy, gradient, or hessian. How did you end up here?")


def _print_nbody_energy(energy_body_dict: Dict[int, float], header: str, nfragments: int, embedding: bool = False):
    """Format output string for user for a single bsse_type. Prints to output and logger.
    Called repeatedly by assemble_nbody_component."""

    info = f"""\n   ==> N-Body: {header} energies <==\n\n"""
    info += f"""  {"n-Body":>12}     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy\n"""
    info += f"""                   [Eh]                    [Eh]                  [kcal/mol]            [Eh]                  [kcal/mol]\n"""
    previous_e = energy_body_dict[1]
    tot_e = (previous_e != 0.0)
    nbody_range = list(energy_body_dict)
    nbody_range.sort()
    for nb in range(1, nfragments + 1):
        lbl = []
        if nb == nfragments:
            lbl.append("FULL")
        if nb == max(nbody_range):
            lbl.append("RTN")
        lbl = "/".join(lbl)

        if nb in nbody_range:
            delta_e = (energy_body_dict[nb] - previous_e)
            delta_e_kcal = delta_e * constants.hartree2kcalmol
            if embedding:
                int_e = np.nan
                int_e_kcal = np.nan
            else:
                int_e = energy_body_dict[nb] - energy_body_dict[1]
                int_e_kcal = int_e * constants.hartree2kcalmol
            if tot_e:
                info += f"""  {lbl:>8} {nb:3}  {energy_body_dict[nb]:20.12f}  {int_e:20.12f}  {int_e_kcal:20.12f}  {delta_e:20.12f}  {delta_e_kcal:20.12f}\n"""
            else:
                info += f"""  {lbl:>8} {nb:3}  {"N/A":20}  {int_e:20.12f}  {int_e_kcal:20.12f}  {delta_e:20.12f}  {delta_e_kcal:20.12f}\n"""
            previous_e = energy_body_dict[nb]
        else:
            info += f"""  {lbl:>8} {nb:3}        {"N/A":20}  {"N/A":20}  {"N/A":20}  {"N/A":20}  {"N/A":20}\n"""

    info += "\n"
    core.print_out(info)
    logger.info(info)


def build_nbody_compute_list(
    bsse_type: List[BsseEnum],
    nbodies: List[Union[int, Literal["supersystem"]]],
    nfragments: int,
    return_total_data: bool,
    verbose: int = 1,
) -> Dict[str, Dict[int, Set[FragBasIndex]]]:
    """Generates lists of N-Body computations needed for requested BSSE treatments.

    Parameters
    ----------
    bsse_type
        Requested BSSE treatments.
    nbodies
        List of n-body levels (e.g., `[2]` or `[1, 2]` or `["supersystem"]`) for which to generate tasks.
        Often this value is an element of self.nbodies_per_mc_level.
        Note the natural 1-indexing, so `[1]` covers one-body contributions.
        Formerly nbody
    nfragments
        Number of distinct fragments comprising the full molecular supersystem. Usually self.nfragments.
        Formerly max_frag
    return_total_data
        Whether the total data (True; energy/gradient/Hessian) of the molecular system has been requested, as opposed to interaction data (False).
    verbose
        Control volume of printing.

    Returns
    -------
    compute_dict
        Dictionary containing subdicts enumerating compute lists for each possible BSSE treatment.
        Subdict keys are n-body levels and values are sets of all the `mc_(frag, bas)` indices
        needed to compute that n-body level. A given index can appear multiple times within a
        subdict and among subdicts.
        Formerly, the subdict values were sets of indices needed for given BSSE treatment _of_ given
        n-body level. See current (left) and former (right) definitions below for CP dimer.

            compute_dict["cp"] = {                  compute_dict["cp"] = {
                1: set(),                               1: {((1,), (1, 2)),
                2: {((1,), (1, 2)),                         ((2,), (1, 2))},
                    ((2,), (1, 2)),                     2: {((1, 2), (1, 2))}
                    ((1, 2), (1, 2))}               }
            }

        Subdicts below are always returned. Any may be empty if not requested through *bsse_type*.

        * ``'all'`` |w---w| full list of computations required
        * ``'cp'`` |w---w| list of computations required for CP procedure
        * ``'nocp'`` |w---w| list of computations required for non-CP procedure
        * ``'vmfc_compute'`` |w---w| list of computations required for VMFC procedure
        * ``'vmfc_levels'`` |w---w| list of levels required for VMFC procedure

    """
    # What levels do we need?
    fragment_range = range(1, nfragments + 1)

    # Need nbodies and all lower-body in full basis
    cp_compute_list = {x: set() for x in nbodies}
    nocp_compute_list = {x: set() for x in nbodies}
    vmfc_compute_list = {x: set() for x in nbodies}
    vmfc_level_list = {x: set() for x in nbodies}  # Need to sum something slightly different

    # Verify proper passing of bsse_type. already validated in Computer
    bsse_type_remainder = set(bsse_type) - {e.value for e in BsseEnum}
    if bsse_type_remainder:
        raise ValidationError("""Unrecognized BSSE type(s): {bsse_type_remainder}""")

    # Build up compute sets
    if 'cp' in bsse_type:
        # Everything is in full n-mer basis
        basis_tuple = tuple(fragment_range)

        for nb in nbodies:
            if nb > 1:
                for sublevel in range(1, nb + 1):
                    for x in itertools.combinations(fragment_range, sublevel):
                        # below was `nbodies`, which would never hit. present is closest to pre-DDD. purpose unclear to me.
                        # if self.max_nbody == 1: break
                        cp_compute_list[nb].add((x, basis_tuple))

    if 'nocp' in bsse_type or return_total_data:
        # Everything in monomer basis
        for nb in nbodies:
            for sublevel in range(1, nb + 1):
                for x in itertools.combinations(fragment_range, sublevel):
                    nocp_compute_list[nb].add((x, x))

    if 'vmfc' in bsse_type:
        # Like a CP for all combinations of pairs or greater
        for nb in nbodies:
            for cp_combos in itertools.combinations(fragment_range, nb):
                basis_tuple = tuple(cp_combos)
                for interior_nbody in range(1, nb + 1):
                    for x in itertools.combinations(cp_combos, interior_nbody):
                        combo_tuple = (x, basis_tuple)
                        vmfc_compute_list[nb].add(combo_tuple)
                        vmfc_level_list[len(basis_tuple)].add(combo_tuple)

    # Build a comprehensive compute range
    # * do not use list length to count number of {nb}-body computations
    compute_list = {x: set() for x in nbodies}
    for nb in nbodies:
        compute_list[nb] |= cp_compute_list[nb]
        compute_list[nb] |= nocp_compute_list[nb]
        compute_list[nb] |= vmfc_compute_list[nb]

    # Rearrange compute_list from key nb having values to compute all of that nb
    #   to key nb including values of that nb. Use for counting.
    compute_list_count = {x: set() for x in nbodies}
    for nb in nbodies:
        for nbset in compute_list.values():
            for item in nbset:
                if len(item[0]) == nb:
                    compute_list_count[nb].add(item)
    if verbose >= 1:
        info = "\n".join([f"        Number of {nb}-body computations:     {len(compute_list_count[nb])}" for nb in nbodies])
        core.print_out(info + "\n")
        logger.info(info)

    compute_dict = {
        'all': compute_list,
        'cp': cp_compute_list,
        'nocp': nocp_compute_list,
        'vmfc_compute': vmfc_compute_list,
        'vmfc_levels': vmfc_level_list
    }
    return compute_dict


def assemble_nbody_components(
    ptype: DriverEnum,
    component_results: Dict[str, Union[float, np.ndarray]],
    metadata: Dict[str, Any],
) -> Dict[str, Any]:
    """Assembles N-body components for a single derivative level and a single model chemistry level into interaction quantities according to requested BSSE treatment(s).

    Parameters
    ----------
    ptype
        Derivative level of component results to assemble. Matters mostly for scalar vs. array and array dimensions.
    component_results
        Dictionary with keys "mc_(frag, bas)" and values e/g/H computed component results according to *ptype*.
    metadata
        Dictionary of N-body metadata. Items described below.
        Later, assemble_nbody_components should become a class function and the below are simply class member data.
    quiet : bool
        See class field. Whether to print/log energy summaries. Default True. False used by multilevel to suppress per-mc-level printing.
    nfragments : int
        See class field. Number of distinct fragments comprising the full molecular supersystem.
        Formerly max_frag
    return_total_data : bool
        See class field. Whether the total data (e/g/H) of the molecular system has been requested, as opposed to interaction data.
    max_nbody : int
        See class field. Maximum number of bodies to include in the many-body treatment."
    embedding_charges : bool
        Whether embedding charges are present. Used to NaN the output printing rather than print bad numbers.
    molecule : psi4.core.Molecule
        See class field. Used to count atoms in fragments.
    nbodies_per_mc_level: List[List[Union[int, Literal["supersystem"]]]]
        See class field. Distribution of active n-body levels among model chemistry levels.
        Formerly nbody_list
    bsse_type : List[BsseEnum]
        See class field. Requested BSSE treatments. First in list determines which interaction or total energy/gradient/Hessian returned.
        Note that this is the only arg that gets RESET. Happens for supersystem "nbody".

    Returns
    -------
    results
        Dictionary of all N-body results. See contents at ManyBodyComputer.prepare_results docstring.

    """
    # which level are we assembling?
    mc_level_labels = {int(i.split("_")[0]) for i in component_results.keys()}

    if len(mc_level_labels) != 1:
        raise ValidationError(f"Something's wrong - this fn handles single-level (e.g., 1- & 2-body w/mp2) not multi-level (e.g., 1-body w/hf & 2-body w/mp2) assembly: len({mc_level_labels}) != 1")

    # get the range of nbodies for this level
    # * modelchem level label (mc_level_lbl) used in qcvars and dict keys is 1-indexed counterpart to 0-indexed modelchem level position (mc_level_idx) used to navigate self.nbodies_per_mc_level
    mc_level_lbl = list(mc_level_labels)[0]
    nbodies = metadata['nbodies_per_mc_level'][mc_level_lbl - 1]
    if nbodies[0] == 'supersystem':
        # range for supersystem sub-components
        nbodies = metadata['nbodies_per_mc_level'][mc_level_lbl]
        metadata['bsse_type'] = ['nocp']

    # regenerate per-bsse required calcs list
    compute_dict = build_nbody_compute_list(
        metadata['bsse_type'], nbodies, metadata['nfragments'], metadata["return_total_data"], verbose=0
    )

    # Build size and slices dictionaries
    fragment_size_dict = {}
    fragment_slice_dict = {}
    iat = 0
    for ifr in range(1, metadata["nfragments"] + 1):
        nat = metadata["molecule"].extract_subsets(ifr).natom()
        fragment_size_dict[ifr] = nat
        fragment_slice_dict[ifr] = slice(iat, iat + nat)
        iat += nat

    def shaped_zero(der: DriverEnum):
        if der == "energy":
            return 0.0
        elif der == "gradient":
            arr_shape = (nat, 3)
            return np.zeros(arr_shape)
        elif der == 'hessian':
            arr_shape = (nat * 3, nat * 3)
            return np.zeros(arr_shape)

    # Final dictionaries
    if ptype == "energy":
        cp_by_level = {n: 0.0 for n in range(1, nbodies[-1] + 1)}
        nocp_by_level = {n: 0.0 for n in range(1, nbodies[-1] + 1)}
        vmfc_by_level = {n: 0.0 for n in range(1, nbodies[-1] + 1)}

        cp_body_dict = {n: 0.0 for n in range(1, nbodies[-1] + 1)}
        nocp_body_dict = {n: 0.0 for n in range(1, nbodies[-1] + 1)}
        vmfc_body_dict = {n: 0.0 for n in range(1, nbodies[-1] + 1)}

    else:
        nat = sum(fragment_size_dict.values())
        if ptype == 'gradient':
            arr_shape = (nat, 3)
        elif ptype == 'hessian':
            arr_shape = (nat * 3, nat * 3)

        cp_by_level = {n: np.zeros(arr_shape) for n in range(1, nbodies[-1] + 1)}
        nocp_by_level = {n: np.zeros(arr_shape) for n in range(1, nbodies[-1] + 1)}
        vmfc_by_level = {n: np.zeros(arr_shape) for n in range(1, nbodies[-1] + 1)}

        cp_body_dict = {n: np.zeros(arr_shape) for n in range(1, nbodies[-1] + 1)}
        nocp_body_dict = {n: np.zeros(arr_shape) for n in range(1, nbodies[-1] + 1)}
        vmfc_body_dict = {n: np.zeros(arr_shape) for n in range(1, nbodies[-1] + 1)}

    # Sum up all of the levels
    # * compute_dict[bt][nb] holds all the computations needed to compute nb
    #   *not* all the nb-level computations, so build the latter
    cp_compute_list = {nb: set() for nb in range(1, nbodies[-1] + 1)}
    nocp_compute_list = {nb: set() for nb in range(1, nbodies[-1] + 1)}

    for nb in nbodies:
        for v in compute_dict["cp"][nb]:
            if len(v[1]) != 1:
                cp_compute_list[len(v[0])].add(v)
        for w in compute_dict["nocp"][nb]:
            nocp_compute_list[len(w[0])].add(w)

    for nb in range(1, nbodies[-1] + 1):
        cp_by_level[nb] = _sum_cluster_ptype_data(
            ptype,
            component_results,
            cp_compute_list[nb],
            fragment_slice_dict,
            fragment_size_dict,
            mc_level_lbl=mc_level_lbl,
        )
        nocp_by_level[nb] = _sum_cluster_ptype_data(
            ptype,
            component_results,
            nocp_compute_list[nb],
            fragment_slice_dict,
            fragment_size_dict,
            mc_level_lbl=mc_level_lbl,
        )
        if nb in compute_dict["vmfc_levels"]:
            vmfc_by_level[nb] = _sum_cluster_ptype_data(
                ptype,
                component_results,
                compute_dict["vmfc_levels"][nb],
                fragment_slice_dict,
                fragment_size_dict,
                vmfc=True,
                nb=nb,
                mc_level_lbl=mc_level_lbl,
            )

    def labeler(item) -> str:
        return str(mc_level_lbl) + "_" + str(item)

    # Extract data for monomers in monomer basis for CP total data
    if 1 in nbodies:
        monomers_in_monomer_basis = [v for v in compute_dict["nocp"][1] if len(v[1]) == 1]

        if ptype == "energy":
            monomer_energy_list = [component_results[labeler(m)] for m in monomers_in_monomer_basis]
            monomer_sum = sum(monomer_energy_list)
        else:
            monomer_sum = _sum_cluster_ptype_data(
                ptype,
                component_results,
                monomers_in_monomer_basis,
                fragment_slice_dict,
                fragment_size_dict,
                mc_level_lbl=mc_level_lbl,
            )
    else:
        monomer_sum = shaped_zero(ptype)

    nbody_dict = {}

    # Compute cp
    if 'cp' in metadata['bsse_type']:
        for nb in range(1, nbodies[-1] + 1):
            if nb == metadata['nfragments']:
                if ptype == "energy":
                    cp_body_dict[nb] = cp_by_level[nb] - bsse
                else:
                    cp_body_dict[nb][:] = cp_by_level[nb] - bsse
                continue

            for k in range(1, nb + 1):
                take_nk = math.comb(metadata['nfragments'] - k - 1, nb - k)
                sign = ((-1)**(nb - k))
                cp_body_dict[nb] += take_nk * sign * cp_by_level[k]

            if nb == 1:
                bsse = cp_body_dict[nb] - monomer_sum
                if ptype == "energy":
                    cp_body_dict[nb] = monomer_sum
                else:
                    cp_body_dict[nb] = monomer_sum.copy()
            else:
                cp_body_dict[nb] -= bsse

        if ptype == "energy":
            if not metadata["quiet"]:
                _print_nbody_energy(cp_body_dict, "Counterpoise Corrected (CP)", metadata["nfragments"], metadata['embedding_charges'])

            if monomer_sum != 0.0:
                nbody_dict["CP-CORRECTED TOTAL ENERGY"] = cp_body_dict[metadata['max_nbody']]
            nbody_dict["CP-CORRECTED INTERACTION ENERGY"] = cp_body_dict[metadata['max_nbody']] - cp_body_dict[1]

            for nb in nbodies[1:]:
                nbody_dict[f"CP-CORRECTED INTERACTION ENERGY THROUGH {nb}-BODY"] = cp_body_dict[nb] - cp_body_dict[1]
                nbody_dict[f"CP-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY"] = cp_body_dict[nb] - cp_body_dict[nb-1]
            for nb in nbodies:
                nbody_dict[f"CP-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY"] = cp_body_dict[nb]

    # Compute nocp
    if 'nocp' in metadata['bsse_type']:
        for nb in range(1, nbodies[-1] + 1):
            if nb == metadata['nfragments']:
                if ptype == "energy":
                    nocp_body_dict[nb] = nocp_by_level[nb]
                else:
                    nocp_body_dict[nb][:] = nocp_by_level[nb]
                continue

            for k in range(1, nb + 1):
                take_nk = math.comb(metadata['nfragments'] - k - 1, nb - k)
                sign = ((-1)**(nb - k))
                nocp_body_dict[nb] += take_nk * sign * nocp_by_level[k]

        if ptype == "energy":
            if not metadata["quiet"]:
                _print_nbody_energy(nocp_body_dict, "Non-Counterpoise Corrected (NoCP)", metadata["nfragments"], metadata['embedding_charges'])

            nbody_dict['NOCP-CORRECTED TOTAL ENERGY'] = nocp_body_dict[metadata['max_nbody']]
            nbody_dict['NOCP-CORRECTED INTERACTION ENERGY'] = nocp_body_dict[metadata['max_nbody']] - nocp_body_dict[1]

            for nb in nbodies[1:]:
                nbody_dict[f"NOCP-CORRECTED INTERACTION ENERGY THROUGH {nb}-BODY"] = nocp_body_dict[nb] - nocp_body_dict[1]
                nbody_dict[f"NOCP-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY"] = nocp_body_dict[nb] - nocp_body_dict[nb-1]
            for nb in nbodies:
                nbody_dict[f"NOCP-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY"] = nocp_body_dict[nb]

    # Compute vmfc
    if 'vmfc' in metadata['bsse_type']:
        for nb in nbodies:
            if ptype == "energy":
                for k in range(1, nb + 1):
                    vmfc_body_dict[nb] += vmfc_by_level[k]

            else:
                if nb > 1:
                    vmfc_body_dict[nb] = vmfc_by_level[nb - 1]
                vmfc_body_dict[nb] += vmfc_by_level[nb]

        if ptype == "energy":
            if not metadata["quiet"]:
                _print_nbody_energy(vmfc_body_dict, "Valiron-Mayer Function Counterpoise (VMFC)", metadata["nfragments"], metadata['embedding_charges'])

            vmfc_interaction_energy = vmfc_body_dict[metadata['max_nbody']] - vmfc_body_dict[1]
            nbody_dict['VMFC-CORRECTED TOTAL ENERGY'] = vmfc_body_dict[metadata['max_nbody']]
            nbody_dict['VMFC-CORRECTED INTERACTION ENERGY'] = vmfc_interaction_energy

            for nb in nbodies[1:]:
                nbody_dict[f"VMFC-CORRECTED INTERACTION ENERGY THROUGH {nb}-BODY"] = vmfc_body_dict[nb] - vmfc_body_dict[1]
                nbody_dict[f"VMFC-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY"] = vmfc_body_dict[nb] - vmfc_body_dict[nb-1]
            for nb in nbodies:
                nbody_dict[f"VMFC-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY"] = vmfc_body_dict[nb]

    # Collect specific and generalized returns
    results = {
        f"cp_{ptype}_body_dict" : {f"{nb}cp": j for nb, j in cp_body_dict.items()},
        f"nocp_{ptype}_body_dict": {f"{nb}nocp": j for nb, j in nocp_body_dict.items()},
        f"vmfc_{ptype}_body_dict": {f"{nb}vmfc": j for nb, j in vmfc_body_dict.items()},
    }

    if ptype == "energy":
        results['nbody'] = nbody_dict

    return_bsse_type = metadata["bsse_type"][0]

    if return_bsse_type == "cp":
        results[f"{ptype}_body_dict"] = cp_body_dict
    elif return_bsse_type == "nocp":
        results[f"{ptype}_body_dict"] = nocp_body_dict
    elif return_bsse_type == "vmfc":
        results[f"{ptype}_body_dict"] = vmfc_body_dict
    else:
        raise ValidationError(
            "N-Body Wrapper: Invalid return type. Should never be here, please post this error on github.")

    if ptype == "energy":
        piece = results[f"{ptype}_body_dict"][metadata['max_nbody']]
    else:
        piece = results[f"{ptype}_body_dict"][metadata['max_nbody']].copy()

    if metadata['return_total_data']:
        results[f"ret_{ptype}"] = piece
    else:
        results[f"ret_{ptype}"] = piece
        results[f"ret_{ptype}"] -= results[f"{ptype}_body_dict"][1]

    results['ret_ptype'] = results[f"ret_{ptype}"]

    return results


class ManyBodyComputer(BaseComputer):
    # user kwargs (all but levels become fields)
    # ------------------------------------------
    # * bsse_type
    # * levels
    # * max_nbody
    # * molecule  -- general
    # * return_total_data
    # * return_wfn -- general

    # fields set in construction
    # --------------------------
    # * nfragments (<- max_frag) -- from molecule

    # fields set in task_planner
    # --------------------------
    # * max_nbody -- from levels
    # * nbodies_per_mc_level -- from levels

    # TODO perhaps rework levels kwarg so that it's processed in class init into nbodies_per_mc_level. Right now, levels resets max_nbody.
    # TODO also, perhaps change nbodies_per_mc_level into dict of lists so that pos'n/label indexing coincides

    molecule: Any = Field(..., description="The target molecule, if not the last molecule defined.")
    basis: str = "(auto)"
    method: str = "(auto)"
    driver: DriverEnum = Field(..., description="The computation driver; i.e., energy, gradient, hessian.")
    keywords: Dict[str, Any] = Field({}, description="The computation keywords/options.")

    bsse_type: List[BsseEnum] = Field([BsseEnum.cp], description="Requested BSSE treatments. First in list determines which interaction or total energy/gradient/Hessian returned.")
    nfragments: int = Field(-1, description="Number of distinct fragments comprising full molecular supersystem.")  # formerly max_frag
    max_nbody: int = Field(-1, description="Maximum number of bodies to include in the many-body treatment. Possible: max_nbody <= nfragments. Default: max_nbody = nfragments.")

    nbodies_per_mc_level: List[List[Union[int, Literal["supersystem"]]]] = Field([], description="Distribution of active n-body levels among model chemistry levels. All bodies in range [1, self.max_nbody] must be present exactly once. Number of items in outer list is how many different modelchems. Each inner list specifies what n-bodies to be run at the corresponding modelchem (e.g., `[[1, 2]]` has max_nbody=2 and 1-body and 2-body contributions computed at the same level of theory; `[[1], [2]]` has max_nbody=2 and 1-body and 2-body contributions computed at different levels of theory. An entry 'supersystem' means all higher order n-body effects up to the number of fragments. The n-body levels are effectively sorted in the outer list, and any 'supersystem' element is at the end.")  # formerly nbody_list

    embedding_charges: Dict[int, List[float]] = Field({}, description="Atom-centered point charges to be used on molecule fragments whose basis sets are not included in the computation. Keys: 1-based index of fragment. Values: list of atom charges for that fragment.")

    return_total_data: Optional[bool] = Field(None, description="When True, returns the total data (energy/gradient/Hessian) of the system, otherwise returns interaction data. Default is False for energies, True for gradients and Hessians. Note that the calculation of total counterpoise corrected energies implies the calculation of the energies of monomers in the monomer basis, hence specifying ``return_total_data = True`` may carry out more computations than ``return_total_data = False``.")
    quiet: bool = Field(False, description="Whether to print/log formatted n-body energy analysis. Presently used by multi to suppress output. Candidate for removal from class once in-class/out-of-class functions sorted.")

    task_list: Dict[str, SubTaskComputers] = {}

    # Note that validation of user fields happens through typing and validator functions, so no class __init__ needed.

    @validator("bsse_type", pre=True)
    def set_bsse_type(cls, v):
        if not isinstance(v, list):
            v = [v]
        # emulate ordered set
        return list(dict.fromkeys([bt.lower() for bt in v]))

    @validator('molecule')
    def set_molecule(cls, mol):
        mol.update_geometry()
        mol.fix_com(True)
        mol.fix_orientation(True)
        return mol

    @validator("nfragments", always=True)
    def set_nfragments(cls, v, values):
        return values["molecule"].nfragments()

    @validator("max_nbody", always=True)
    def set_max_nbody(cls, v, values):
        if v == -1:
            return values["nfragments"]
        else:
            return min(v, values["nfragments"])

    @validator("embedding_charges")
    def set_embedding_charges(cls, v, values):
        if len(v) != values["nfragments"]:
            raise ValueError("embedding_charges dict should have entries for each 1-indexed fragment.")

        return v

    @validator("return_total_data", always=True)
    def set_return_total_data(cls, v, values):
        if v is not None:
            rtd = v
        elif values["driver"] in ["gradient", "hessian"]:
            rtd = True
        else:
            rtd = False

        if values.get("embedding_charges", False) and rtd is False:
            raise ValueError("Cannot return interaction data when using embedding scheme.")

        return rtd

    def build_tasks(
        self,
        mb_computer: SubTaskComputers,
        mc_level_idx: int,
        **kwargs: Dict[str, Any],
    ) -> int:
        """Adds to the task_list as many new unique tasks as necessary to treat a single model chemistry level at one or several n-body levels.
        New tasks are of type *mb_computer* with model chemistry level specified in *kwargs* and n-body levels accessed through *mc_level_idx*.

        Parameters
        ----------
        mb_computer
            Class of task computers to instantiate and add to self.task_list. Usually :class:`~psi4.driver.AtomicComputer` but may be other when wrappers are layered.
        mc_level_idx
            Position in field self.nbodies_per_mc_level used to obtain ``nbodies``, the list of n-body
            levels (e.g., `[1]` or `[1, 2]` or `["supersystem"]`) to which the modelchem specified in **kwargs** applies.
            That is, `nbodies = self.nbodies_per_mc_level[mc_level_idx]`.
            Note the natural 1-indexing of ``nbodies`` _contents_, so `[1]` covers one-body contributions.
            The corresponding user label is the 1-indexed counterpart, `mc_level_lbl = mc_level_idx + 1`
            Formerly nlevel as in `nbody = self.nbody_list[nbody_level=nlevel]`.
        kwargs
            Other arguments for initializing **mb_computer**. In particular, specifies model chemistry.

        Returns
        -------
        count : int
            Number of new tasks planned by this call.
            Formerly, didn't include supersystem in count.

        """
        # Get the n-body orders for this level
        nbodies = self.nbodies_per_mc_level[mc_level_idx]

        info = "\n" + p4util.banner(f" ManyBody Setup: N-Body Levels {nbodies}", strNotOutfile=True) + "\n"
        core.print_out(info)
        logger.info(info)

        for kwg in ['dft_functional']:
            if kwg in kwargs:
                kwargs['keywords']['function_kwargs'][kwg] = kwargs.pop(kwg)

        count = 0
        template = copy.deepcopy(kwargs)

        # Get compute list
        if nbodies == ["supersystem"]:
            # Add supersystem computation if requested -- always nocp
            data = template
            data["molecule"] = self.molecule
            key = f"supersystem_{self.nfragments}"
            self.task_list[key] = mb_computer(**data)
            count += 1

            compute_dict = build_nbody_compute_list(
                ["nocp"], list(range(1, self.max_nbody + 1)), self.nfragments, self.return_total_data
            )
        else:
            compute_dict = build_nbody_compute_list(self.bsse_type, nbodies, self.nfragments, self.return_total_data)

        def labeler(item) -> str:
            mc_level_lbl = mc_level_idx + 1
            return str(mc_level_lbl) + "_" + str(item)

        # Add current compute list to the master task list
        # * `pair` looks like `((1,), (1, 3))` where first is real (not ghost) fragment indices
        #    and second is basis set fragment indices, all 1-indexed
        for nb in compute_dict["all"]:
            for pair in compute_dict["all"][nb]:
                lbl = labeler(pair)
                if lbl in self.task_list:
                    continue

                data = template
                ghost = list(set(pair[1]) - set(pair[0]))
                data["molecule"] = self.molecule.extract_subsets(list(pair[0]), ghost)
                if self.embedding_charges:
                    embedding_frags = list(set(range(1, self.nfragments + 1)) - set(pair[1]))
                    charges = []
                    for frag in embedding_frags:
                        positions = self.molecule.extract_subsets(frag).geometry().np.tolist()
                        charges.extend([[chg, i] for i, chg in zip(positions, self.embedding_charges[frag])])
                    data['keywords']['function_kwargs'].update({'external_potentials': charges})

                self.task_list[lbl] = mb_computer(**data)
                count += 1

        return count

    def plan(self):
        # uncalled function
        return [t.plan() for t in self.task_list.values()]

    def compute(self, client: Optional["qcportal.FractalClient"] = None):
        """Run quantum chemistry."""

        info = "\n" + p4util.banner(f" ManyBody Computations ", strNotOutfile=True) + "\n"
        #core.print_out(info)
        logger.info(info)

        with p4util.hold_options_state():
            for t in self.task_list.values():
                t.compute(client=client)

    def prepare_results(
        self,
        results: Optional[Dict[str, SubTaskComputers]] = None,
        client: Optional["qcportal.FractalClient"] = None,
    ) -> Dict[str, Any]:
        """Process the results from all n-body component molecular systems and model chemistry levels into final quantities.

        Parameters
        ----------
        results
            A set of tasks to process instead of self.task_list. Used in multilevel processing to pass a subset of
            self.task_list filtered to only one modelchem level.
        client
            QCFractal client if using QCArchive for distributed compute.

        Returns
        -------
        nbody_results
            When the ManyBodyComputer specifies a single model chemistry level (see self.nbodies_per_mc_level), the
            return is a dictionary, nbody_results, described in the table below. Many of the items are actually filled
            by successive calls to assemble_nbody_components(). When multiple model chemistry levels are specified, this
            function diverts its return to driver_nbody_multilevel.prepare_results() wherein each mc level calls this
            function again and collects separate nbody_results dictionaries and processes them into a final return that
            is a small subset of the table below.


                                       ptype_size = (1,)/(nat, 3)/(3 * nat, 3 * nat)
                                        e/g/h := energy or gradient or Hessian
                                        rtd := return_total_data

        .. |em| unicode:: U+02003 .. em space

        .. _`table:nbody_return`:

        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | item                                                          | size                 | present / zeroed                                                   | contents / interpretation                                                                                          |
        +===============================================================+======================+====================================================================+====================================================================================================================+
        | ret_ptype                                                     | ptype_size           | always                                                             | interaction data requested: IE or total (depending on return_total_data) e/g/h (depending on driver)               |
        |                                                               |                      |                                                                    |   with cp/nocp/vmfc treatment (depending on 1st of bsse_type)                                                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | ret_energy                                                    | 1                    | always                                                             | interaction energy: IE or total (depending on return_total_data) w/ cp/nocp/vmfc treat. (dep. on 1st of bsse_type) |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | ret_gradient                                                  | (nat, 3)             | when driver is g/h                                                 | interaction gradient: IE or total (depending on return_total_data) w/ cp/nocp/vmfc treat. (dep. on 1st of bsse_type|
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | ret_hessian                                                   | (nat * 3, nat * 3)   | when driver is h                                                   | interaction Hessian: IE or total (depending on return_total_data) w/ cp/nocp/vmfc treat. (dep. on 1st of bsse_type)|
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | nbody                                                         | >=1                  | always                                                             | energy n-body QCVariables to be set                                                                                |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY               |  |em| 1              | when cp in bsse_type                                               | MBE sum of subsystems of 1-body. summed are total energies with cp treatment                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY               |  |em| 1              | when cp in bsse_type & max_nbody>=2                                | MBE sum of subsystems of 2-body or fewer (cumulative); summed are total energies with cp treatment                 |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY            |  |em| 1              | when cp in bsse_type                                               | MBE sum of subsystems of {max_nbody}-body or fewer (cumulative); summed are total energies w/ cp treatment         |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED TOTAL ENERGY                              |  |em| 1              | when cp in bsse_type & rtd=T                                       | best available total energy with cp treatment: CP-CORRECTED TOTAL ENERGY THROUGH {max_nbody}-BODY                  |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY         |  |em| 1              | when cp in bsse_type & max_nbody>=2                                | 2-body total data less 1-body total data for cumulative IE; inputs are total energies with cp treatment            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED INTERACTION ENERGY THROUGH {nb}-BODY      |  |em| 1              | when cp in bsse_type                                               | {max_nbody}-body total data less 1-body total data for cumulative IE; inputs are total energies with cp treatment  |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED INTERACTION ENERGY                        |  |em| 1              | when cp in bsse_type                                               | best available interaction energy with cp treatment: CP-CORRECTED INTERACTION ENERGY THROUGH {max_nbody}-BODY      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY             |  |em| 1              | when cp in bsse_type & max_nbody>=2                                | 2-body total data less (2-1)-body total data for partial IE; inputs are total energies w/ cp treatment             |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY          |  |em| 1              | when cp in bsse_type                                               | {max_nbody}-body total data less ({max_nbody}-1)-body data for partial IE; inputs are total energies w/ cp treat.  |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY             |  |em| 1              | when nocp in bsse_type                                             | MBE sum of subsystems of 1-body. summed are total energies without cp treatment                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY             |  |em| 1              | when nocp in bsse_type & max_nbody>=2                              | MBE sum of subsystems of 2-body or fewer (cumulative); summed are total energies without cp treatment              |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY          |  |em| 1              | when nocp in bsse_type                                             | MBE sum of subsystems of {max_nbody}-body or fewer (cumulative); summed are total energies w/o cp treatment        |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED TOTAL ENERGY                            |  |em| 1              | when nocp in bsse_type                                             | best available total energy without cp treatment: NOCP-CORRECTED TOTAL ENERGY THROUGH {max_nbody}-BODY             |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY       |  |em| 1              | when nocp in bsse_type & max_nbody>=2                              | 2-body total data less 1-body total data for cumulative IE; inputs are total energies w/o cp treatment             |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED INTERACTION ENERGY THROUGH {nb}-BODY    |  |em| 1              | when nocp in bsse_type                                             | {max_nbody}-body total data less 1-body total data for cumulative IE; inputs are total energies w/o cp treatment   |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED INTERACTION ENERGY                      |  |em| 1              | when nocp in bsse_type                                             | best available interaction energy without cp treatment: NOCP-CORRECTED INTERACTION ENERGY THROUGH {max_nbody}-BODY |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY           |  |em| 1              | when nocp in bsse_type & max_nbody>=2                              | 2-body total data less (2-1)-body total data for partial IE; inputs are total energies w/o cp treatment            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY        |  |em| 1              | when nocp in bsse_type                                             | {max_nbody}-body total data less ({max_nbody}-1)-body data for partial IE; inputs are total energies w/o cp treat. |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY             |  |em| 1              | when vmfc in bsse_type                                             | MBE sum of subsystems of 1-body. summed are total energies with vmfc treatment                                     |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY             |  |em| 1              | when vmfc in bsse_type & max_nbody>=2                              | MBE sum of subsystems of 2-body or fewer (cumulative); summed are total energies with vmfc treatment               |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY          |  |em| 1              | when vmfc in bsse_type                                             | MBE sum of subsystems of {max_nbody}-body or fewer (cumulative); summed are total energies w/ vmfc treatment       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED TOTAL ENERGY                            |  |em| 1              | when vmfc in bsse_type                                             | best available total energy with vmfc treatment: VMFC-CORRECTED TOTAL ENERGY THROUGH {max_nbody}-BODY              |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY       |  |em| 1              | when vmfc in bsse_type & max_nbody>=2                              | 2-body total data less 1-body total data for cumulative IE; inputs are total energies w/ vmfc treatment            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED INTERACTION ENERGY THROUGH {nb}-BODY    |  |em| 1              | when vmfc in bsse_type                                             | {max_nbody}-body total data less 1-body total data for cumulative IE; inputs are total energies w/ vmfc treatment  |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED INTERACTION ENERGY                      |  |em| 1              | when vmfc in bsse_type                                             | best available interaction energy with vmfc treatment: VMFC-CORRECTED INTERACTION ENERGY THROUGH {max_nbody}-BODY  |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED 2-BODY CONTRIBUTION TO ENERGY           |  |em| 1              | when vmfc in bsse_type & max_nbody>=2                              | 2-body total data less (2-1)-body total data for partial IE; inputs are total energies w/ vmfc treatment           |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY        |  |em| 1              | when vmfc in bsse_type                                             | {max_nbody}-body total data less ({max_nbody}-1)-body data for partial IE; inputs are total energies w/ vmfc treat.|
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | energy_body_dict                                              | max_nbody            | always                                                             | total energies at each n-body level                                                                                |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| 1                                                       |  |em| 1              | always; zeroed if cp & rtd=F                                       | cumulative through 1-body total energies w/ cp/nocp/vmfc treatment (dep. on 1st of bsse_type)                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| 2                                                       |  |em| 1              | max_nbody>=2                                                       | cumulative through 2-body total energies w/ cp/nocp/vmfc treatment (dep. on 1st of bsse_type)                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| {max_nbody}                                             |  |em| 1              | always                                                             | cumulative through {max_nbody}-body total energies w/ cp/nocp/vmfc treatment (dep. on 1st of bsse_type)            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | gradient_body_dict                                            | max_nbody            | when driver is g/h                                                 |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| 1                                                       |  |em| (nat, 3)       | when driver is g/h                                                 | cumulative through 1-body total gradients with cp/nocp/vmfc treatment (dep. on 1st of bsse_type)                   |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| 2                                                       |  |em| (nat, 3)       | when driver is g/h & max_nbody>=2                                  | cumulative through 2-body total gradients with cp/nocp/vmfc treatment (dep. on 1st of bsse_type)                   |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| {max_nbody}                                             |  |em| (nat, 3)       | when driver is g/h                                                 | cumulative through {max_nbody}-body total gradients w/ cp/nocp/vmfc treatment (dep. on 1st of bsse_type)           |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | hessian_body_dict                                             | max_nbody            | when driver is h                                                   |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| 1                                                       |  |em| (nat*3, nat*3) | when driver is h                                                   | cumulative through 1-body total Hessians w/ cp/nocp/vmfc treatment (dep. on 1st of bsse_type)                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| 2                                                       |  |em| (nat*3, nat*3) | when driver is h & max_nbody>=2                                    | cumulative through 2-body total Hessians w/ cp/nocp/vmfc treatment (dep. on 1st of bsse_type)                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| {max_nbody}                                             |  |em| (nat*3, nat*3) | when driver is h                                                   | cumulative through {max_nbody}-body total Hessians w/ cp/nocp/vmfc treatment (dep. on 1st of bsse_type)            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | cp_energy_body_dict                                           | max_nbody            | always; zeroed if cp not in bsse_type                              | total energies at each n-body level with cp treatment                                                              |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1cp                                                    |  |em| 1              | always; zeroed if cp not in bsse_type or rtd=F                     | cumulative through 1-body total energies with cp treatment                                                         |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2cp                                                    |  |em| 1              | when max_nbody>=2; zeroed if cp not in bsse_type                   | cumulative through 2-body total energies with cp treatment                                                         |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}cp                                          |  |em| 1              | always; zeroed if cp not in bsse_type                              | cumulative through {max_nbody}-body total energies with cp treatment                                               |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | cp_gradient_body_dict                                         | max_nbody            | when driver is g/h; zeroed if cp not in bsse_type                  | total gradients at each n-body level with cp treatment                                                             |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1cp                                                    |  |em| (nat, 3)       | when driver is g/h; zeroed if cp not in bsse_type or rtd=F         | cumulative through 1-body total gradients with cp treatment                                                        |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2cp                                                    |  |em| (nat, 3)       | when driver is g/h & max_nbody>=2; zeroed if cp not in bsse_type   | cumulative through 2-body total gradients with cp treatment                                                        |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}cp                                          |  |em| (nat, 3)       | when driver is g/h; zeroed if cp not in bsse_type                  | cumulative through {max_nbody}-body total gradients with cp treatment                                              |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | cp_hessian_body_dict                                          | max_nbody            | when driver is h; zeroed if cp not in bsse_type                    | total Hessians at each n-body level with cp treatment                                                              |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1cp                                                    |  |em| (nat*3, nat*3) | when driver is h; zeroed if cp not in bsse_type or rtd=F           | cumulative through 1-body total Hessians with cp treatment                                                         |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2cp                                                    |  |em| (nat*3, nat*3) | when driver is h & max_nbody>=2; zeroed if cp not in bsse_type     | cumulative through 2-body total Hessians with cp treatment                                                         |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}cp                                          |  |em| (nat*3, nat*3) | when driver is h; zeroed if cp not in bsse_type                    | cumulative through {max_nbody}-body total Hessians with cp treatment                                               |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | nocp_energy_body_dict                                         | max_nbody            | always; zeroed if nocp not in bsse_type                            | total energies at each n-body level with nocp treatment                                                            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1nocp                                                  |  |em| 1              | always; zeroed if nocp not in bsse_type                            | cumulative through 1-body total energies with nocp treatment                                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2nocp                                                  |  |em| 1              | when max_nbody>=2; zeroed if nocp not in bsse_type                 | cumulative through 2-body total energies with nocp treatment                                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}nocp                                        |  |em| 1              | always; zeroed if nocp not in bsse_type                            | cumulative through {max_nbody}-body total energies with nocp treatment                                             |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | nocp_gradient_body_dict                                       | max_nbody            | when driver is g/h; zeroed if nocp not in bsse_type                | total gradients at each n-body level with nocp treatment                                                           |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1nocp                                                  |  |em| (nat, 3)       | when driver is g/h; zeroed if nocp not in bsse_type                | cumulative through 1-body total gradients with nocp treatment                                                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2nocp                                                  |  |em| (nat, 3)       | when driver is g/h & max_nbody>=2; zeroed if nocp not in bsse_type | cumulative through 2-body total gradients with nocp treatment                                                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}nocp                                        |  |em| (nat, 3)       | when driver is g/h; zeroed if nocp not in bsse_type                | cumulative through {max_nbody}-body total gradients with nocp treatment                                            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | nocp_hessian_body_dict                                        | max_nbody            | when driver is h; zeroed if nocp not in bsse_type                  | total Hessians at each n-body level with nocp treatment                                                            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1nocp                                                  |  |em| (nat*3, nat*3) | when driver is h; zeroed if nocp not in bsse_type                  | cumulative through 1-body total Hessians with nocp treatment                                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2nocp                                                  |  |em| (nat*3, nat*3) | when driver is h & max_nbody>=2; zeroed if nocp not in bsse_type   | cumulative through 2-body total Hessians with nocp treatment                                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}nocp                                        |  |em| (nat*3, nat*3) | when driver is h; zeroed if nocp not in bsse_type                  | cumulative through {max_nbody}-body total Hessians with nocp treatment                                             |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | vmfc_energy_body_dict                                         | max_nbody            | always; zeroed if vmfc not in bsse_type                            | total energies at each n-body level with vmfc treatment                                                            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1vmfc                                                  |  |em| 1              | always; zeroed if vmfc not in bsse_type                            | cumulative through 1-body total energies with vmfc treatment                                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2vmfc                                                  |  |em| 1              | when max_nbody>=2; zeroed if vmfc not in bsse_type                 | cumulative through 2-body total energies with vmfc treatment                                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}vmfc                                        |  |em| 1              | always; zeroed if vmfc not in bsse_type                            | cumulative through {max_nbody}-body total energies with vmfc treatment                                             |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | vmfc_gradient_body_dict                                       | max_nbody            | when driver is g/h; zeroed if vmfc not in bsse_type                | total gradients at each n-body level with vmfc treatment                                                           |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1vmfc                                                  |  |em| (nat, 3)       | when driver is g/h; zeroed if vmfc not in bsse_type                | cumulative through 1-body total gradients with vmfc treatment                                                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2vmfc                                                  |  |em| (nat, 3)       | when driver is g/h & max_nbody>=2; zeroed if vmfc not in bsse_type | cumulative through 2-body total gradients with vmfc treatment                                                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}vmfc                                        |  |em| (nat, 3)       | when driver is g/h; zeroed if vmfc not in bsse_type                | cumulative through {max_nbody}-body total gradients with vmfc treatment                                            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | vmfc_hessian_body_dict                                        | max_nbody            | when driver is h; zeroed if vmfc not in bsse_type                  | total Hessians at each n-body level with vmfc treatment                                                            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1vmfc                                                  |  |em| (nat*3, nat*3) | when driver is h; zeroed if vmfc not in bsse_type                  | cumulative through 1-body total Hessians with vmfc treatment                                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2vmfc                                                  |  |em| (nat*3, nat*3) | when driver is h & max_nbody>=2; zeroed if vmfc not in bsse_type   | cumulative through 2-body total Hessians with vmfc treatment                                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}vmfc                                        |  |em| (nat*3, nat*3) | when driver is h; zeroed if vmfc not in bsse_type                  | cumulative through {max_nbody}-body total Hessians with vmfc treatment                                             |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | intermediates                                                 | ntasks               | always                                                             | all individual energies with nice labels                                                                           |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| N-BODY (1, 2)@(1, 2) TOTAL ENERGY                      |  |em| 1              | always                                                             | total energy for 1st modelchem, 1st & 2nd fragments in basis of 1st & 2nd fragments                                |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| N-BODY (3)@(2, 3) TOTAL ENERGY                         |  |em| 1              | always                                                             | total energy for 2nd modelchem, 3rd fragment in basis of 2nd and 3rd fragments                                     |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| ...                                                    |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | intermediates_energy                                          | ntasks               | always                                                             | all individual energies                                                                                            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1_((1, 2), (1, 2))                                     |  |em| 1              | always                                                             | total energy for 1st modelchem, 1st & 2nd fragments in basis of 1st & 2nd fragments                                |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2_((3,), (2, 3))                                       |  |em| 1              | always                                                             | total energy for 2nd modelchem, 3rd fragment in basis of 2nd and 3rd fragments                                     |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| ...                                                    |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | intermediates_gradient                                        | ntasks               | when driver is g/h                                                 | all individual gradients                                                                                           |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1_((1, 2), (1, 2))                                     |  |em| (nat, 3)       | when driver is g/h                                                 | total gradient for 1st modelchem, 1st & 2nd fragments in basis of 1st & 2nd fragments                              |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2_((3,), (2, 3))                                       |  |em| (nat, 3)       | when driver is g/h                                                 | total gradient for 2nd modelchem, 3rd fragment in basis of 2nd and 3rd fragments                                   |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| ...                                                    |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | intermediates_hessian                                         | ntasks               | when driver is h                                                   | all individual Hessians                                                                                            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1_((1, 2), (1, 2))                                     |  |em| (nat*3, nat*3) | when driver is h                                                   | total Hessian for 1st modelchem, 1st & 2nd fragments in basis of 1st & 2nd fragments                               |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2_((3,), (2, 3))                                       |  |em| (nat*3, nat*3) | when driver is h                                                   | total Hessian for 2nd modelchem, 3rd fragment in basis of 2nd and 3rd fragments                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| ...                                                    |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+

        """
        if results is None:
            results = {}

        # formerly nlevels
        mc_level_labels = {i.split("_")[0] for i in self.task_list}
        if len(mc_level_labels) > 1 and not results:
            return driver_nbody_multilevel.prepare_results(self, client)

        results_list = {k: v.get_results(client=client) for k, v in (results.items() or self.task_list.items())}
        trove = {  # AtomicResult.properties return None if missing
            "energy": {k: v.properties.return_energy for k, v in results_list.items()},
            "gradient": {k: v.properties.return_gradient for k, v in results_list.items()},
            "hessian": {k: v.properties.return_hessian for k, v in results_list.items()},
        }

        # TODO: make assemble_nbody_components and driver_nbody_multilevel.prepare_results into class functions.
        #   note that the former uses metadata as read-only (except for one solveable case) while the latter overwrites self (!).
        metadata = {
            "quiet": self.quiet,
            "nbodies_per_mc_level": self.nbodies_per_mc_level,
            "bsse_type": self.bsse_type,
            "nfragments": self.nfragments,
            "return_total_data": self.return_total_data,
            "molecule": self.molecule,
            "embedding_charges": bool(self.embedding_charges),
            "max_nbody": self.max_nbody,
        }
        if self.driver.name == "energy":
            nbody_results = assemble_nbody_components("energy", trove["energy"], metadata.copy())

        elif self.driver.name == "gradient":
            nbody_results = assemble_nbody_components("energy", trove["energy"], metadata.copy())
            nbody_results.update(assemble_nbody_components("gradient", trove["gradient"], metadata.copy()))

        elif self.driver.name == "hessian":
            nbody_results = assemble_nbody_components("energy", trove["energy"], metadata.copy())
            nbody_results.update(assemble_nbody_components("gradient", trove["gradient"], metadata.copy()))
            nbody_results.update(assemble_nbody_components("hessian", trove["hessian"], metadata.copy()))

        def delabeler(item: str, return_obj: bool = False) -> Union[Tuple[str, str, str], Tuple[int, Tuple[int], Tuple[int]]]:
            """Transform labels like string "1_((2,), (1, 2))" into string tuple ("1", "2", "1, 2") or object tuple (1, (2,), (1, 2))."""

            mc, _, fragbas = item.partition("_")
            frag, bas = literal_eval(fragbas)

            if return_obj:
                return int(mc), frag, bas
            else:
                return mc, ", ".join(map(str, frag)), ", ".join(map(str, bas))

        # save some mc_(frag, bas) component results
        # * formerly, intermediates_energy was intermediates2
        # * formerly, intermediates_gradient was intermediates_ptype
        # * formerly, intermediates_hessian was intermediates_ptype

        nbody_results["intermediates"] = {}
        for idx, task in results_list.items():
            mc, frag, bas = delabeler(idx)
            nbody_results["intermediates"][f"N-BODY ({frag})@({bas}) TOTAL ENERGY"] = task.properties.return_energy

        nbody_results["intermediates_energy"] = trove["energy"]

        if not all(x is None for x in trove["gradient"].values()):
            nbody_results["intermediates_gradient"] = trove["gradient"]

        if not all(x is None for x in trove["hessian"].values()):
            nbody_results["intermediates_hessian"] = trove["hessian"]

        debug = False
        if debug:
            for k, v in nbody_results.items():
                if isinstance(v, np.ndarray):
                    print(f"CLS-prepared results >>> {k} {v.size}")
                elif isinstance(v, dict):
                    print(f"CLS-prepared results >>> {k} {len(v)}")
                    for k2, v2 in v.items():
                        if isinstance(v2, np.ndarray):
                            print(f"CLS-prepared results      >>> {k2} {v2.size}")
                        else:
                            print(f"CLS-prepared results      >>> {k2} {v2}")
                else:
                    print(f"CLS-prepared results >>> {k} {v}")

        return nbody_results

    def get_results(self, client: Optional["qcportal.FractalClient"] = None) -> AtomicResult:
        """Return results as ManyBody-flavored QCSchema."""

        info = "\n" + p4util.banner(f" ManyBody Results ", strNotOutfile=True) + "\n"
        core.print_out(info)
        logger.info(info)

        results = self.prepare_results(client=client)
        ret_energy = results.pop("ret_energy")
        ret_ptype = results.pop("ret_ptype")
        ret_gradient = results.pop("ret_gradient", None)

        # load QCVariables
        qcvars = {
            'NUCLEAR REPULSION ENERGY': self.molecule.nuclear_repulsion_energy(),
            'NBODY NUMBER': len(self.task_list),
        }

        properties = {
            "calcinfo_natom": self.molecule.natom(),
            "nuclear_repulsion_energy": self.molecule.nuclear_repulsion_energy(),
            "return_energy": ret_energy,
        }

        for k, val in results.items():
            qcvars[k] = val

        qcvars['CURRENT ENERGY'] = ret_energy
        if self.driver == 'gradient':
            qcvars['CURRENT GRADIENT'] = ret_ptype
            properties["return_gradient"] = ret_ptype
        elif self.driver == 'hessian':
            qcvars['CURRENT GRADIENT'] = ret_gradient
            qcvars['CURRENT HESSIAN'] = ret_ptype
            properties["return_gradient"] = ret_gradient
            properties["return_hessian"] = ret_ptype

        component_results = self.dict()['task_list']
        for k, val in component_results.items():
            val['molecule'] = val['molecule'].to_schema(dtype=2)

        nbody_model = AtomicResult(
            **{
                'driver': self.driver,
                'model': {
                    'method': self.method,
                    'basis': self.basis,
                },
                'molecule': self.molecule.to_schema(dtype=2),
                'properties': properties,
                'provenance': p4util.provenance_stamp(__name__),
                'extras': {
                    'qcvars': qcvars,
                    'component_results': component_results,
                },
                'return_result': ret_ptype,
                'success': True,
            })

        logger.debug('\nNBODY QCSchema:\n' + pp.pformat(nbody_model.dict()))

        return nbody_model

    def get_psi_results(
        self,
        client: Optional["qcportal.FractalClient"] = None,
        *,
        return_wfn: bool = False) -> EnergyGradientHessianWfnReturn:
        """Called by driver to assemble results into ManyBody-flavored QCSchema,
        then reshape and return them in the customary Psi4 driver interface: ``(e/g/h, wfn)``.

        Parameters
        ----------
        return_wfn
            Whether to additionally return the dummy :py:class:`~psi4.core.Wavefunction`
            calculation result as the second element of a tuple. Contents are:

            - supersystem molecule
            - dummy basis, def2-svp
            - e/g/h member data
            - QCVariables

        Returns
        -------
        ret
            Energy, gradient, or Hessian according to self.driver.
        wfn
            Wavefunction described above when *return_wfn* specified.

        """
        nbody_model = self.get_results(client=client)
        ret = nbody_model.return_result

        wfn = core.Wavefunction.build(self.molecule, "def2-svp", quiet=True)

        # TODO all besides nbody may be better candidates for extras than qcvars. energy/gradient/hessian_body_dict in particular are too simple for qcvars (e.g., "2")
        dicts = [
            #"energies",  # retired
            #"ptype",     # retired
            "intermediates",
            "intermediates_energy",  #"intermediates2",
            "intermediates_gradient",  #"intermediates_ptype",
            "intermediates_hessian",  #"intermediates_ptype",
            "energy_body_dict",
            "gradient_body_dict",  # ptype_body_dict
            "hessian_body_dict",  # ptype_body_dict
            "nbody",
            "cp_energy_body_dict",
            "nocp_energy_body_dict",
            "vmfc_energy_body_dict",
            "cp_gradient_body_dict",
            "nocp_gradient_body_dict",
            "vmfc_gradient_body_dict",
            "cp_hessian_body_dict",
            "nocp_hessian_body_dict",
            "vmfc_hessian_body_dict",
        ]

        for qcv, val in nbody_model.extras['qcvars'].items():
            if isinstance(val, dict):
                if qcv in dicts:
                    for qcv2, val2 in val.items():
                        for obj in [core, wfn]:
                            try:
                                obj.set_variable(str(qcv2), val2)
                            except ValidationError:
                                obj.set_variable(f"{self.driver.name} {qcv2}", val2)
            else:
                for obj in [core, wfn]:
                    obj.set_variable(qcv, val)

        if self.driver == 'gradient':
            ret = core.Matrix.from_array(ret)
            wfn.set_gradient(ret)
        elif self.driver == 'hessian':
            ret = core.Matrix.from_array(ret)
            grad = core.Matrix.from_array(nbody_model.properties.return_gradient)
            wfn.set_hessian(ret)
            wfn.set_gradient(grad)

        if return_wfn:
            return (ret, wfn)
        else:
            return ret


# TODO questions to check:
# * can work with supersystem and embedding_charges?
# * can levels work with same method, different basis?
