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

__all__ = [
    "expand_cbs_methods",
    "task_planner",
    "TaskComputers",
]

import copy
import logging
import os
from typing import Dict, Tuple, Union

from qcelemental.models import DriverEnum

from psi4 import core

from . import p4util
from .constants import pp
from .driver_cbs import CompositeComputer, cbs_text_parser, composite_procedures
from .driver_findif import FiniteDifferenceComputer
from .driver_nbody import ManyBodyComputer
from .driver_util import negotiate_convergence_criterion, negotiate_derivative_type
from .task_base import AtomicComputer

logger = logging.getLogger(__name__)

TaskComputers = Union[AtomicComputer, CompositeComputer, FiniteDifferenceComputer, ManyBodyComputer]


def expand_cbs_methods(method: str, basis: str, driver: DriverEnum, **kwargs) -> Tuple[str, str, Dict]:
    """Sort out the user input method string into recognized fields.
    Handles cases like:

    (i) ``"mp2"`` -- passes through;
    (ii) ``"mp2/cc-pvdz"`` -- broken into method and basis fields;
    (iii) ``"mp2/cc-pv[d,t]z"`` -- processed into method="cbs" & CBSMetadata spec;
    (iv) ``method="cbs", cbsmeta=CBSMetadata`` -- passes through.

    Parameters
    ----------
    method
        User first argument to driver function. A string hint of the method --
        see cases above.
    basis
        User basis hint.
    driver
        The calling driver function. Note for finite difference that this is
        the target driver, not the means driver.

    """
    if method == 'cbs' and kwargs.get('cbsmeta', False):
        return method, basis, kwargs['cbsmeta']

    # Expand CBS methods
    if "/" in method:
        kwargs["ptype"] = driver
        cbsmeta = cbs_text_parser(method, **kwargs)

        # Single call detected
        if "cbs_metadata" not in cbsmeta:
            method = cbsmeta["method"]
            basis = cbsmeta["basis"]
        else:
            method = "cbs"
    else:
        cbsmeta = {}

    return method, basis, cbsmeta


def task_planner(driver: DriverEnum, method: str, molecule: core.Molecule, **kwargs) -> TaskComputers:
    """Plans a task graph of a complex computation.

    Canonical Task layering:
     - ManyBody - BSSE treatment, many-body expansion
     - FiniteDifference - derivatives through stencils
     - Composite - basis set extrapolation, focal-point methods
     - Atomic - analytic single-points

    Parameters
    ----------
    driver
        The resulting type of computation: e/g/h. Note for finite difference
        that this should be the target driver, not the means driver.
    method
        A string representation of the method such as "HF" or "B3LYP". Special
        cases are: "cbs".
    molecule
        A Psi4 base molecule to use.
    kwargs
        User keyword arguments, often used to configure task computers.

    Returns
    -------
    Union[AtomicComputer, CompositeComputer, FiniteDifferenceComputer, ManyBodyComputer]
        A simple (:class:`~psi4.driver.AtomicComputer`) or layered (:class:`~psi4.driver.driver_cbs.CompositeComputer`, :class:`~psi4.driver.driver_findif.FiniteDifferenceComputer`, :class:`~psi4.driver.driver_nbody.ManyBodyComputer`) task object. Layered objects contain many and multiple types of computers in a graph.

    """

    # Only pull the changed options
    keywords = p4util.prepare_options_for_set_options()

    keywords["function_kwargs"] = {}
    if "external_potentials" in kwargs:
        keywords["function_kwargs"].update({"external_potentials": kwargs.pop("external_potentials")})

    # Need to add full path to pcm file
    if "PCM__PCMSOLVER_PARSED_FNAME" in keywords.keys():
        fname = keywords["PCM__PCMSOLVER_PARSED_FNAME"]
        keywords["PCM__PCMSOLVER_PARSED_FNAME"] = os.path.join(os.getcwd(), fname)

    # Pull basis out of kwargs, override globals if user specified
    basis = kwargs.pop("basis", keywords.pop("BASIS", "(auto)"))
    method = method.lower()

    # Expand CBS methods
    method, basis, cbsmeta = expand_cbs_methods(method, basis, driver, **kwargs)
    if method in composite_procedures:
        kwargs.update({'cbs_metadata': composite_procedures[method](**kwargs)})
        method = 'cbs'

    pertinent_findif_kwargs = ['findif_irrep', 'findif_stencil_size', 'findif_step_size', 'findif_verbose']
    current_findif_kwargs = {kw: kwargs.pop(kw) for kw in pertinent_findif_kwargs if kw in kwargs}
    # explicit: 'findif_mode'

    # Build a packet
    packet = {"molecule": molecule, "driver": driver, "method": method, "basis": basis, "keywords": keywords}

    # First check for BSSE type
    if kwargs.get("bsse_type", None) is not None:
        levels = kwargs.pop('levels', None)

        plan = ManyBodyComputer(**packet, **kwargs)
        original_molecule = packet.pop("molecule")

        # Add tasks for every nbody level requested
        if levels is None:
            levels = {plan.max_nbody: method}
        else:
            # rearrange bodies in order with supersystem last lest body count fail in organization loop below
            levels = dict(sorted(levels.items(), key=lambda item: 1000 if item[0] == "supersystem" else item[0]))

            # We define cp as being a correction to only interaction energies
            # If only doing cp, we need to ignore any user-specified 1st (monomer) level
            if 'cp' in kwargs.get("bsse_type", None) and 'nocp' not in kwargs.get("bsse_type", None): 
                if 1 in levels.keys():
                    removed_level = levels.pop(1)
                    logger.info("NOTE: User specified exclusively 'cp' correction, but provided level 1 details") 
                    logger.info(f"NOTE: Removing level {removed_level}")
                    logger.info("NOTE: For total energies, add 'nocp' to bsse_list")


        # Organize nbody calculations into modelchem levels
        # * expand keys of `levels` into full lists of nbodies covered. save to plan, resetting max_nbody accordingly
        # * below, process values of `levels`, which are modelchem strings, into kwargs specs
        nbodies_per_mc_level = []
        prev_body = 0
        for nb in levels:
            nbodies = []
            if nb == "supersystem":
                nbodies.append(nb)
            elif nb != (prev_body + 1):
                for m in range(prev_body + 1, nb + 1):
                    nbodies.append(m)
            else:
                nbodies.append(nb)
            nbodies_per_mc_level.append(nbodies)
            prev_body += 1

        plan.max_nbody = max(nb for nb in levels if nb != "supersystem")
        plan.nbodies_per_mc_level = nbodies_per_mc_level

        for mc_level_idx, mtd in enumerate(levels.values()):
            method, basis, cbsmeta = expand_cbs_methods(mtd, basis, driver, cbsmeta=cbsmeta, **kwargs)
            packet.update({'method': method, 'basis': basis})

            # Tell the task builder which level to add a task list for
            # * see https://github.com/psi4/psi4/pull/1351#issuecomment-549948276 for discussion of where build_tasks logic should live
            if method == "cbs":
                # This CompositeComputer is discarded after being used for dermode.
                simplekwargs = copy.deepcopy(kwargs)
                simplekwargs.pop('dertype', None)
                simplecbsmeta = copy.deepcopy(cbsmeta)
                simplecbsmeta['verbose'] = 0
                dummyplan = CompositeComputer(**packet, **simplecbsmeta, molecule=original_molecule, **simplekwargs)

                methods = [sr.method for sr in dummyplan.task_list]
                # TODO: pass more info, so fn can use for managed_methods -- ref, qc_module, fc/ae, conv/df
                dermode = negotiate_derivative_type(driver, methods, kwargs.pop('dertype', None), verbose=1)

                if dermode[0] == dermode[1]:  # analytic
                    logger.info("PLANNING MB(CBS):  {mc_level_idx=} {packet=} {cbsmeta=} kw={kwargs}")
                    plan.build_tasks(CompositeComputer, **packet, mc_level_idx=mc_level_idx, **cbsmeta, **kwargs)

                else:
                    logger.info(
                        f"PLANNING MB(FD(CBS):  {mc_level_idx=} {packet=} {cbsmeta=} findif_kw={current_findif_kwargs} kw={kwargs}"
                    )
                    plan.build_tasks(FiniteDifferenceComputer,
                                     **packet,
                                     mc_level_idx=mc_level_idx,
                                     findif_mode=dermode,
                                     computer=CompositeComputer,
                                     **cbsmeta,
                                     **current_findif_kwargs,
                                     **kwargs)

            else:
                dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
                if dermode[0] == dermode[1]:  # analytic
                    logger.info(f"PLANNING MB:  {mc_level_idx=} {packet=}")
                    plan.build_tasks(AtomicComputer, **packet, mc_level_idx=mc_level_idx, **kwargs)
                else:
                    logger.info(
                        f"PLANNING MB(FD):  {mc_level_idx=} {packet=} findif_kw={current_findif_kwargs} kw={kwargs}"
                    )
                    plan.build_tasks(FiniteDifferenceComputer,
                                     **packet,
                                     mc_level_idx=mc_level_idx,
                                     findif_mode=dermode,
                                     **current_findif_kwargs,
                                     **kwargs)

        return plan

    # Check for CBS
    elif method == "cbs":
        kwargs.update(cbsmeta)
        # This CompositeComputer is discarded after being used for dermode. Could have used directly for analytic except for excess printing with FD
        simplekwargs = copy.deepcopy(kwargs)
        simplekwargs['verbose'] = 0
        dummyplan = CompositeComputer(**packet, **simplekwargs)

        methods = [sr.method for sr in dummyplan.task_list]
        # TODO: pass more info, so fn can use for managed_methods -- ref, qc_module, fc/ae, conv/df
        dermode = negotiate_derivative_type(driver, methods, kwargs.pop('dertype', None), verbose=1)

        if dermode[0] == dermode[1]:  # analytic
            logger.info('PLANNING CBS:  packet={packet} kw={kwargs}')
            plan = CompositeComputer(**packet, **kwargs)
            return plan
        else:
            # For FD(CBS(Atomic)), the CompositeComputer above is discarded after being used for dermode.
            logger.info(
                f'PLANNING FD(CBS):  dermode={dermode} packet={packet} findif_kw={current_findif_kwargs} kw={kwargs}')
            plan = FiniteDifferenceComputer(**packet,
                                            findif_mode=dermode,
                                            computer=CompositeComputer,
                                            **current_findif_kwargs,
                                            **kwargs)
            return plan

    # Done with Wrappers -- know we want E, G, or H -- but may still be FD or AtomicComputer
    else:
        dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
        convcrit = negotiate_convergence_criterion(dermode, method, return_optstash=False)

        if dermode[0] == dermode[1]:  # analytic
            logger.info(f'PLANNING Atomic:  keywords={keywords}')
            return AtomicComputer(**packet, **kwargs)
        else:
            keywords.update(convcrit)
            logger.info(
                f'PLANNING FD:  dermode={dermode} keywords={keywords} findif_kw={current_findif_kwargs} kw={kwargs}')
            return FiniteDifferenceComputer(**packet,
                                            findif_mode=dermode,
                                            **current_findif_kwargs,
                                            **kwargs)
