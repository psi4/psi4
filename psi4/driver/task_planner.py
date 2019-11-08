#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
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

import abc
import copy
import math
import itertools
import pydantic
import pprint
pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)

from typing import Dict, List, Any, Union

import numpy as np

from psi4.driver import p4util
from psi4.driver.p4util.exceptions import UpgradeHelper
from psi4.driver.task_base import SingleComputer
from psi4.driver.driver_findif import FinDifComputer
from psi4.driver.driver_nbody import NBodyComputer
from psi4.driver.driver_cbs import CBSComputer, _cbs_text_parser
from psi4.driver.driver_util import negotiate_derivative_type, negotiate_convergence_criterion

__all__ = ["task_planner"]


def _expand_cbs_methods(method, basis, driver, **kwargs):
    if method == 'cbs' and kwargs.get('cbsmeta', None):
        return method, basis, kwargs['cbsmeta']

    # Expand CBS methods
    if "/" in method:
        kwargs["ptype"] = driver
        cbsmeta = _cbs_text_parser(method, **kwargs)

        # Single call detected
        if "cbs_metadata" not in cbsmeta:
            method = cbsmeta["method"]
            basis = cbsmeta["basis"]
        else:
            method = "cbs"
    else:
        cbsmeta = {}

    return method, basis, cbsmeta


def task_planner(driver, method, molecule, **kwargs):
    """Plans a task graph of a complex computations.

    Canonical Task layering:
     - NBody
     - FinDif
     - CBS
     - Single

    Parameters
    ----------
    driver : {"energy", "gradient", "hessian"}
        The resulting type of computation. That is, target, not means.
    method : str
        A string representation of the method such as "HF" or "B3LYP". Special cases are:
        "cbs"
    molecule : psi4.core.Molecule
        A Psi4 base molecule to use
    **kwargs
        Description

    Returns
    -------
    Task
        A task object

    """

    # Only pull the changed options
    keywords = p4util.prepare_options_for_set_options()
    keywords["function_kwargs"] = {}

    # Pull basis out of kwargs, override globals if user specified
    basis = keywords.pop("BASIS", "(auto)") #None)
    basis = kwargs.pop("basis", basis)
    method = method.lower()

    # Expand CBS methods
    method, basis, cbsmeta = _expand_cbs_methods(method, basis, driver, **kwargs)

    pertinent_findif_kwargs = ['findif_irrep', 'findif_stencil_size', 'findif_step_size', 'findif_verbose']
    current_findif_kwargs = {kw: kwargs.pop(kw) for kw in pertinent_findif_kwargs if kw in kwargs}
    # explicit: 'findif_mode'

    # Build a packet
    packet = {"molecule": molecule, "driver": driver, "method": method, "basis": basis, "keywords": keywords}

    # First check for BSSE type
    if 'bsse_type' in kwargs:
        levels = kwargs.pop('levels', None)

        plan = NBodyComputer(**packet, **kwargs)
        original_molecule = packet.pop("molecule")

        # Add tasks for every nbody level requested
        if levels is None:
            levels = {plan.max_nbody: method}
        
        # Organize nbody calculations into levels
        nbody_list = []
        prev_body = 0
        for n in levels:
            level = []
            if n == 'supersystem':
                level.append(n)
            elif n != (prev_body + 1):
                for m in range(prev_body + 1, n + 1):
                    level.append(m)
            else:
                level.append(n)
            nbody_list.append(level)
            prev_body += 1

        plan.max_nbody = max([i for i in levels if isinstance(i, int)])
        plan.nbody_list = nbody_list

        for nlevel, [n, level] in enumerate(levels.items()):
            method, basis, cbsmeta = _expand_cbs_methods(level, basis, driver, cbsmeta=cbsmeta, **kwargs)
            packet.update({'method': method, 'basis': basis})

      #      if n == 'supersytem':
      #          nlevel = plan.max_nbody

            # Tell the task bulider which level to add a task list for
            if method == "cbs":
                # This CBSComputer is discarded after being used for dermode.
                dummyplan = CBSComputer(**packet, **cbsmeta, molecule=original_molecule)

                methods = [sr.method for sr in dummyplan.task_list]
                # TODO: pass more info, so fn can use for managed_methods -- ref, qc_module, fc/ae, conv/df
                dermode = negotiate_derivative_type(driver, methods, kwargs.pop('dertype', None), verbose=1)

                if dermode[0] == dermode[1]:  # analytic
                    print('PLANNING NBody(CBS)', 'n=', n, nlevel, 'packet=', packet, 'cbsmeta=', cbsmeta, 'kw=', kwargs)
                    plan.build_tasks(CBSComputer, **packet, nlevel=nlevel, level=n, **cbsmeta, **kwargs)

                else:
                    print('PLANNING NBody(FinDif(CBS)', 'n=', n, nlevel, 'packet=', packet, 'cbsmeta=', cbsmeta ,'findif_kw=', current_findif_kwargs, 'kw=', kwargs)
                    plan.build_tasks(FinDifComputer, **packet, nlevel=nlevel, level=n, findif_mode=dermode, computer=CBSComputer, **cbsmeta, **current_findif_kwargs, **kwargs)

            else:
                dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
                if dermode[0] == dermode[1]:  # analytic
                    print('PLANNING NBody', 'n=', n, nlevel, 'packet=', packet)
                    plan.build_tasks(SingleComputer, **packet, nlevel=nlevel, level=n)
                else:
                    print('PLANNING NBody(FinDif)', 'n=', n, nlevel, 'packet=', packet, 'findif_kw=', current_findif_kwargs, 'kw=', kwargs)
                    plan.build_tasks(FinDifComputer, **packet, nlevel=nlevel, level=n, findif_mode=dermode, **current_findif_kwargs, **kwargs)

        return plan

    # Check for CBS
    elif method == "cbs":
        kwargs.update(cbsmeta)
        print('PLANNING CBS', 'keywords=', keywords)
        plan = CBSComputer(**packet, **kwargs)

        methods = [sr.method for sr in plan.task_list]
        # TODO: pass more info, so fn can use for managed_methods -- ref, qc_module, fc/ae, conv/df
        dermode = negotiate_derivative_type(driver, methods, kwargs.pop('dertype', None), verbose=1)

        if dermode[0] == dermode[1]:  # analytic
            return plan
        else:
            # For FD(CBS(Single)), the CBSComputer above is discarded after being used for dermode.
            print('PLANNING FinDif(CBS)', 'dermode=', dermode, 'keywords=', keywords, 'findif_kw=', current_findif_kwargs)
            plan = FinDifComputer(**packet, findif_mode=dermode, computer=CBSComputer, **current_findif_kwargs, **kwargs)
            return plan

    # Done with Wrappers -- know we want E, G, or H -- but may still be FinDif or SingleComputer
    else:
        dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
        convcrit = negotiate_convergence_criterion(dermode, method, return_optstash=False)

        if dermode[0] == dermode[1]:  # analytic
            print('PLANNING Single', 'keywords=', keywords)
            return SingleComputer(**packet, **kwargs)
        else:
            keywords.update(convcrit)
            print('PLANNING FinDif', 'dermode=', dermode, 'keywords=', keywords, 'findif_kw=', current_findif_kwargs, 'kw=', kwargs)
            return FinDifComputer(**packet, findif_mode=dermode, **current_findif_kwargs, **kwargs)
