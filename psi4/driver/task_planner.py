#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2018 The Psi4 Developers.
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
import math
import itertools
import pydantic

from typing import Dict, List, Any, Union

import numpy as np

from psi4.driver import p4util
from psi4.driver.p4util.exceptions import UpgradeHelper
from psi4.driver.task_base import SingleResult
from psi4.driver.driver_findif import FinDifComputer
from psi4.driver.driver_nbody import NBodyComputer
from psi4.driver.driver_cbs import CBSComputer, _cbs_text_parser
from psi4.driver.driver_util import negotiate_derivative_type

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
     - CBS
     - FinDif
     - SingleResult

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
    print('task_planner()', keywords)

    try:
        method.lower()
    except AttributeError as e:
        if method.__name__ == 'cbs':
            raise UpgradeHelper(method, repr(method.__name__), 1.4, ' Replace cbs or complete_basis_set function with cbs string.')
        else:
            raise e

    # Pull basis out of kwargs, override globals if user specified
    basis = keywords.pop("BASIS", None)
    basis = kwargs.pop("basis", basis)
    method = method.lower()

    # Expand CBS methods
    method, basis, cbsmeta = _expand_cbs_methods(method, basis, driver, **kwargs)

    # Build a packet
    packet = {"molecule": molecule, "driver": driver, "method": method, "basis": basis, "keywords": keywords}

    # First check for BSSE type
    if 'bsse_type' in kwargs:
        plan = NBodyComputer(**packet, **kwargs)
        del packet["molecule"]

        # Add tasks for every nbody level requested
        levels = kwargs.pop('levels', {plan.max_nbody: method})
        plan.max_nbody = max([i for i in levels if isinstance(i, int)])

        for n, level in levels.items():
            method, basis, cbsmeta = _expand_cbs_methods(level, basis, driver, cbsmeta=cbsmeta, **kwargs)
            packet.update({'method': method, 'basis': basis})

            if method == "cbs":
                print('PLANNING NBody(CBS)')
                plan.build_tasks(CBSComputer, **packet, **cbsmeta, max_nbody=n)
            else:
                print('PLANNING NBody()')
                plan.build_tasks(SingleResult, **packet, max_nbody=n)

        return plan

    # Check for CBS
    elif method == "cbs":
        kwargs.update(cbsmeta)
        print('PLANNING CBS')
        return CBSComputer(**packet, **kwargs)

    # Done with Wrappers -- know we want E, G, or H -- but may still be FinDif or SingleResult
    else:
        dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1, return_strategy=True)

        if dermode == 'analytic':
            print('PLANNING ()')
            return SingleResult(**packet, **kwargs)
        else:
            print('PLANNING FinDif', dermode)
            return FinDifComputer(**packet, findif_mode=dermode, **kwargs)
