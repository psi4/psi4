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

from __future__ import print_function
from __future__ import absolute_import
import abc
import math
import itertools
import pydantic

from typing import Dict, List, Any, Union

import numpy as np

from psi4 import core
from psi4.driver import p4util


# __all__ = ["BaseTask", "SingleResult", "planner"]


class BaseTask(pydantic.BaseModel, abc.ABC):
    @abc.abstractmethod
    def compute(self):
        pass

    @abc.abstractmethod
    def plan(self):
        pass


class SingleResult(BaseTask):

    molecule: Any
    basis: str
    method: str
    driver: str
    keywords: Dict[str, Any] = None
    computed: bool = False
    result: Dict[str, Any] = None

    def plan(self):

        data = {
            "schema_name": "qc_schema_input",
            "schema_version": 1,
            "molecule": self.molecule.to_schema(dtype=1)["molecule"],
            "driver": self.driver,
            "model": {
                "method": self.method,
                "basis": self.basis
            },
            "keywords": self.keywords,
        }

        return data

    def compute(self):
        if self.computed:
            return

        print(json.dumps(self.plan(), indent=2))
        from psi4.driver import json_wrapper
        self.result = json_wrapper.run_json(self.plan())
        self.computed = True

    def get_results(self):
        return self.result


def planner(driver, method, molecule, **kwargs):
    """Plans a task graph of a complex computations.

    Canonical Task layering:
     - N-Body
     - CBS
     - Finite Difference
     - Single Result

    Parameters
    ----------
    driver :
        The resulting type of computation {"energy", "gradient", "hessian"}
    method :
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
    from psi4.driver.driver_nbody import NBodyComputer, nbody_gufunc
    from psi4.driver.driver_cbs import CBSComputer, cbs_gufunc, _cbs_text_parser

    # Only pull the changed options
    # keywords = {key: v['value'] for key, v in p4util.prepare_options_for_modules(changedOnly=True)['GLOBALS'].items()}
    keywords = {}
    if not len(keywords):
        keywords = None
        basis = None
    else:
        basis = keywords.pop("BASIS", None)

    # Pull basis out of kwargs, override globals if user specified
    basis = kwargs.pop("basis", basis)

    # Expand CBS methods
    if "/" in method:
        tmp_kwargs = kwargs.copy()
        tmp_kwargs["ptype"] = driver
        cbsmeta = _cbs_text_parser(method, **tmp_kwargs)

        # Single call detected
        if "cbs_metadata" not in cbsmeta:
            method = cbsmeta["method"]
            basis = cbsmeta["basis"]
        else:
            method = "cbs"
    else:
        cbsmeta = {}


    # Build a packet
    packet = {"molecule": molecule, "driver": driver, "method": method, "basis": basis, "keywords": keywords}

    # First check for BSSE type
    if 'bsse_type' in kwargs:
        raise Exception("BSSE not yet ready")
        NBodyComputer

    # Check for CBS
    elif method.lower() == "cbs":

        return CBSComputer(**packet, **cbsmeta)

    # Finally a single result
    else:

        return SingleResult(**packet, **kwargs)

