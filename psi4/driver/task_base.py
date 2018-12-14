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
    keywords: Dict[str, Any] = {}

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
        from psi4.driver import json_wrapper

        return json_wrapper.run_json(self.plan())


def planner(name, **kwargs):

    keywords = {i: j['value'] for i, j in p4util.prepare_options_for_modules(changedOnly=True)['GLOBALS'].items()}
    data = {'driver': kwargs['ptype'], 'method': name, 'basis': core.get_global_option('BASIS'), 'keywords': keywords}

    comp_plan = {}
    if 'bsse_type' in kwargs:
        # Call nbody wrapper
        from psi4.driver.driver_nbody import nbody_gufunc
        from psi4.driver.driver_nbody import NBodyComputer

        comp_plan.update({'function': nbody_gufunc, 'computer': SingleResult, 'data': data})
        ComputeInstance = NBodyComputer

    if hasattr(name, '__call__') and name.__name__ in ['cbs', 'complete_basis_set']:
        function = comp_plan.get('function', name)

        if function != name:
            # Use CBSComputer inside the nbody wrapper
            from psi4.driver.driver_cbs import CBSComputer

            data.update({k: v for k, v in kwargs.items() if k not in ComputeInstance.__fields__})
            comp_plan.update({'computer': CBSComputer, 'data': data})

        else:
            # Call CBS wrapper
            comp_plan.update({'function': function})
            comp_plan.update({'name': kwargs.pop('label', 'custom function')})

    elif '/' in name:
        from psi4.driver.driver_cbs import _cbs_text_parser
        from psi4.driver.driver_cbs import _cbs_gufunc
        from psi4.driver.driver_cbs import CBSComputer

        tmp = _cbs_text_parser(name, **kwargs)

        function = comp_plan.get('function', _cbs_gufunc)

        if function != _cbs_gufunc:
            if 'cbs_metadata' in tmp:
                # Use CBSComputer inside the nbody wrapper
                data.update({'cbs_metadata': tmp['cbs_metadata']})
                comp_plan.update({'computer': CBSComputer, 'data': data})

            else:
                # Use SingleResult inside the nbody wrapper
                data.update(tmp)

        else:
            # Call _cbs_gufunc
            comp_plan.update({'function': function})


    return comp_plan

