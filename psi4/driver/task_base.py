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
import json
import itertools
import pydantic

from typing import Dict, List, Any, Union

import numpy as np

from psi4 import core
from psi4.driver import p4util
from psi4.driver.p4util import exceptions

__all__ = ["BaseTask", "SingleResult"]


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
    computed: bool = False
    result: Dict[str, Any] = None

    @pydantic.validator('basis')
    def set_basis(cls, basis):
        return basis.lower()

    @pydantic.validator('method')
    def set_method(cls, method):
        return method.lower()

    @pydantic.validator('driver')
    def set_driver(cls, driver):
        driver = driver.lower()
        if driver not in ["energy", "gradient", "hessian"]:
            raise exceptions.ValidationError(
                f"Driver must be either energy, gradient, or hessian. Found {driver}.")

        return driver

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

    def get_json_results(self):
        return self.result


# use from qcel once settled
def unnp(dicary, flat=False, _path=None):
    """Return `dicary` with any ndarray values replaced by lists.

    Parameters
    ----------
    dicary: dict
        Dictionary where any internal iterables are dict or list.
    flat : bool, optional
        Whether the returned lists are flat or nested.

    Returns
    -------
    dict
        Input with any ndarray values replaced by lists.

    """
    if _path is None:
        _path = []

    ndicary = {}
    for k, v in dicary.items():
        if isinstance(v, dict):
            ndicary[k] = unnp(v, flat, _path + [str(k)])
        elif isinstance(v, list):
            # relying on Py3.6+ ordered dict here
            fakedict = {kk: vv for kk, vv in enumerate(v)}
            tolisted = unnp(fakedict, flat, _path + [str(k)])
            ndicary[k] = list(tolisted.values())
        else:
            try:
                v.shape
            except AttributeError:
                ndicary[k] = v
            else:
                if flat:
                    ndicary[k] = v.ravel().tolist()
                else:
                    ndicary[k] = v.tolist()
    return ndicary


def plump_qcvar(val, shape_clue, ret='np'):
    """Convert flat arra

    Parameters
    ----------
    val : list or scalar
        flat (?, ) list or scalar, probably from JSON storage.
    shape_clue : str
        Label that includes (case insensitive) one of the following as
        a clue to the array's natural dimensions: 'gradient', 'hessian'
    ret : {'np', 'psi4'}
        Whether to return `np.ndarray` or `psi4.core.Matrix`.

    Returns
    -------
    np.ndarray or psi4.core.Matrix
        Reshaped array of type `ret` with natural dimensions of `shape_clue`.

    Raises
    ------
    TODO

    """
    if isinstance(val, (np.ndarray, core.Matrix)):
        raise TypeError
    elif isinstance(val, list):
        tgt = np.asarray(val)
    else:
        # presumably scalar
        return val

    if 'gradient' in shape_clue.lower():
        reshaper = (-1, 3)
    elif 'hessian' in shape_clue.lower():
        ndof = int(math.sqrt(len(tgt)))
        reshaper = (ndof, ndof)
    else:
        raise ValidationError(f'Uncertain how to reshape array: {shape_clue}')

    if ret == 'np':
        return tgt.reshape(reshaper)
    elif ret == 'psi4':
        return core.Matrix.from_array(tgt.reshape(reshaper))
#wfn.gradient().np.ravel().tolist()
    else:
        raise ValidationError(f'Return type not among [np, psi4]: {ret}')
