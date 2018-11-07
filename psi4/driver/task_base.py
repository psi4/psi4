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

# from psi4 import core
# from psi4.driver import p4util
# from psi4.driver import constants
# from psi4.driver.p4util.exceptions import *


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
