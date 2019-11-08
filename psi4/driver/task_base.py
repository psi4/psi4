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
import math
import json
import pprint
pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)
from typing import Any, Dict, List, Optional, Union
import itertools

import numpy as np
import pydantic
import qcelemental as qcel
from qcelemental.models import DriverEnum, ResultInput
qcel.models.molecule.GEOMETRY_NOISE = 13  # need more precision in geometries for high-res findif
import qcengine as qcng

from psi4 import core
from psi4.driver import p4util
from psi4.driver.p4util import exceptions

__all__ = ["BaseComputer", "SingleComputer"]


class BaseComputer(qcel.models.ProtoModel):
    @abc.abstractmethod
    def compute(self):
        pass

    @abc.abstractmethod
    def plan(self):
        pass

    class Config(qcel.models.ProtoModel.Config):
        extra = 'allow'
        allow_mutation = True


class SingleComputer(BaseComputer):

    molecule: Any
    basis: str
    method: str
    driver: DriverEnum
    keywords: Dict[str, Any] = {}
    computed: bool = False
    result: Any = {}

    result_id: str = None

    class Config(qcel.models.ProtoModel.Config):
        pass

    @pydantic.validator('basis')
    def set_basis(cls, basis):
        return basis.lower()

    @pydantic.validator('method')
    def set_method(cls, method):
        return method.lower()

    def plan(self):

        data = ResultInput(**{
            "molecule": self.molecule.to_schema(dtype=2),
            "driver": self.driver,
            "model": {
                "method": self.method,
                "basis": self.basis
            },
            "keywords": self.keywords,
            "protocols": {
                "stdout": True,
            },
        })

        return data

    def compute(self, client=None):

        if self.computed:
            return

        if client:

            self.computed = True
            from qcfractal.interface.models import KeywordSet
            from qcfractal.interface import Molecule

            # Build the keywords
            keyword_id = client.add_keywords([KeywordSet(values=self.keywords)])[0]

            # Build the molecule
            mol = Molecule(**self.molecule.to_schema(dtype=2))

            r = client.add_compute("psi4", self.method, self.basis, self.driver, keyword_id, [mol])
            self.result_id = r.ids[0]
            # NOTE: The following will re-run errored jobs by default
            if self.result_id in r.existing:
                ret = client.query_tasks(base_result=self.result_id)
                if ret:
                    if ret[0].status == "ERROR":
                        upd = client.modify_tasks("restart",base_result=self.result_id)
                        print("Resubmitting Errored Job {}".format(self.result_id))
                    elif ret[0].status == "COMPLETE":
                        print("Job already completed {}".format(self.result_id))
                else:
                    print("Job already completed {}".format(self.result_id))
            else:
                print("Submitting Single Result {}".format(self.result_id))

            return

        # gof = core.get_output_file()
        # core.close_outfile()

        print('<<< JSON launch ...', self.molecule.schoenflies_symbol(), self.molecule.nuclear_repulsion_energy())
        #pp.pprint(self.plan().dict())

        # EITHER ...
        #from psi4.driver import schema_wrapper
        #self.result = schema_wrapper.run_qcschema(self.plan())
        # ... OR ...
        self.result = qcng.compute(self.plan(), 'psi4', raise_error=True,
                                   # local_options below suitable for continuous mode
                                   local_options={"memory": core.get_memory()/1073741824, "ncores": core.get_num_threads()}
                                  )
        # ... END

        #pp.pprint(self.result.dict())
        #print('... JSON returns >>>')
        self.computed = True

        # core.set_output_file(gof, True)

    def get_results(self, client=None):
        if self.result:
            return self.result

        if client:
            result = client.query_results(id=self.result_id)
            print("Querying Single Result {}".format(self.result_id))
            if len(result) == 0:
                return self.result

            self.result = result[0].dict(encoding='msgpack-ext')
            return self.result

    def get_json_results(self):
        return self.result
