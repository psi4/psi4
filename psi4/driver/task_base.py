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
import logging
from typing import Any, Dict

import pydantic
import qcelemental as qcel
from qcelemental.models import DriverEnum, AtomicInput
qcel.models.molecule.GEOMETRY_NOISE = 13  # need more precision in geometries for high-res findif
import qcengine as qcng

from psi4 import core

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

__all__ = ["BaseComputer", "AtomicComputer"]


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


class AtomicComputer(BaseComputer):

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

    @pydantic.validator('keywords')
    def set_keywords(cls, keywords):
        return copy.deepcopy(keywords)

    def plan(self):

        data = AtomicInput(**{
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
            "extras": {
                "psiapi": True,
            },
        })

        return data

    def compute(self, client=None):
        from psi4.driver import pp

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
                        client.modify_tasks("restart",base_result=self.result_id)
                        print("Resubmitting Errored Job {}".format(self.result_id))
                    elif ret[0].status == "COMPLETE":
                        print("Job already completed {}".format(self.result_id))
                else:
                    print("Job already completed {}".format(self.result_id))
            else:
                print("Submitting AtomicResult {}".format(self.result_id))

            return

        print('<<< JSON launch ...', self.molecule.schoenflies_symbol(), self.molecule.nuclear_repulsion_energy())
        logger.info(f'<<< JSON launch ... {self.molecule.schoenflies_symbol()} {self.molecule.nuclear_repulsion_energy()}')
        #pp.pprint(self.plan().dict())

        # EITHER ...
        #from psi4.driver import schema_wrapper
        #self.result = schema_wrapper.run_qcschema(self.plan())
        # ... OR ...
        self.result = qcng.compute(self.plan(), 'psi4', raise_error=True,
                                   # local_options below suitable for continuous mode
                                   local_options={"memory": core.get_memory()/1000000000, "ncores": core.get_num_threads()}
                                  )
        # ... END

        logger.debug(pp.pformat(self.result.dict()))
        #pp.pprint(self.result.dict())
        #print('... JSON returns >>>')
        self.computed = True

    def get_results(self, client=None):
        if self.result:
            return self.result

        if client:
            result = client.query_results(id=self.result_id)
            print("Querying AtomicResult {}".format(self.result_id))
            if len(result) == 0:
                return self.result

            self.result = result[0].dict(encoding='msgpack-ext')
            return self.result

    def get_json_results(self):
        return self.result
