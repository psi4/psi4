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

__all__ = [
    "AtomicComputer",
    "BaseComputer",
    "EnergyGradientHessianWfnReturn",
]

import abc
import copy
import logging
import pprint
from typing import Any, Dict, Optional, Tuple, Union, TYPE_CHECKING

try:
    from pydantic.v1 import Field, validator
except ImportError:
    from pydantic import Field, validator
import qcelemental as qcel
from qcelemental.models import DriverEnum, AtomicInput, AtomicResult
qcel.models.molecule.GEOMETRY_NOISE = 13  # need more precision in geometries for high-res findif
import qcengine as qcng

from psi4 import core
from . import p4util

if TYPE_CHECKING:
    import qcportal

logger = logging.getLogger(__name__)

EnergyGradientHessianWfnReturn = Union[float, core.Matrix, Tuple[Union[float, core.Matrix], core.Wavefunction]]


class BaseComputer(qcel.models.ProtoModel):
    """Base class for "computers" that plan, run, and process QC tasks."""

    @abc.abstractmethod
    def compute(self):
        pass

    @abc.abstractmethod
    def plan(self):
        pass

    class Config(qcel.models.ProtoModel.Config):
        extra = "allow"
        allow_mutation = True


class AtomicComputer(BaseComputer):
    """Computer for analytic single-geometry computations."""

    molecule: Any = Field(..., description="The molecule to use in the computation.")
    basis: str = Field(..., description="The quantum chemistry basis set to evaluate (e.g., 6-31g, cc-pVDZ, ...).")
    method: str = Field(..., description="The quantum chemistry method to evaluate (e.g., B3LYP, MP2, ...).")
    driver: DriverEnum = Field(..., description="The resulting type of computation: energy, gradient, hessian, properties."
        "Note for finite difference that this should be the target driver, not the means driver.")
    keywords: Dict[str, Any] = Field(default_factory=dict, description="The keywords to use in the computation.")
    computed: bool = Field(False, description="Whether quantum chemistry has been run on this task.")
    result: Any = Field(default_factory=dict, description=":py:class:`~qcelemental.models.AtomicResult` return.")
    result_id: Optional[str] = Field(None, description="The optional ID for the computation.")

    class Config(qcel.models.ProtoModel.Config):
        pass

    @validator("basis")
    def set_basis(cls, basis):
        return basis.lower()

    @validator("method")
    def set_method(cls, method):
        return method.lower()

    @validator("keywords")
    def set_keywords(cls, keywords):
        return copy.deepcopy(keywords)

    def plan(self) -> AtomicInput:
        """Form QCSchema input from member data."""

        atomic_model = AtomicInput(**{
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
                "wfn_qcvars_only": True,
            },
        })

        return atomic_model

    def compute(self, client: Optional["qcportal.client.FractalClient"] = None):
        """Run quantum chemistry."""
        from psi4.driver import pp

        if self.computed:
            return

        if client:
            self.computed = True

            try:
                # QCFractal v0.15.8
                from qcportal.models import KeywordSet, Molecule
                qca_next_branch = False
            except ImportError:
                # QCFractal `next`
                from qcelemental.models import Molecule
                qca_next_branch = True

            # Build the molecule
            mol = Molecule(**self.molecule.to_schema(dtype=2))

            if not qca_next_branch:
                # QCFractal v0.15.8

                # Build the keywords
                keyword_id = client.add_keywords([KeywordSet(values=self.keywords)])[0]

                r = client.add_compute("psi4", self.method, self.basis, self.driver, keyword_id, [mol])
                self.result_id = r.ids[0]
                # NOTE: The following will re-run errored jobs by default
                if self.result_id in r.existing:
                    ret = client.query_tasks(base_result=self.result_id)
                    if ret:
                        if ret[0].status == "ERROR":
                            client.modify_tasks("restart", base_result=self.result_id)
                            logger.info("Resubmitting Errored Job {}".format(self.result_id))
                        elif ret[0].status == "COMPLETE":
                            logger.debug("Job already completed {}".format(self.result_id))
                    else:
                        logger.debug("Job already completed {}".format(self.result_id))
                else:
                    logger.debug("Submitting AtomicResult {}".format(self.result_id))

            else:
                # QCFractal `next`

                meta, ids = client.add_singlepoints(
                    molecules=mol,
                    program="psi4",
                    driver=self.driver,
                    method=self.method,
                    basis=self.basis,
                    keywords=self.keywords,
                    # protocols,
                )
                self.result_id = ids[0]
                # NOTE: The following will re-run errored jobs by default
                if meta.existing_idx:
                    rec = client.get_singlepoints(self.result_id)
                    if rec.status == "error":
                        client.reset_records(self.result_id)
                        logger.info("Resubmitting Errored Job {}".format(self.result_id))
                    elif rec.status == "complete":
                        logger.debug("Job already completed {}".format(self.result_id))
                else:
                    logger.debug("Submitting AtomicResult {}".format(self.result_id))

            return

        logger.info(f'<<< JSON launch ... {self.molecule.schoenflies_symbol()} {self.molecule.nuclear_repulsion_energy()}')
        gof = core.get_output_file()

        # EITHER ...
        # from psi4.driver import schema_wrapper
        # self.result = schema_wrapper.run_qcschema(self.plan())
        # ... OR ...
        self.result = qcng.compute(
            self.plan(),
            "psi4",
            raise_error=True,
            # local_options below suitable for serial mode where each job takes all the resources of the parent Psi4 job.
            #   distributed runs through QCFractal will likely need a different setup.
            task_config={
                # B -> GiB
                "memory": core.get_memory() / (2 ** 30),
                "ncores": core.get_num_threads(),
            },
        )
        # ... END

        #pp.pprint(self.result.dict())
        #print("... JSON returns >>>")
        core.set_output_file(gof, True)
        core.reopen_outfile()
        logger.debug(pp.pformat(self.result.dict()))
        core.print_out(_drink_filter(self.result.dict()["stdout"]))
        self.computed = True

    def get_results(self, client: Optional["qcportal.FractalClient"] = None) -> AtomicResult:
        """Return results as Atomic-flavored QCSchema."""

        if self.result:
            return self.result

        if client:
            try:
                # QCFractal/QCPortal v0.15.8
                result = client.query_results(id=self.result_id)
                qca_next_branch = False
            except AttributeError:
                # QCFractal/QCPortal `next`
                record = client.get_singlepoints(record_ids=self.result_id)
                qca_next_branch = True

            logger.debug(f"Querying AtomicResult {self.result_id}")

            if not qca_next_branch:
                # QCFractal v0.15.8
                if len(result) == 0:
                    return self.result

                self.result = result[0]

            else:
                # QCFractal `next`
                if record.status != "complete":
                    return self.result

                self.result = _singlepointrecord_to_atomicresult(record)

            return self.result


def _singlepointrecord_to_atomicresult(spr: "qcportal.singlepoint.SinglepointRecord") -> AtomicResult:
    atres = spr.to_qcschema_result()

    # QCFractal `next` database stores return_result, properties, and extras["qcvars"] merged
    #   together and with lowercase keys. `to_qcschema_result` partitions properties back out,
    #   but we need to restore qcvars keys, types, and dimensions.
    qcvars = atres.extras.pop("extra_properties")
    qcvars.pop("return_result")
    qcvars = {k.upper(): p4util.plump_qcvar(k, v) for k, v in qcvars.items()}
    atres.extras["qcvars"] = qcvars

    return atres


def _drink_filter(stdout: str) -> str:
    """Don't mess up the widespread ``grep beer`` test of Psi4 doneness by printing multiple drinks per outfile."""

    stdout = stdout.replace("\n*** Psi4 exiting successfully. Buy a developer a beer!", "")
    stdout = stdout.replace("\n*** Psi4 encountered an error. Buy a developer more coffee!", "")
    return stdout
