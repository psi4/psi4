"""
The qcore QCEngine Harness
"""

from typing import TYPE_CHECKING, Any, Dict, List, Set

import numpy as np
from qcelemental.models import AtomicResult, BasisSet
from qcelemental.util import parse_version, safe_version, which_import

from ..exceptions import InputError, UnknownError
from .model import ProgramHarness
from .util import (
    cca_ao_order_spherical,
    get_ao_conversion,
    reorder_column_ao_indices,
    reorder_row_and_column_ao_indices,
)

if TYPE_CHECKING:
    from qcelemental.models import AtomicInput

    from ..config import TaskConfig


def qcore_ao_order_spherical(max_angular_momentum: int) -> Dict[int, List[int]]:
    ao_order = {}
    for ang_mom in range(max_angular_momentum):
        ao_order[ang_mom] = [x for x in range(ang_mom, -1 * ang_mom - 1, -1)]
    return ao_order


class QcoreHarness(ProgramHarness):
    _defaults: Dict[str, Any] = {
        "name": "Qcore",
        "scratch": False,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": False,
    }
    version_cache: Dict[str, str] = {}

    # List of DFT functionals
    _dft_functionals: Set[str] = {
        "SLATER",
        "DIRAC",
        "SLATERD3",
        "DIRACD3",
        "VWN5",
        "VWN",
        "VWN1",
        "SVWN",
        "LDA",
        "BLYP",
        "BPW91",
        "BLYPD3",
        "B88",
        "PBEX",
        "PBERX",
        "PBEC",
        "LYP",
        "PW91",
        "P86",
        "PBE",
        "PBER",
        "PBED3",
        "PBERD3",
        "B3LYP3",
        "B3LYP",
        "B3LYP5",
        "PBE0",
        "PBE1PBE",
        "B3LYP3D3",
        "B3LYPD3",
        "B3LYP5D3",
        "PBE0D3",
        "PBE1PBED3",
        "CAMB3LYP",
        "WB97X",
        "CAMB3LYPD3",
        "WB97XD3",
    }

    _xtb_models: Set[str] = {"GFN1", "GFN0"}

    # This map order converts qcore ordering to CCA ordering
    # Entos spherical basis ordering for each angular momentum. Follows reverse order of CCA.
    _qcore_to_cca_ao_order = {"spherical": get_ao_conversion(cca_ao_order_spherical(10), qcore_ao_order_spherical(10))}

    class Config(ProgramHarness.Config):
        pass

    def found(self, raise_error: bool = False) -> bool:
        return which_import(
            "qcore",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install py-qcore -c entos -c conda-forge`.",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which_import("qcore")
        if which_prog not in self.version_cache:
            import qcore

            self.version_cache[which_prog] = safe_version(qcore.__version__)

        return self.version_cache[which_prog]

    def compute(self, input_model: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        """
        Run qcore
        """
        # Check if qcore executable is found
        self.found(raise_error=True)

        # Check qcore version
        if parse_version(self.get_version()) < parse_version("0.8.9"):
            raise TypeError(f"qcore version {self.get_version()} not supported")

        import qcore

        if isinstance(input_model.model.basis, BasisSet):
            raise InputError("QCSchema BasisSet for model.basis not implemented. Use string basis name.")

        method = input_model.model.method.upper()
        if method in self._dft_functionals:
            method = {"kind": "dft", "xc": method, "ao": input_model.model.basis}
        elif method == "HF":
            method = {"kind": "hf", "ao": input_model.model.basis}
        elif method in self._xtb_models:
            method = {"kind": "xtb", "model": method}
        else:
            raise InputError(f"Method is not valid: {method}")

        method["details"] = input_model.keywords

        qcore_input = {
            # "schema_name": "single_input",
            "molecule": {
                "geometry": input_model.molecule.geometry,
                "atomic_numbers": input_model.molecule.atomic_numbers,
                "charge": input_model.molecule.molecular_charge,
                "multiplicity": input_model.molecule.molecular_multiplicity,
            },
            "method": method,
            "result_contract": {"wavefunction": "all"},
            "result_type": input_model.driver,
        }
        try:
            result = qcore.run(qcore_input, ncores=config.ncores)
        except Exception as exc:
            return UnknownError(str(exc))

        return self.parse_output(result.dict(), input_model)

    def parse_output(self, output: Dict[str, Any], input_model: "AtomicInput") -> "AtomicResult":

        wavefunction_map = {
            "orbitals_alpha": "scf_orbitals_a",
            "orbitals_beta": "scf_orbitals_b",
            "density_alpha": "scf_density_a",
            "density_beta": "scf_density_b",
            "fock_alpha": "scf_fock_a",
            "fock_beta": "scf_fock_b",
            "eigenvalues_alpha": "scf_eigenvalues_a",
            "eigenvalues_beta": "scf_eigenvalues_b",
            "occupations_alpha": "scf_occupations_a",
            "occupations_beta": "scf_occupations_b",
        }

        output_data = input_model.dict()

        output_data["return_result"] = output[input_model.driver.value]

        # Always build a wavefunction, it will be stripped
        obas = output["wavefunction"]["ao_basis"]
        for k, center in obas["center_data"].items():
            # Convert basis set, cannot handle arrays
            for shell in center["electron_shells"]:
                shell.pop("normalized_primitives", None)
                for el_k in ["coefficients", "exponents", "angular_momentum"]:
                    shell[el_k] = shell[el_k].tolist()

            if center["ecp_potentials"] is not None:
                for shell in center["ecp_potentials"]:
                    shell.pop("ecp_potentials", None)
                    for ecp_k in ["angular_momentum", "r_exponents", "gaussian_exponents", "coefficients"]:
                        shell[ecp_k] = shell[ecp_k].tolist()

        basis_set = BasisSet(
            name=str(input_model.model.basis), center_data=obas["center_data"], atom_map=obas["atom_map"]
        )

        wavefunction = {"basis": basis_set}
        for key, qcschema_key in wavefunction_map.items():
            qcore_data = output["wavefunction"].get(key, None)
            if qcore_data is None:
                continue

            if ("density" in key) or ("fock" in key):
                qcore_data = reorder_row_and_column_ao_indices(qcore_data, basis_set, self._qcore_to_cca_ao_order)
            # Handles orbitals and 1D
            elif "orbitals" in key:
                qcore_data = reorder_column_ao_indices(qcore_data, basis_set, self._qcore_to_cca_ao_order)
            elif "eigenvalues" in key:
                qcore_data = reorder_column_ao_indices(
                    qcore_data.reshape(1, -1), basis_set, self._qcore_to_cca_ao_order
                ).ravel()

            elif "occupations" in key:
                tmp = np.zeros(basis_set.nbf)
                tmp[: qcore_data.shape[0]] = qcore_data
                qcore_data = reorder_column_ao_indices(
                    tmp.reshape(1, -1), basis_set, self._qcore_to_cca_ao_order
                ).ravel()
            else:
                raise KeyError("Wavefunction conversion key not understood")

            wavefunction[qcschema_key] = qcore_data

        wavefunction["restricted"] = True
        if "scf_eigenvalues_b" in wavefunction:
            wavefunction["restricted"] = False

        output_data["wavefunction"] = wavefunction

        # Handle remaining top level keys
        properties = {
            "calcinfo_nbasis": basis_set.nbf,
            "calcinfo_nmo": basis_set.nbf,
            "calcinfo_nalpha": np.sum(wavefunction["scf_occupations_a"] > 0),
            "calcinfo_natom": input_model.molecule.symbols.shape[0],
            "return_energy": output["energy"],
        }
        if wavefunction["restricted"]:
            properties["calcinfo_nbeta"] = properties["calcinfo_nalpha"]
        else:
            properties["calcinfo_nbeta"] = np.sum(wavefunction["scf_occupations_b"] > 0)

        output_data["properties"] = properties

        output_data["schema_name"] = "qcschema_output"
        output_data["success"] = True

        return AtomicResult(**output_data)


class EntosHarness(QcoreHarness):
    _defaults: Dict[str, Any] = {
        "name": "Entos",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
