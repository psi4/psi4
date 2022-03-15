"""
Calls the TorchANI package.
"""

from typing import TYPE_CHECKING, Dict

from qcelemental.models import AtomicResult, Provenance
from qcelemental.util import parse_version, safe_version, which_import

from ..exceptions import InputError, ResourceError
from ..units import ureg
from .model import ProgramHarness

if TYPE_CHECKING:
    from qcelemental.models import AtomicInput

    from ..config import TaskConfig


class TorchANIHarness(ProgramHarness):

    _CACHE = {}

    _defaults = {
        "name": "TorchANI",
        "scratch": False,
        "thread_safe": True,
        "thread_parallel": False,
        "node_parallel": False,
        "managed_memory": False,
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which_import(
            "torchani",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `pip install torchani`.",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which_import("torchani")
        if which_prog not in self.version_cache:
            import torchani

            self.version_cache[which_prog] = safe_version(torchani.__version__)

        return self.version_cache[which_prog]

    def get_model(self, name: str) -> "torchani.models.BuiltinModels":
        name = name.lower()

        if name in self._CACHE:
            return self._CACHE[name]

        import torch
        import torchani

        # graft a custom forward pass to Ensemble class
        def ensemble_forward(self, species_input):
            outputs = torch.cat([x(species_input)[1] for x in self])
            species, _ = species_input
            return species, outputs

        torchani.nn.Ensemble.forward = ensemble_forward

        ani_models = {
            "ani1x": torchani.models.ANI1x,
            "ani1ccx": torchani.models.ANI1ccx,
        }

        if parse_version(self.get_version()) >= parse_version("2.0"):
            ani_models["ani2x"] = torchani.models.ANI2x

        try:
            self._CACHE[name] = ani_models[name]()
        except KeyError:
            raise InputError(f"TorchANI only accepts methods: {ani_models.keys()}")

        return self._CACHE[name]

    def compute(self, input_data: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        """
        Runs TorchANI in FF typing
        """

        # Check if existings and version
        self.found(raise_error=True)
        if parse_version(self.get_version()) < parse_version("0.9"):
            raise ResourceError("QCEngine's TorchANI wrapper requires version 0.9 or greater.")

        import numpy as np
        import torch
        import torchani

        device = torch.device("cpu")

        # Failure flag
        ret_data = {"success": False}

        # Build model
        method = input_data.model.method
        model = self.get_model(method)

        # Build species
        species = input_data.molecule.symbols

        known_sym = {"H", "C", "N", "O"}
        if method.lower() == "ani2x":
            known_sym.update({"S", "F", "Cl"})

        unknown_sym = set(species) - known_sym
        if unknown_sym:
            raise InputError(f"TorchANI model '{method}' does not support symbols: {unknown_sym}.")

        num_atoms = len(species)
        species = model.species_to_tensor(species).to(device).unsqueeze(0)

        # Build coord array
        geom_array = input_data.molecule.geometry.reshape(1, -1, 3) * ureg.conversion_factor("bohr", "angstrom")
        coordinates = torch.tensor(geom_array.tolist(), requires_grad=True, device=device)

        _, energy_array = model((species, coordinates))
        energy = energy_array.mean()
        ensemble_std = energy_array.std()
        ensemble_scaled_std = ensemble_std / np.sqrt(num_atoms)

        ret_data["properties"] = {"return_energy": energy.item()}

        if input_data.driver == "energy":
            ret_data["return_result"] = ret_data["properties"]["return_energy"]
        elif input_data.driver == "gradient":
            derivative = torch.autograd.grad(energy.sum(), coordinates)[0].squeeze()
            ret_data["return_result"] = (
                np.asarray(derivative * ureg.conversion_factor("angstrom", "bohr")).ravel().tolist()
            )
        elif input_data.driver == "hessian":
            hessian = torchani.utils.hessian(coordinates, energies=energy)
            ret_data["return_result"] = np.asarray(hessian)
        else:
            raise InputError(
                f"TorchANI can only compute energy, gradient, and hessian driver methods. Found {input_data.driver}."
            )

        #######################################################################
        # Description of the quantities stored in `extras`
        #
        # ensemble_energies:
        #   An energy array of all members (models) in an ensemble of models
        #
        # ensemble_energy_avg:
        #   The average value of energy array which is also recorded with as
        #   `energy` in QCEngine
        #
        # ensemble_energy_std:
        #   The standard deviation of energy array
        #
        # ensemble_per_root_atom_disagreement:
        #   The standard deviation scaled by the square root of N, with N being
        #   the number of atoms in the molecule. This is the quantity used in
        #   the query-by-committee (QBC) process in active learning to infer
        #   the reliability of the models in an ensemble, and produce more data
        #   points in the regions where this quantity is below a certain
        #   threshold (inclusion criteria)
        ret_data["extras"] = input_data.extras.copy()
        ret_data["extras"].update(
            {
                "ensemble_energies": energy_array.detach().numpy(),
                "ensemble_energy_avg": energy.item(),
                "ensemble_energy_std": ensemble_std.item(),
                "ensemble_per_root_atom_disagreement": ensemble_scaled_std.item(),
            }
        )

        ret_data["provenance"] = Provenance(
            creator="torchani", version="unknown", routine="torchani.builtin.aev_computer"
        )

        ret_data["schema_name"] = "qcschema_output"
        ret_data["success"] = True

        # Form up a dict first, then sent to BaseModel to avoid repeat kwargs which don't override each other
        return AtomicResult(**{**input_data.dict(), **ret_data})
