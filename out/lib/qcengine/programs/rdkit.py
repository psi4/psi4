"""
Calls the RDKit package.
"""

from typing import TYPE_CHECKING, Dict

from qcelemental.models import AtomicResult, Provenance
from qcelemental.util import safe_version, which_import

from ..exceptions import InputError
from ..units import ureg
from .model import ProgramHarness

if TYPE_CHECKING:
    from qcelemental.models import AtomicInput

    from ..config import TaskConfig


class RDKitHarness(ProgramHarness):

    _defaults = {
        "name": "RDKit",
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
    def _process_molecule_rdkit(jmol):
        from rdkit import Chem

        # Handle errors
        if abs(jmol.molecular_charge) > 1.0e-6:
            raise InputError("RDKit does not currently support charged molecules.")

        if not jmol.connectivity:  # Check for empty list
            raise InputError("RDKit requires molecules to have a connectivity graph.")

        # Build out the base molecule
        base_mol = Chem.Mol()
        rw_mol = Chem.RWMol(base_mol)
        for sym in jmol.symbols:
            rw_mol.AddAtom(Chem.Atom(sym.title()))

        # Add in connectivity
        bond_types = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
        for atom1, atom2, bo in jmol.connectivity:
            rw_mol.AddBond(atom1, atom2, bond_types[bo])

        mol = rw_mol.GetMol()

        # Write out the conformer
        natom = len(jmol.symbols)
        conf = Chem.Conformer(natom)
        bohr2ang = ureg.conversion_factor("bohr", "angstrom")
        for line in range(natom):
            conf.SetAtomPosition(
                line,
                (
                    bohr2ang * jmol.geometry[line, 0],
                    bohr2ang * jmol.geometry[line, 1],
                    bohr2ang * jmol.geometry[line, 2],
                ),
            )

        mol.AddConformer(conf)
        Chem.rdmolops.SanitizeMol(mol)

        return mol

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which_import(
            "rdkit",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install rdkit -c conda-forge`.",
        )

    def get_version(self) -> str:
        """Return the currently used version of RDKit."""
        self.found(raise_error=True)

        which_prog = which_import("rdkit")
        if which_prog not in self.version_cache:
            import rdkit

            self.version_cache[which_prog] = safe_version(rdkit.__version__)

        return self.version_cache[which_prog]

    def compute(self, input_data: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        """
        Runs RDKit in FF typing
        """

        self.found(raise_error=True)
        import rdkit
        from rdkit.Chem import AllChem

        # Failure flag
        ret_data = {"success": False}

        # Build the Molecule
        jmol = input_data.molecule
        mol = self._process_molecule_rdkit(jmol)

        if input_data.model.method.lower() == "uff":
            ff = AllChem.UFFGetMoleculeForceField(mol)
            all_params = AllChem.UFFHasAllMoleculeParams(mol)
        elif input_data.model.method.lower() in ["mmff94", "mmff94s"]:
            props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=input_data.model.method)
            ff = AllChem.MMFFGetMoleculeForceField(mol, props)
            all_params = AllChem.MMFFHasAllMoleculeParams(mol)
        else:
            raise InputError("RDKit only supports the UFF, MMFF94, and MMFF94s methods currently.")

        if all_params is False:
            raise InputError("RDKit parameters not found for all atom types in molecule.")

        ff.Initialize()

        ret_data["properties"] = {"return_energy": ff.CalcEnergy() * ureg.conversion_factor("kJ / mol", "hartree")}

        if input_data.driver == "energy":
            ret_data["return_result"] = ret_data["properties"]["return_energy"]
        elif input_data.driver == "gradient":
            coef = ureg.conversion_factor("kJ / mol", "hartree") * ureg.conversion_factor("angstrom", "bohr")
            ret_data["return_result"] = [x * coef for x in ff.CalcGrad()]
        else:
            raise InputError(f"Driver {input_model.driver} not implemented for RDKit.")

        ret_data["provenance"] = Provenance(
            creator="rdkit", version=rdkit.__version__, routine="rdkit.Chem.AllChem.UFFGetMoleculeForceField"
        )

        ret_data["schema_name"] = "qcschema_output"
        ret_data["success"] = True

        # Form up a dict first, then sent to BaseModel to avoid repeat kwargs which don't override each other
        return AtomicResult(**{**input_data.dict(), **ret_data})
