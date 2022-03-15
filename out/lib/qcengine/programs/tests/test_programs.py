"""
Tests the DQM compute dispatch module
"""


import numpy as np
import pytest
from qcelemental.models import AtomicInput, Molecule

import qcengine as qcng
from qcengine.testing import failure_engine, using


def test_missing_key():
    ret = qcng.compute({"hello": "hi"}, "bleh")
    assert ret.success is False
    assert "hello" in ret.input_data


def test_missing_key_raises():
    with pytest.raises(qcng.exceptions.InputError):
        ret = qcng.compute({"hello": "hi"}, "bleh", raise_error=True)


@using("psi4")
def test_psi4_task():
    input_data = {
        "molecule": qcng.get_molecule("water"),
        "driver": "energy",
        "model": {"method": "SCF", "basis": "sto-3g"},
        "keywords": {"scf_type": "df"},
    }

    ret = qcng.compute(input_data, "psi4", raise_error=True)

    assert ret.driver == "energy"
    assert "Final Energy" in ret.stdout

    prov_keys = {"cpu", "hostname", "username", "wall_time"}
    assert ret.provenance.dict().keys() >= prov_keys
    assert "retries" not in ret.provenance.dict()

    assert ret.success is True


@using("psi4")
@using("gcp")
def test_psi4_hf3c_task():
    input_data = {
        "molecule": qcng.get_molecule("water"),
        "driver": "energy",
        "model": {"method": "HF3c"},
        "keywords": {"scf_type": "df"},
    }

    ret = qcng.compute(input_data, "psi4", raise_error=True)

    assert ret.success is True
    assert ret.model.basis is None


@using("psi4_runqcsk")
def test_psi4_interactive_task():
    input_data = {
        "molecule": qcng.get_molecule("water"),
        "driver": "energy",
        "model": {"method": "SCF", "basis": "sto-3g"},
        "keywords": {"scf_type": "df"},
        "extras": {"psiapi": True},
    }

    ret = qcng.compute(input_data, "psi4", raise_error=True)

    assert "Final Energy" in ret.stdout
    assert ret.success
    assert ret.extras.pop("psiapi_evaluated", False)


@using("psi4_runqcsk")
def test_psi4_wavefunction_task():
    input_data = {
        "molecule": qcng.get_molecule("water"),
        "driver": "energy",
        "model": {"method": "SCF", "basis": "sto-3g"},
        "keywords": {"scf_type": "df"},
        "protocols": {"wavefunction": "orbitals_and_eigenvalues"},
    }

    ret = qcng.compute(input_data, "psi4", raise_error=True)
    assert ret.success, ret.error.error_message
    assert ret.wavefunction.scf_orbitals_a.shape == (7, 7)


@using("psi4")
def test_psi4_internal_failure():

    mol = Molecule.from_data(
        """0 3
     O    0.000000000000     0.000000000000    -0.068516245955
    """
    )

    psi4_task = {
        "molecule": mol,
        "driver": "energy",
        "model": {"method": "ccsd", "basis": "6-31g"},
        "keywords": {"reference": "rhf"},
    }
    with pytest.raises(qcng.exceptions.InputError) as exc:
        ret = qcng.compute(psi4_task, "psi4", raise_error=True)

    assert "reference is only" in str(exc.value)


@using("psi4")
def test_psi4_ref_switch():
    inp = AtomicInput(
        **{
            "molecule": {"symbols": ["Li"], "geometry": [0, 0, 0], "molecular_multiplicity": 2},
            "driver": "energy",
            "model": {"method": "B3LYP", "basis": "sto-3g"},
            "keywords": {"scf_type": "df"},
        }
    )

    ret = qcng.compute(inp, "psi4", raise_error=True, return_dict=False)

    assert ret.success is True
    assert ret.properties.calcinfo_nalpha == 2
    assert ret.properties.calcinfo_nbeta == 1


@using("rdkit")
@pytest.mark.parametrize("method", ["UFF", "MMFF94", "MMFF94s"])
def test_rdkit_task(method):
    input_data = {
        "molecule": qcng.get_molecule("water"),
        "driver": "gradient",
        "model": {"method": method},
        "keywords": {},
    }

    ret = qcng.compute(input_data, "rdkit", raise_error=True)

    assert ret.success is True


@using("rdkit")
def test_rdkit_connectivity_error():
    input_data = {
        "molecule": qcng.get_molecule("water").dict(),
        "driver": "energy",
        "model": {"method": "UFF", "basis": ""},
        "keywords": {},
    }
    del input_data["molecule"]["connectivity"]

    ret = qcng.compute(input_data, "rdkit")
    assert ret.success is False
    assert "connectivity" in ret.error.error_message

    with pytest.raises(qcng.exceptions.InputError):
        qcng.compute(input_data, "rdkit", raise_error=True)


@using("torchani")
def test_torchani_task():
    input_data = {
        "molecule": qcng.get_molecule("water"),
        "driver": "gradient",
        "model": {"method": "ANI1x", "basis": None},
        "keywords": {},
    }

    ret = qcng.compute(input_data, "torchani", raise_error=True)

    assert ret.success is True
    assert ret.driver == "gradient"


@using("mopac")
def test_mopac_task():
    input_data = {
        "molecule": qcng.get_molecule("water"),
        "driver": "gradient",
        "model": {"method": "PM6", "basis": None},
        "keywords": {"pulay": False},
    }

    ret = qcng.compute(input_data, "mopac", raise_error=True)
    assert ret.extras.keys() >= {"heat_of_formation", "energy_electronic", "dip_vec"}
    energy = pytest.approx(-0.08474117913025125, rel=1.0e-5)

    # Check gradient
    ret = qcng.compute(input_data, "mopac", raise_error=True)
    assert ret.extras.keys() >= {"heat_of_formation", "energy_electronic", "dip_vec"}
    assert np.linalg.norm(ret.return_result) == pytest.approx(0.03543560156912385, rel=1.0e-4)
    assert ret.properties.return_energy == energy

    # Check energy
    input_data["driver"] = "energy"
    ret = qcng.compute(input_data, "mopac", raise_error=True)
    assert ret.return_result == energy
    assert "== MOPAC DONE ==" in ret.stdout


def test_random_failure_no_retries(failure_engine):

    failure_engine.iter_modes = ["input_error"]
    ret = qcng.compute(failure_engine.get_job(), failure_engine.name, raise_error=False)
    assert ret.error.error_type == "input_error"
    assert "retries" not in ret.input_data["provenance"].keys()

    failure_engine.iter_modes = ["random_error"]
    ret = qcng.compute(failure_engine.get_job(), failure_engine.name, raise_error=False)
    assert ret.error.error_type == "random_error"
    assert "retries" not in ret.input_data["provenance"].keys()


def test_random_failure_with_retries(failure_engine):

    failure_engine.iter_modes = ["random_error", "random_error", "random_error"]
    ret = qcng.compute(failure_engine.get_job(), failure_engine.name, raise_error=False, local_options={"retries": 2})
    assert ret.input_data["provenance"]["retries"] == 2
    assert ret.error.error_type == "random_error"

    failure_engine.iter_modes = ["random_error", "input_error"]
    ret = qcng.compute(failure_engine.get_job(), failure_engine.name, raise_error=False, local_options={"retries": 4})
    assert ret.input_data["provenance"]["retries"] == 1
    assert ret.error.error_type == "input_error"


def test_random_failure_with_success(failure_engine):

    failure_engine.iter_modes = ["random_error", "pass"]
    failure_engine.ncalls = 0
    ret = qcng.compute(failure_engine.get_job(), failure_engine.name, raise_error=False, local_options={"retries": 1})

    assert ret.success, ret.error.error_message
    assert ret.provenance.retries == 1
    assert ret.extras["ncalls"] == 2


@using("openmm")
def test_openmm_task_smirnoff():
    from qcengine.programs.openmm import OpenMMHarness

    input_data = {
        "molecule": qcng.get_molecule("water"),
        "driver": "energy",
        "model": {"method": "openff-1.0.0", "basis": "smirnoff"},
        "keywords": {},
    }

    ret = qcng.compute(input_data, "openmm", raise_error=True)

    cachelength = len(OpenMMHarness._CACHE)

    assert cachelength > 0
    assert ret.success is True

    ret = qcng.compute(input_data, "openmm", raise_error=True)

    # ensure cache has not grown
    assert len(OpenMMHarness._CACHE) == cachelength
    assert ret.success is True


@pytest.mark.skip("`basis` must be explicitly specified at this time")
@using("openmm")
def test_openmm_task_url_basis():
    from qcengine.programs.openmm import OpenMMHarness

    input_data = {
        "molecule": qcng.get_molecule("water"),
        "driver": "energy",
        "model": {
            "method": "openmm",
            "basis": "openff-1.0.0",
            "url": "https://raw.githubusercontent.com/openforcefield/openff-forcefields/1.0.0/openforcefields/offxml/openff-1.0.0.offxml",
        },
        "keywords": {},
    }

    ret = qcng.compute(input_data, "openmm", raise_error=True)

    cachelength = len(OpenMMHarness._CACHE)

    assert cachelength > 0
    assert ret.success is True

    ret = qcng.compute(input_data, "openmm", raise_error=True)

    # ensure cache has not grown
    assert len(OpenMMHarness._CACHE) == cachelength
    assert ret.success is True


@using("openmm")
def test_openmm_cmiles_gradient():
    program = "openmm"

    water = qcng.get_molecule("water")

    water_dict = water.dict()
    # add water cmiles to the molecule
    water_dict["extras"] = {"cmiles": {"canonical_isomeric_explicit_hydrogen_mapped_smiles": "[H:2][O:1][H:3]"}}

    molecule = Molecule.from_data(water_dict)

    model = {"method": "openff-1.0.0", "basis": "smirnoff"}

    inp = AtomicInput(molecule=molecule, driver="gradient", model=model)
    ret = qcng.compute(inp, program, raise_error=False)

    assert ret.success is True


@using("openmm")
def test_openmm_cmiles_gradient_nomatch():
    program = "openmm"

    water = qcng.get_molecule("water")

    water_dict = water.dict()
    # add ethane cmiles to the molecule
    water_dict["extras"] = {
        "cmiles": {
            "canonical_isomeric_explicit_hydrogen_mapped_smiles": "[H:3][C:1]([H:4])([H:5])[C:2]([H:6])([H:7])[H:8]"
        }
    }

    molecule = Molecule.from_data(water_dict)

    model = {"method": "openff-1.0.0", "basis": "smirnoff"}

    inp = AtomicInput(molecule=molecule, driver="gradient", model=model)
    ret = qcng.compute(inp, program, raise_error=False)

    # if we correctly find the cmiles this should fail as the molecule and cmiles are different
    assert ret.success is False
    assert (
        "molecule.add_conformer given input of the wrong shape: Given (3, 3), expected (8, 3)"
        in ret.error.error_message
    )


@using("openmm")
@pytest.mark.parametrize(
    "gaff_settings",
    [
        pytest.param(({}, None, 0.0013904199062156914), id="gaff no keywords"),
        pytest.param(({"constraints": "ALLBONDS"}, None, 8.108238580315493e-05), id="constraints allbonds"),
        pytest.param(({"nonbondedMethod": "LjPmE"}, None, 0.0013904199062156914), id="nonbonded ljpme"),
        pytest.param(
            ({"nonbondedMethod": "PME", "constraints": "Hbonds"}, None, 8.108238580315493e-05),
            id="nonbonded pme constraints hbonds",
        ),
        pytest.param(({"constraints": "badmethod"}, ValueError, 0), id="incorrect constraint"),
        pytest.param(({"nonbondedMethod": "badmethod"}, ValueError, 0), id="incorrect nonbondedmethod"),
    ],
)
def test_openmm_gaff_keywords(gaff_settings):
    """
    Test the different running settings with gaff.
    """
    program = "openmm"
    water = qcng.get_molecule("water")

    water_dict = water.dict()
    # add water cmiles to the molecule
    water_dict["extras"] = {"cmiles": {"canonical_isomeric_explicit_hydrogen_mapped_smiles": "[H:2][O:1][H:3]"}}

    molecule = Molecule.from_data(water_dict)
    keywords, error, expected_result = gaff_settings
    model = {"method": "gaff-2.1", "basis": "antechamber"}
    inp = AtomicInput(molecule=molecule, driver="energy", model=model, keywords=keywords)
    if error is not None:
        with pytest.raises(error):
            _ = qcng.compute(inp, program, raise_error=True)
    else:
        ret = qcng.compute(inp, program, raise_error=False)
        assert ret.success is True
        assert ret.return_result == pytest.approx(expected_result, rel=1e-6)
