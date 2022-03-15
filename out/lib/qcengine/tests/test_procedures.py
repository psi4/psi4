"""
Tests the DQM compute dispatch module
"""

import pytest
from qcelemental.models import DriverEnum, OptimizationInput
from qcelemental.models.common_models import Model
from qcelemental.models.procedures import OptimizationSpecification, QCInputSpecification, TDKeywords, TorsionDriveInput

import qcengine as qcng
from qcengine.testing import failure_engine, using


@pytest.fixture(scope="function")
def input_data():
    return {
        "keywords": {"coordsys": "tric", "maxiter": 100, "program": None},
        "input_specification": {"driver": "gradient", "model": None, "keywords": {}},
        "initial_molecule": None,
    }


@using("psi4")
@pytest.mark.parametrize(
    "optimizer",
    [
        pytest.param("geometric", marks=using("geometric")),
        pytest.param("optking", marks=using("optking")),
        pytest.param("berny", marks=using("berny")),
    ],
)
def test_geometric_psi4(input_data, optimizer):

    input_data["initial_molecule"] = qcng.get_molecule("hydrogen")
    input_data["input_specification"]["model"] = {"method": "HF", "basis": "sto-3g"}
    input_data["input_specification"]["keywords"] = {"scf_properties": ["wiberg_lowdin_indices"]}
    input_data["keywords"]["program"] = "psi4"

    input_data = OptimizationInput(**input_data)

    ret = qcng.compute_procedure(input_data, optimizer, raise_error=True)
    assert 10 > len(ret.trajectory) > 1

    assert pytest.approx(ret.final_molecule.measure([0, 1]), 1.0e-4) == 1.3459150737
    assert ret.provenance.creator.lower() == optimizer
    assert ret.trajectory[0].provenance.creator.lower() == "psi4"

    # Check keywords passing
    for single in ret.trajectory:
        assert "scf_properties" in single.keywords
        assert "WIBERG_LOWDIN_INDICES" in single.extras["qcvars"]


@using("psi4")
@using("geometric")
def test_geometric_local_options(input_data):

    input_data["initial_molecule"] = qcng.get_molecule("hydrogen")
    input_data["input_specification"]["model"] = {"method": "HF", "basis": "sto-3g"}
    input_data["keywords"]["program"] = "psi4"

    input_data = OptimizationInput(**input_data)

    # Set some extremely large number to test
    ret = qcng.compute_procedure(input_data, "geometric", raise_error=True, local_options={"memory": "5000"})
    assert pytest.approx(ret.trajectory[0].provenance.memory, 1) == 4900

    # Make sure we cleaned up
    assert "_qcengine_local_config" not in ret.input_specification
    assert "_qcengine_local_config" not in ret.trajectory[0].extras


@using("rdkit")
@using("geometric")
def test_geometric_stdout(input_data):

    input_data["initial_molecule"] = qcng.get_molecule("water")
    input_data["input_specification"]["model"] = {"method": "UFF", "basis": ""}
    input_data["keywords"]["program"] = "rdkit"

    input_data = OptimizationInput(**input_data)

    ret = qcng.compute_procedure(input_data, "geometric", raise_error=True)
    assert ret.success is True
    assert "Converged!" in ret.stdout


@using("psi4")
@using("berny")
def test_berny_stdout(input_data):

    input_data["initial_molecule"] = qcng.get_molecule("water")
    input_data["input_specification"]["model"] = {"method": "HF", "basis": "sto-3g"}
    input_data["keywords"]["program"] = "psi4"

    input_data = OptimizationInput(**input_data)

    ret = qcng.compute_procedure(input_data, "berny", raise_error=True)
    assert ret.success is True
    assert "All criteria matched" in ret.stdout


@using("geometric")
@using("rdkit")
def test_geometric_rdkit_error(input_data):

    input_data["initial_molecule"] = qcng.get_molecule("water").copy(exclude={"connectivity_"})
    input_data["input_specification"]["model"] = {"method": "UFF", "basis": ""}
    input_data["keywords"]["program"] = "rdkit"

    input_data = OptimizationInput(**input_data)

    ret = qcng.compute_procedure(input_data, "geometric")
    assert ret.success is False
    assert isinstance(ret.error.error_message, str)


@using("rdkit")
@using("geometric")
def test_optimization_protocols(input_data):

    input_data["initial_molecule"] = qcng.get_molecule("water")
    input_data["input_specification"]["model"] = {"method": "UFF"}
    input_data["keywords"]["program"] = "rdkit"
    input_data["protocols"] = {"trajectory": "initial_and_final"}

    input_data = OptimizationInput(**input_data)

    ret = qcng.compute_procedure(input_data, "geometric", raise_error=True)
    assert ret.success, ret.error.error_message

    assert len(ret.trajectory) == 2
    assert ret.initial_molecule.get_hash() == ret.trajectory[0].molecule.get_hash()
    assert ret.final_molecule.get_hash() == ret.trajectory[1].molecule.get_hash()


@using("geometric")
@using("rdkit")
def test_geometric_retries(failure_engine, input_data):

    failure_engine.iter_modes = ["random_error", "pass", "random_error", "random_error", "pass"]  # Iter 1  # Iter 2
    failure_engine.iter_modes.extend(["pass"] * 20)

    input_data["initial_molecule"] = {
        "symbols": ["He", "He"],
        "geometry": [0, 0, 0, 0, 0, failure_engine.start_distance],
    }
    input_data["input_specification"]["model"] = {"method": "something"}
    input_data["keywords"]["program"] = failure_engine.name

    input_data = OptimizationInput(**input_data)

    ret = qcng.compute_procedure(input_data, "geometric", local_options={"ncores": 13}, raise_error=True)
    assert ret.success is True
    assert ret.trajectory[0].provenance.retries == 1
    assert ret.trajectory[0].provenance.ncores == 13
    assert ret.trajectory[1].provenance.retries == 2
    assert ret.trajectory[1].provenance.ncores == 13
    assert "retries" not in ret.trajectory[2].provenance.dict()

    # Ensure we still fail
    failure_engine.iter_modes = ["random_error", "pass", "random_error", "random_error", "pass"]  # Iter 1  # Iter 2
    ret = qcng.compute_procedure(input_data, "geometric", local_options={"ncores": 13, "retries": 1})
    assert ret.success is False
    assert ret.input_data["trajectory"][0]["provenance"]["retries"] == 1
    assert len(ret.input_data["trajectory"]) == 2


@using("geometric")
@pytest.mark.parametrize(
    "program, model, bench",
    [
        pytest.param(
            "rdkit", {"method": "UFF"}, [1.87130923886072, 2.959448636243545, 104.5099642579023], marks=using("rdkit")
        ),
        pytest.param(
            "rdkit",
            {"method": "mmff94"},
            [1.8310842343589573, 2.884612338953529, 103.93822919865106],
            marks=using("rdkit"),
        ),
        pytest.param(
            "rdkit",
            {"method": "MMFF94s"},
            [1.8310842343589573, 2.884612338953529, 103.93822919865106],
            marks=using("rdkit"),
        ),
        pytest.param(
            "torchani",
            {"method": "ANI1x"},
            [1.82581873750194, 2.866376526793269, 103.4332610730292],
            marks=using("torchani"),
        ),
        pytest.param(
            "mopac",
            {"method": "PM6"},
            [1.7927843431811934, 2.893333237502448, 107.60441967992045],
            marks=using("mopac"),
        ),
        pytest.param(
            "openmm",
            {"method": "openff-1.0.0", "basis": "smirnoff"},
            [1.8344994291796748, 3.010099477501204, 110.25177977849998],
            marks=using("openmm"),
        ),
        pytest.param(
            "openmm",
            {"method": "openff_unconstrained-1.0.0", "basis": "smirnoff"},
            [1.8344994291195869, 3.0100994772976124, 110.25259556886984],
            marks=using("openmm"),
        ),
        pytest.param(
            "openmm",
            {"method": "smirnoff99Frosst-1.1.0", "basis": "smirnoff"},
            [1.814137087600702, 3.025566213038376, 112.9999999990053],
            marks=using("openmm"),
        ),
        pytest.param(
            "qcore",
            {"method": "GFN1"},
            [1.8104763949897031, 2.9132449420655213, 107.13403040879244],
            marks=using("qcore"),
        ),
    ],
)
def test_geometric_generic(input_data, program, model, bench):

    input_data["initial_molecule"] = qcng.get_molecule("water")
    input_data["input_specification"]["model"] = model
    input_data["keywords"]["program"] = program
    input_data["input_specification"]["extras"] = {"_secret_tags": {"mysecret_tag": "data1"}}

    ret = qcng.compute_procedure(input_data, "geometric", raise_error=True)
    assert ret.success is True
    assert "Converged!" in ret.stdout

    r01, r02, r12, a102 = ret.final_molecule.measure([[0, 1], [0, 2], [1, 2], [1, 0, 2]])

    assert pytest.approx(r01, 1.0e-4) == bench[0]
    assert pytest.approx(r02, 1.0e-4) == bench[0]
    assert pytest.approx(r12, 1.0e-4) == bench[1]
    assert pytest.approx(a102, 1.0e-4) == bench[2]

    assert "_secret_tags" in ret.trajectory[0].extras
    assert "data1" == ret.trajectory[0].extras["_secret_tags"]["mysecret_tag"]


@using("nwchem")
@pytest.mark.parametrize("linopt", [0, 1])
def test_nwchem_relax(linopt):
    # Make the input file
    input_data = {
        "input_specification": {
            "model": {"method": "HF", "basis": "sto-3g"},
            "keywords": {"set__driver:linopt": linopt},
        },
        "initial_molecule": qcng.get_molecule("hydrogen"),
    }
    input_data = OptimizationInput(**input_data)

    # Run the relaxation
    ret = qcng.compute_procedure(input_data, "nwchemdriver", raise_error=True)
    assert 10 > len(ret.trajectory) > 1

    assert pytest.approx(ret.final_molecule.measure([0, 1]), 1.0e-4) == 1.3459150737


@using("rdkit")
@using("torsiondrive")
def test_torsiondrive_generic():

    input_data = TorsionDriveInput(
        keywords=TDKeywords(dihedrals=[(2, 0, 1, 5)], grid_spacing=[180]),
        input_specification=QCInputSpecification(driver=DriverEnum.gradient, model=Model(method="UFF", basis=None)),
        initial_molecule=qcng.get_molecule("ethane"),
        optimization_spec=OptimizationSpecification(
            procedure="geomeTRIC",
            keywords={
                "coordsys": "dlc",
                "maxiter": 300,
                "program": "rdkit",
            },
        ),
    )

    ret = qcng.compute_procedure(input_data, "torsiondrive", raise_error=True)

    assert ret.error is None
    assert ret.success

    expected_grid_ids = {"180", "0"}

    assert {*ret.optimization_history} == expected_grid_ids

    assert {*ret.final_energies} == expected_grid_ids
    assert {*ret.final_molecules} == expected_grid_ids

    assert (
        pytest.approx(ret.final_molecules["180"].measure([2, 0, 1, 5]), abs=1.0e-2) == 180.0
        or pytest.approx(ret.final_molecules["180"].measure([2, 0, 1, 5]), abs=1.0e-2) == -180.0
    )
    assert pytest.approx(ret.final_molecules["0"].measure([2, 0, 1, 5]), abs=1.0e-2) == 0.0

    assert ret.provenance.creator.lower() == "torsiondrive"
    assert ret.optimization_history["180"][0].provenance.creator.lower() == "geometric"
    assert ret.optimization_history["180"][0].trajectory[0].provenance.creator.lower() == "rdkit"

    assert ret.stdout == "All optimizations converged at lowest energy. Job Finished!\n"
