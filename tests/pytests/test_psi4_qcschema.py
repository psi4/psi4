import sys
import json
import numpy as np
import pytest
import pprint
import os
from shutil import copytree

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick, pytest.mark.smoke]

# Generating
# * see test_psi4.py for generation
# * equivalent to test_psi4. copy over the job, then run below to generate atomicinput
#    atin = psi4.driver.p4util.state_to_atomicinput(driver="energy", method="ccsd", molecule=ethene_ethyne)
#    print(f'    jatin = """{atin.serialize("json")}"""')
#    assert 0
# * switch to json running
#    atres = psi4.schema_wrapper.run_qcschema(jatin)
#    pprint.pprint(atres.dict())


_sch_name = {1: "qcschema_output", 2: "qcschema_atomic_result"}
_ispy314 = sys.version_info >= (3, 14)


@pytest.fixture
def datadir(tmpdir, request):
    """
    from: https://stackoverflow.com/a/29631801
    Fixture responsible for searching a folder with the same name of test
    module and, if available, moving all contents to a temporary directory so
    tests can use them freely.
    """
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)

    if os.path.isdir(test_dir):
        copytree(test_dir, str(tmpdir), dirs_exist_ok=True)
    else:
        raise FileNotFoundError("Test folder not found.")

    return tmpdir


@pytest.mark.parametrize("schver", [pytest.param(1, marks=[pytest.mark.skipif(_ispy314, reason="Py314+QCSk1")]), pytest.param(2)])
def test_psi4_basic(datadir, schver):
    """tu1-h2o-energy"""
    #! Sample HF/cc-pVDZ H2O computation

    ref_file = datadir.join(f"jatin1-{schver}.ref")
    with open(ref_file) as f:
        jatin = json.load(f)

    atres = psi4.schema_wrapper.run_qcschema(jatin)
    assert atres.success, pprint.pprint(atres.dict(), width=200)
    assert atres.schema_version == schver
    assert atres.schema_name == _sch_name[schver]
    assert psi4.compare_values(-76.0266327341067125, atres.return_result, 6, 'SCF energy')


@pytest.mark.parametrize("schver", [pytest.param(1, marks=[pytest.mark.skipif(_ispy314, reason="Py314+QCSk1")]), pytest.param(2)])
def test_psi4_cc(datadir, schver):
    """cc1"""
    #! RHF-CCSD 6-31G** all-electron optimization of the H2O molecule

    ref_file = datadir.join(f"jatin2-{schver}.ref")
    with open(ref_file) as f:
        jatin = json.load(f)

    atres = psi4.schema_wrapper.run_qcschema(jatin)

    if "qcengine.exceptions.InputError" in atres.error.error_message:
        pytest.xfail("no AtomicInput optimization")

    # psi4.optimize('ccsd')

    refnuc = 9.1654609427539
    refscf = -76.0229427274435
    refccsd = -0.20823570806196
    reftotal = -76.2311784355056

    assert psi4.compare_values(refnuc, atres.properties.nuclear_repulsion_energy, 3, "Nuclear repulsion energy")
    assert psi4.compare_values(refscf, atres.extras["qcvars"]["SCF total energy"], 5, "SCF energy")
    assert psi4.compare_values(refccsd, psi4.variable("CCSD correlation energy"), 4, "CCSD contribution")
    assert psi4.compare_values(reftotal, psi4.variable("Current energy"), 7, "Total energy")


@pytest.mark.parametrize("schver", [pytest.param(1, marks=[pytest.mark.skipif(_ispy314, reason="Py314+QCSk1")]), pytest.param(2)])
def test_psi4_cas(datadir, schver):
    """casscf-sp"""
    #! CASSCF/6-31G** energy point

    ref_file = datadir.join(f"jatin3-{schver}.ref")
    with open(ref_file) as f:
        jatin = json.load(f)

    atres = psi4.schema_wrapper.run_qcschema(jatin)
    assert atres.success, pprint.pprint(atres.dict(), width=200)
    assert atres.schema_version == schver
    assert atres.schema_name == _sch_name[schver]
    assert psi4.compare_values(-76.2198474477531, atres.return_result, 6, 'CISD Energy')

    ref_file = datadir.join(f"jatin4-{schver}.ref")
    with open(ref_file) as f:
        jatin = json.load(f)

    atres = psi4.schema_wrapper.run_qcschema(jatin)
    assert atres.success, pprint.pprint(atres.dict(), width=200)
    assert atres.schema_version == schver
    assert atres.schema_name == _sch_name[schver]
    assert psi4.compare_values(-76.073865006902, atres.return_result, 6, 'CASSCF Energy')


@pytest.mark.nbody
@pytest.mark.parametrize("schver", [pytest.param(1, marks=[pytest.mark.skipif(_ispy314, reason="Py314+QCSk1")]), pytest.param(2)])
def test_psi4_dfmp2(datadir, schver):
    """dfmp2-1"""
    #! Density fitted MP2 cc-PVDZ/cc-pVDZ-RI computation of formic acid dimer binding energy
    #! using automatic counterpoise correction.  Monomers are specified using Cartesian coordinates.

    Enuc = 235.94662124
    Ecp = -0.0224119246

    ref_file = datadir.join(f"jatin5-{schver}.ref")
    with open(ref_file) as f:
        jatin = json.load(f)

    atres = psi4.schema_wrapper.run_qcschema(jatin)
    assert atres.success, pprint.pprint(atres.dict(), width=200)
    assert atres.schema_version == schver
    assert atres.schema_name == _sch_name[schver]
    assert psi4.compare_values(Enuc, atres.properties.nuclear_repulsion_energy, 7, "Nuclear Repulsion Energy")
    assert psi4.compare_values(Ecp, atres.return_result, 5, "CP Corrected cc-pVDZ/cc-pVDZ-RI DFMP2")


@pytest.mark.parametrize("schver", [pytest.param(1, marks=[pytest.mark.skipif(_ispy314, reason="Py314+QCSk1")]), pytest.param(2)])
def test_psi4_sapt(datadir, schver):
    """sapt1"""
    #! SAPT0 cc-pVDZ computation of the ethene-ethyne interaction energy, using the cc-pVDZ-JKFIT RI basis for SCF
    #! and cc-pVDZ-RI for SAPT.  Monomer geometries are specified using Cartesian coordinates.

    # TODO: add `"extras": {"wfn_qcvars_only": true}` when SAPT stored on wfn properly

    Eref = [85.1890645313, -0.00359915058, 0.00362911158, -0.00083137117, -0.00150542374, -0.00230683391]

    ref_file = datadir.join(f"jatin6-{schver}.ref")
    with open(ref_file) as f:
        jatin = json.load(f)

    atres = psi4.schema_wrapper.run_qcschema(jatin)
    assert atres.success, pprint.pprint(atres.dict(), width=200)
    assert atres.schema_version == schver
    assert atres.schema_name == _sch_name[schver]

    assert psi4.compare_values(Eref[0], atres.properties.nuclear_repulsion_energy, 9, "Nuclear Repulsion Energy")
    assert psi4.compare_values(Eref[1], atres.extras["qcvars"]["SAPT ELST ENERGY"], 6, "SAPT0 Eelst")
    assert psi4.compare_values(Eref[2], atres.extras["qcvars"]["SAPT EXCH ENERGY"], 6, "SAPT0 Eexch")
    assert psi4.compare_values(Eref[3], atres.extras["qcvars"]["SAPT IND ENERGY"], 6, "SAPT0 Eind")
    assert psi4.compare_values(Eref[4], atres.extras["qcvars"]["SAPT DISP ENERGY"], 6, "SAPT0 Edisp")
    assert psi4.compare_values(Eref[5], atres.extras["qcvars"]["SAPT0 TOTAL ENERGY"], 6, "SAPT0 Etotal")


@pytest.mark.parametrize("schver", [pytest.param(1, marks=[pytest.mark.skipif(_ispy314, reason="Py314+QCSk1")]), pytest.param(2)])
def test_psi4_scfproperty(datadir, schver):
    """scf-property"""
    #! UFH and B3LYP cc-pVQZ properties for the CH2 molecule.

    ref_hf_di_au = np.array([0.0, 0.0, 0.22531665104559076])
    ref_hf_quad_au = np.array([[-5.69804565317, 0.0, 0.0], [0.0, -4.53353128969, 0.0], [0.0, 0.0, -5.25978856037]])
    ref_b3lyp_di_au = np.array([0.0, 0.0, 0.252480541747])
    ref_b3lyp_quad_au = np.array([[-5.66266837697, 0.0, 0.0], [0.0, -4.46523692003, 0.0], [0.0, 0.0, -5.22054902407]])

    props = [
        'DIPOLE', 'QUADRUPOLE', 'MULLIKEN_CHARGES', 'LOWDIN_CHARGES', 'LOWDIN_SPINS', 'WIBERG_LOWDIN_INDICES', 'MAYER_INDICES',
        'MAYER_INDICES', 'MO_EXTENTS', 'GRID_FIELD', 'GRID_ESP', 'ESP_AT_NUCLEI', 'MULTIPOLE(5)', 'NO_OCCUPATIONS'
    ]

    ref_file = datadir.join(f"jatin7-{schver}.ref")
    with open(ref_file) as f:
        jatin = json.load(f)

    atres = psi4.schema_wrapper.run_qcschema(jatin)
    assert atres.success, pprint.pprint(atres.dict(), width=200)
    assert atres.schema_version == schver
    assert atres.schema_name == _sch_name[schver]
    assert psi4.compare_values(-38.91591819679808, atres.extras["qcvars"]["CURRENT ENERGY"], 6, "SCF energy")
    assert psi4.compare_values(ref_hf_di_au, atres.extras["qcvars"]["SCF DIPOLE"], 4, "SCF DIPOLE")
    assert psi4.compare_values(ref_hf_quad_au, atres.extras["qcvars"]["SCF QUADRUPOLE"], 4, "SCF QUADRUPOLE")
    if schver == 1:
        assert "grid.dat" in atres.extras["extra_infiles"]
    elif schver == 2:
        assert "grid.dat" in atres.input_data.specification.extras["extra_infiles"]
        assert "extra_infiles" not in atres.extras
    print(atres.native_files["grid_field.dat"])

    ref_file = datadir.join(f"jatin8-{schver}.ref")
    with open(ref_file) as f:
        jatin = json.load(f)

    atres = psi4.schema_wrapper.run_qcschema(jatin)
    assert atres.success, pprint.pprint(atres.dict(), width=200)
    assert atres.schema_version == schver
    assert atres.schema_name == _sch_name[schver]
    assert psi4.compare_values(-39.14134740550916, atres.extras["qcvars"]["CURRENT ENERGY"], 6, "B3LYP energy")
    assert psi4.compare_values(ref_b3lyp_di_au, atres.extras["qcvars"]["B3LYP DIPOLE"], 4, "B3LYP DIPOLE")
    assert psi4.compare_values(ref_b3lyp_quad_au, atres.extras["qcvars"]["B3LYP QUADRUPOLE"], 4, "B3LYP QUADRUPOLE")
    print(atres.native_files["grid_esp.dat"])
