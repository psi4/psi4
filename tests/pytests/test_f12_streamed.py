import pytest
import psi4

from addons import uusing


pytestmark = [pytest.mark.api, pytest.mark.quick]

WATER = "O 0 0 0\nH 0 0.757 0.586\nH 0 -0.757 0.586"
COMPONENTS = (
    "HF-CABS TOTAL ENERGY",
    "MP2 CORRELATION ENERGY",
    "MP2-F12 CORRECTION ENERGY",
    "MP2-F12 CORRELATION ENERGY",
    "MP2-F12 OPPOSITE-SPIN CORRELATION ENERGY",
    "MP2-F12 SAME-SPIN CORRELATION ENERGY",
)


@uusing("einsums")
@pytest.mark.parametrize(
    "geometry,basis,frozen,beta,singles,puream,threads,block",
    [
        pytest.param(WATER, "cc-pvdz-f12", True, 1.0, True, True, 1, 64, id="water"),
        pytest.param(WATER, "cc-pvdz-f12", False, 1.4, False, True, 4, 7, id="unfrozen-no-singles"),
        pytest.param("H 0 0 0\nF 0 0 0.92", "aug-cc-pvdz", True, 1.2, True, True, 1, 7, id="augmented-basis"),
        pytest.param(WATER, "cc-pvdz-f12", True, 1.0, True, False, 4, 64, id="cartesian-truncated-ri"),
        pytest.param("H 0 0 0\nH 0 0 0.74", "cc-pvdz-f12", False, 1.0, True, True, 2, 1, id="one-occupied-pair"),
        pytest.param("H 0 0 0\nH 0 0 0.74", "cc-pvtz-f12", False, 1.0, True, True, 1, 7, id="triple-zeta"),
    ],
)
def test_streamed_matches_incore(geometry, basis, frozen, beta, singles, puream, threads, block, s_tolerance=1e-7):
    molecule = psi4.geometry("0 1\n" + geometry + "\nunits angstrom\nsymmetry c1\nno_com\nno_reorient")
    psi4.set_memory("2 GB")
    psi4.set_num_threads(threads)
    psi4.set_options({
        "basis": basis,
        "cabs_basis": basis + "-optri",
        "df_basis_f12": "aug-cc-pvdz-ri",
        "df_basis_scf": "aug-cc-pvdz-ri",
        "reference": "rhf",
        "scf_type": "df",
        "mp2_type": "df",
        "freeze_core": frozen,
        "f12_beta": beta,
        "cabs_singles": singles,
        "puream": puream,
        "s_tolerance": s_tolerance,
        "e_convergence": 1e-10,
        "d_convergence": 1e-10,
        "f12_subtype": "incore",
    })
    _, reference = psi4.energy("scf", molecule=molecule, return_wfn=True)
    expected_energy, expected = psi4.energy("mp2-f12", ref_wfn=reference, return_wfn=True)
    if s_tolerance > 1e-7:
        assert expected.nmo() < expected.basisset().nbf()
    if not puream:
        combined = expected.get_basisset("CABS")
        ri = psi4.core.OrbitalSpace.build_ri_space(combined, 1e-8)
        assert ri.dim()[0] < combined.nbf()
    components = {name: expected.variable(name) for name in COMPONENTS}
    psi4.core.clean()
    psi4.set_options({"f12_subtype": "streamed", "f12_aux_block_size": block})
    screening = psi4.core.get_global_option("SCREENING")
    energy, result = psi4.energy("mp2-f12", ref_wfn=reference, return_wfn=True)
    assert energy == pytest.approx(expected_energy, abs=1e-8, rel=0)
    assert result.energy() == pytest.approx(energy, abs=1e-12, rel=0)
    for name, value in components.items():
        assert result.variable(name) == pytest.approx(value, abs=1e-8, rel=0), name
    assert psi4.core.get_global_option("SCREENING") == screening
    if not singles:
        assert result.variable("F12 CABS CORRECTION ENERGY") == 0.0


@uusing("einsums")
def test_streamed_truncated_orbital_space():
    # An elevated overlap cutoff forces a rectangular AO-to-MO transform.
    test_streamed_matches_incore(WATER, "cc-pvdz-f12", True, 1.0, True, True, 4, 7, s_tolerance=1e-2)


@uusing("einsums")
@pytest.mark.parametrize("option,value,message", [
    ("MP2_TYPE", "CONV", "requires MP2_TYPE=DF"),
    ("F12_READ_INTS", True, "does not read saved F12 integrals"),
    ("F12_AUX_BLOCK_SIZE", 0, "must be positive"),
])
def test_streamed_rejects_incompatible_options(option, value, message):
    molecule = psi4.geometry("0 1\nH 0 0 0\nH 0 0 0.74\nsymmetry c1")
    psi4.set_options({"basis": "sto-3g", "reference": "rhf", "scf_type": "pk"})
    _, reference = psi4.energy("scf", molecule=molecule, return_wfn=True)
    psi4.core.set_local_option("F12", "F12_SUBTYPE", "STREAMED")
    psi4.core.set_local_option("F12", option, value)
    psi4.core.prepare_options_for_module("F12")
    with pytest.raises(RuntimeError, match=message):
        psi4.core.f12(reference)


@uusing("einsums")
@uusing("ecpint")
def test_streamed_rejects_ecp_reference():
    molecule = psi4.geometry("-1 1\nI 0 0 0\nsymmetry c1")
    psi4.set_num_threads(1)
    psi4.set_options({"basis": "def2-svp", "reference": "rhf", "scf_type": "pk"})
    _, reference = psi4.energy("scf", molecule=molecule, return_wfn=True)
    psi4.core.set_local_option("F12", "F12_SUBTYPE", "STREAMED")
    psi4.core.prepare_options_for_module("F12")
    with pytest.raises(RuntimeError, match="all-electron basis sets only"):
        psi4.core.f12(reference)
