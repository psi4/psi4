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


@uusing("einsums")
@pytest.mark.parametrize("failure", ["missing-cabs", "memory"])
def test_f12_recovers_after_failure(failure, monkeypatch):
    from psi4.driver.procrouting import proc

    molecule = psi4.geometry("0 1\nH 0 0 0\nH 0 0 0.74\nsymmetry c1")
    psi4.set_memory("2 GB")
    psi4.set_num_threads(1)
    psi4.set_options({
        "basis": "sto-3g" if failure == "missing-cabs" else "cc-pvdz-f12",
        "df_basis_f12": "aug-cc-pvdz-ri",
        "f12_subtype": "streamed",
        "mp2_type": "df",
        "reference": "rhf",
    })
    # Allow the F12 driver to select DF-SCF and the default CABS itself.
    psi4.core.set_global_option("SCF_TYPE", "PK")
    psi4.core.revoke_global_option_changed("SCF_TYPE")

    def options():
        values = [(psi4.core.get_global_option(key), psi4.core.has_global_option_changed(key))
                  for key in ("SCF_TYPE", "MP2_TYPE", "DF_BASIS_F12", "SCREENING",
                              "CABS_BASIS", "DF_BASIS_SCF", "DF_BASIS_MP2")]
        values.extend((psi4.core.get_local_option(module, key),
                       psi4.core.has_local_option_changed(module, key))
                      for module, key in (("F12", "CABS_BASIS"), ("SCF", "DF_BASIS_SCF"),
                                          ("DFMP2", "DF_BASIS_MP2")))
        return values

    before = options()
    memory = psi4.core.get_memory()
    factory = psi4.core.f12
    calls = []
    start, stop = psi4.core.tstart, psi4.core.tstop

    def tstart():
        calls.append("start")
        return start()

    def tstop():
        calls.append("stop")
        return stop()

    def small_budget(reference):
        result = factory(reference)
        # Lower only the F12 budget, after the real SCF and MP2 steps finish.
        psi4.core.set_memory_bytes(1, quiet=True)
        return result

    monkeypatch.setattr(psi4.core, "tstart", tstart)
    monkeypatch.setattr(psi4.core, "tstop", tstop)
    with monkeypatch.context() as patch:
        if failure == "memory":
            patch.setattr(psi4.core, "f12", small_budget)
        error = "No CABS_BASIS given" if failure == "missing-cabs" else "STREAMED arrays exceed"
        try:
            with pytest.raises((RuntimeError, psi4.ValidationError), match=error):
                proc.run_mp2f12("mp2-f12", molecule=molecule)
        finally:
            psi4.core.set_memory_bytes(memory, quiet=True)

    assert options() == before
    assert calls.count("start") == calls.count("stop")
    psi4.core.timer_on("after failed F12")
    psi4.core.timer_off("after failed F12")
    probe = next(record for record in psi4.core.get_timer_records().values()
                 if record["timer_name"] == "after failed F12")
    assert "MP2-F12 Compute Energy" not in probe["timer_path"]
    assert "OBS and CABS" not in probe["timer_path"]

    # No clean(), clean_options(), or clean_timers() between failure and retry.
    psi4.set_options({"basis": "cc-pvdz-f12"})
    recovered = proc.run_mp2f12("mp2-f12", molecule=molecule)
    assert recovered.variable("MP2-F12 TOTAL ENERGY") < -1.0
    assert options() == before
    assert calls.count("start") == calls.count("stop")
    psi4.set_options({"f12_subtype": "incore"})
    expected = proc.run_mp2f12("mp2-f12", molecule=molecule)
    assert recovered.variable("MP2-F12 TOTAL ENERGY") == pytest.approx(
        expected.variable("MP2-F12 TOTAL ENERGY"), abs=1e-8, rel=0)


@uusing("einsums")
def test_streamed_repeated_thread_counts():
    molecule = psi4.geometry("0 1\n" + WATER + "\nsymmetry c1")
    psi4.set_memory("2 GB")
    psi4.set_num_threads(1)
    psi4.set_options({
        "basis": "cc-pvdz-f12",
        "df_basis_scf": "aug-cc-pvdz-ri",
        "df_basis_f12": "aug-cc-pvdz-ri",
        "scf_type": "df",
        "mp2_type": "df",
        "f12_subtype": "streamed",
        "freeze_core": True,
        "e_convergence": 1e-10,
        "d_convergence": 1e-10,
    })
    _, reference = psi4.energy("scf", molecule=molecule, return_wfn=True)
    expected = None
    for threads in (1, 4, 4):
        psi4.set_num_threads(threads)
        energy, result = psi4.energy("mp2-f12", ref_wfn=reference, return_wfn=True)
        values = [energy] + [result.variable(key) for key in COMPONENTS]
        if expected is None:
            expected = values
        else:
            assert values == pytest.approx(expected, abs=1e-9, rel=0)


@uusing("einsums")
def test_streamed_reduces_pair_workers(monkeypatch, tmp_path):
    molecule = psi4.geometry("0 1\n" + WATER + "\nsymmetry c1")
    psi4.set_memory("2 GB")
    psi4.set_num_threads(4)
    # Small orbital/CABS spaces keep this a cheap scheduling test. With these
    # spaces, 1 MiB fits the serial stages and two pair workspaces, but not four.
    psi4.set_options({
        "basis": "sto-3g",
        "cabs_basis": "3-21g",
        "df_basis_scf": "aug-cc-pvdz-ri",
        "df_basis_f12": "aug-cc-pvdz-ri",
        "scf_type": "df",
        "mp2_type": "df",
        "f12_subtype": "streamed",
        "f12_aux_block_size": 1,
        "freeze_core": True,
        "e_convergence": 1e-10,
        "d_convergence": 1e-10,
    })
    output = tmp_path / "pair-workers.out"
    psi4.core.set_output_file(str(output), False)
    _, reference = psi4.energy("scf", molecule=molecule, return_wfn=True)
    energy, result = psi4.energy("mp2-f12", ref_wfn=reference, return_wfn=True)
    expected = [energy] + [result.variable(key) for key in COMPONENTS]
    psi4.core.flush_outfile()
    assert "pair workers: 4 of 4" in output.read_text()

    factory = psi4.core.f12
    memory = psi4.core.get_memory()

    def small_budget(reference):
        result = factory(reference)
        # SCF and MP2 retain their normal budget; only F12 is constrained.
        psi4.core.set_memory_bytes(1024 * 1024, quiet=True)
        return result

    with monkeypatch.context() as patch:
        patch.setattr(psi4.core, "f12", small_budget)
        try:
            energy, result = psi4.energy("mp2-f12", ref_wfn=reference, return_wfn=True)
        finally:
            psi4.core.set_memory_bytes(memory, quiet=True)
    psi4.core.flush_outfile()
    assert "pair workers: 2 of 4" in output.read_text()
    assert psi4.core.get_num_threads() == 4
    values = [energy] + [result.variable(key) for key in COMPONENTS]
    assert values == pytest.approx(expected, abs=1e-9, rel=0)

    psi4.set_options({"f12_subtype": "incore"})
    energy, result = psi4.energy("mp2-f12", ref_wfn=reference, return_wfn=True)
    values = [energy] + [result.variable(key) for key in COMPONENTS]
    assert values == pytest.approx(expected, abs=1e-9, rel=0)
