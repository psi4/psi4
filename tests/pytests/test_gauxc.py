from pathlib import Path

import pytest

import psi4

from addons import using


GAUXC_EXEC_PARAMS = [
    pytest.param({"use_gpu": False, "label": "Host (CPU)"}, id="cpu", marks=using("gauxc")),
    pytest.param({"use_gpu": True, "label": "Device (GPU)"}, id="gpu", marks=[*using("gauxc_gpu"), *using("cuda")]),
]


def assert_gauxc_execution_space(gauxc_exec):
    if hasattr(psi4.core, "flush_outfile"):
        psi4.core.flush_outfile()
    output = Path("pytest_output.dat").read_text(errors="replace")
    assert f"Execution Space        = {gauxc_exec['label']:>14}" in output


@pytest.mark.parametrize("gauxc_exec", GAUXC_EXEC_PARAMS)
@pytest.mark.parametrize("inp", [
    pytest.param({'name': 'svwn'}, id='svwn'),
    pytest.param({'name': 'pbe'}, id='pbe'),
    pytest.param({'name': 'b3lyp'}, id='b3lyp'),
    pytest.param({'name': 'wb97x'}, id='wb97x'),
    pytest.param({'name': 'b2plyp'}, id='b2plyp')
])
@pytest.mark.parametrize("symmetry", [
    pytest.param({'on': True}, id='sym-true'),
    pytest.param({'on': False}, id='sym-false'),
])
@pytest.mark.parametrize("basis", [
    pytest.param({'name': "sto-6g"}, id='sto6g'),
    pytest.param({'name': "cc-pvdz"}, id='dz'),
])
def test_rks_energy(inp, symmetry, basis, gauxc_exec):
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
    """)

    if not symmetry["on"]: h2o.reset_point_group("c1")

    psi4.set_options({
        "gauxc_dft_enable": False,
        "basis": basis["name"],
        "d_convergence": 10,
        "dft_radial_points": 80,
        "dft_spherical_points": 590
    })

    enPsi = psi4.energy(inp["name"])

    psi4.set_options({
        "gauxc_dft_enable": True,
        "gauxc_use_gpu": gauxc_exec["use_gpu"],
        "gauxc_radial_points": 80,
        "gauxc_spherical_points": 590,
        "dft_enable_psi": False,
    })

    enGau = psi4.energy(inp["name"])
    assert_gauxc_execution_space(gauxc_exec)

    assert psi4.compare_values(enPsi, enGau, 6, f"{inp['name']} energies")


@pytest.mark.parametrize("gauxc_exec", GAUXC_EXEC_PARAMS)
@pytest.mark.parametrize("inp", [
    pytest.param({'name': 'svwn'}, id='svwn'),
    pytest.param({'name': 'pbe'}, id='pbe'),
    pytest.param({'name': 'b3lyp'}, id='b3lyp'),
    pytest.param({'name': 'wb97x'}, id='wb97x'),
    #pytest.param({'name': 'b2plyp'}, id='b2plyp')
])
@pytest.mark.parametrize("symmetry", [
    pytest.param({'on': True}, id='sym-true'),
    pytest.param({'on': False}, id='sym-false'),
])
@pytest.mark.parametrize("basis", [
    pytest.param({'name': "sto-6g"}, id='sto6g'),
    pytest.param({'name': "cc-pvdz"}, id='dz'),
])
def test_rks_grad(inp, symmetry, basis, gauxc_exec):
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
    """)

    if not symmetry["on"]: h2o.reset_point_group("c1")

    psi4.set_options({
        "gauxc_dft_enable": False,
        "basis": basis["name"],
        "d_convergence": 10,
        "dft_radial_points": 80,
        "dft_spherical_points": 590
    })

    psi4.set_options({
        "gauxc_dft_enable": True,
        "gauxc_use_gpu": gauxc_exec["use_gpu"],
        "gauxc_radial_points": 80,
        "gauxc_spherical_points": 590,
        "dft_enable_psi": False,
    })

    findif_gradient = psi4.gradient(inp["name"], dertype=0)
    analytic_gradient = psi4.gradient(inp["name"], dertype=1)
    assert_gauxc_execution_space(gauxc_exec)

    assert psi4.compare_values(findif_gradient, analytic_gradient, 5, "analytic vs. findif gradient")


@pytest.mark.parametrize("gauxc_exec", GAUXC_EXEC_PARAMS)
@pytest.mark.parametrize("inp", [
    pytest.param({'name': 'svwn'}, id='svwn'),
    pytest.param({'name': 'pbe'}, id='pbe'),
    pytest.param({'name': 'b3lyp'}, id='b3lyp'),
    pytest.param({'name': 'wb97x'}, id='wb97x'),
    #pytest.param({'name': 'b2plyp'}, id='b2plyp') # NotYetImplemented
])
@pytest.mark.parametrize("symmetry", [
    pytest.param({'on': True}, id='sym-true'),
    pytest.param({'on': False}, id='sym-false'),
])
@pytest.mark.parametrize("tda", [
    pytest.param({'on': True}, id='tda-true'),
    pytest.param({'on': False}, id='tda-false'),
])
@pytest.mark.parametrize("basis", [
    pytest.param({'name': "sto-6g"}, id='sto6g'),
    pytest.param({'name': "cc-pvdz"}, id='dz'),
])
def test_rks_tdscf(inp, symmetry, tda, basis, gauxc_exec):
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
    """)

    if not symmetry["on"]: h2o.reset_point_group("c1")

    psi4.set_options({
        "gauxc_dft_enable": False,
        "basis": basis["name"],
        "d_convergence": 10,
        "dft_radial_points": 80,
        "dft_spherical_points": 590,
        "tdscf_tda": tda["on"],
        "tdscf_states": [2, 1, 1, 2] if symmetry["on"] else [6],
        # LRC df integrals are non-hermitian due to the use of a
        # non-erf metric, so we need pk integrals for tight convergence.
        "scf_type": "pk" if "wb97x" == inp["name"] else "df"
    })

    psi_wfn = psi4.energy(inp["name"], return_wfn=True)[1]
    psi4.tdscf(psi_wfn)

    psi4.set_options({
        "gauxc_dft_enable": True,
        "gauxc_use_gpu": gauxc_exec["use_gpu"],
        "gauxc_radial_points": 80,
        "gauxc_spherical_points": 590,
        "dft_enable_psi": False,
    })

    gau_wfn = psi4.energy(inp["name"], return_wfn=True)[1]
    psi4.tdscf(gau_wfn)

    for n in range(6):
        ref_en = psi_wfn.variable(f"TD-{inp['name']} ROOT 0 -> ROOT {n+1} EXCITATION ENERGY")
        gau_en = gau_wfn.variable(f"TD-{inp['name']} ROOT 0 -> ROOT {n+1} EXCITATION ENERGY")
        assert psi4.compare_values(ref_en, gau_en, 4, f"TD-{inp['name']}/cc-pvdz ROOT 0 -> ROOT {n+1} EXCITATION ENERGY")


@pytest.mark.parametrize("gauxc_exec", GAUXC_EXEC_PARAMS)
@pytest.mark.parametrize("inp", [
    pytest.param({'name': 'svwn'}, id='svwn'),
    pytest.param({'name': 'pbe'}, id='pbe'),
    pytest.param({'name': 'b3lyp'}, id='b3lyp'),
    pytest.param({'name': 'wb97x'}, id='wb97x'),
    pytest.param({'name': 'b2plyp'}, id='b2plyp')
])
@pytest.mark.parametrize("symmetry", [
    pytest.param({'on': True}, id='sym-true'),
    pytest.param({'on': False}, id='sym-false'),
])
@pytest.mark.parametrize("basis", [
    pytest.param({'name': "sto-6g"}, id='sto6g'),
    pytest.param({'name': "cc-pvdz"}, id='dz'),
])
def test_uks_energy(inp, symmetry, basis, gauxc_exec):
    h2o = psi4.geometry("""
        1 2
        O
        H 1 1.0
        H 1 1.0 2 101.5
    """)

    if not symmetry["on"]: h2o.reset_point_group("c1")

    psi4.set_options({
        "gauxc_dft_enable": False,
        "basis": basis["name"],
        "d_convergence": 10,
        "dft_radial_points": 80,
        "dft_spherical_points": 590,
        "reference": "uks"
    })

    enPsi = psi4.energy(inp["name"])

    psi4.set_options({
        "gauxc_dft_enable": True,
        "gauxc_use_gpu": gauxc_exec["use_gpu"],
        "gauxc_radial_points": 80,
        "gauxc_spherical_points": 590,
        "dft_enable_psi": False,
    })

    enGau = psi4.energy(inp["name"])
    assert_gauxc_execution_space(gauxc_exec)

    assert psi4.compare_values(enPsi, enGau, 6, f"{inp['name']} energies")


@pytest.mark.parametrize("gauxc_exec", GAUXC_EXEC_PARAMS)
@pytest.mark.parametrize("inp", [
    pytest.param({'name': 'svwn'}, id='svwn'),
    pytest.param({'name': 'pbe'}, id='pbe'),
    pytest.param({'name': 'b3lyp'}, id='b3lyp'),
    pytest.param({'name': 'wb97x'}, id='wb97x'),
    #pytest.param({'name': 'b2plyp'}, id='b2plyp')
])
@pytest.mark.parametrize("symmetry", [
    pytest.param({'on': True}, id='sym-true'),
    pytest.param({'on': False}, id='sym-false'),
])
@pytest.mark.parametrize("basis", [
    pytest.param({'name': "sto-6g"}, id='sto6g'),
    pytest.param({'name': "cc-pvdz"}, id='dz'),
])
def test_uks_gradient(inp, symmetry, basis, gauxc_exec):
    h2o = psi4.geometry("""
        1 2
        O
        H 1 1.0
        H 1 1.0 2 101.5
    """)

    if not symmetry["on"]: h2o.reset_point_group("c1")

    psi4.set_options({
        "gauxc_dft_enable": False,
        "basis": basis["name"],
        "d_convergence": 10,
        "dft_radial_points": 80,
        "dft_spherical_points": 590,
        "reference": "uks"
    })

    psi4.set_options({
        "gauxc_dft_enable": True,
        "gauxc_use_gpu": gauxc_exec["use_gpu"],
        "gauxc_radial_points": 80,
        "gauxc_spherical_points": 590,
        "dft_enable_psi": False,
    })

    findif_gradient = psi4.gradient(inp["name"], dertype=0)
    analytic_gradient = psi4.gradient(inp["name"], dertype=1)
    assert_gauxc_execution_space(gauxc_exec)

    assert psi4.compare_values(findif_gradient, analytic_gradient, 5, "analytic vs. findif gradient")


@pytest.mark.parametrize("gauxc_exec", GAUXC_EXEC_PARAMS)
@pytest.mark.parametrize("inp", [
    pytest.param({'name': 'svwn'}, id='svwn'),
    pytest.param({'name': 'pbe'}, id='pbe'),
    pytest.param({'name': 'b3lyp'}, id='b3lyp'),
    pytest.param({'name': 'wb97x'}, id='wb97x'),
    #pytest.param({'name': 'b2plyp'}, id='b2plyp') NotYetImplemented
])
@pytest.mark.parametrize("symmetry", [
    pytest.param({'on': True}, id='sym-true'),
    pytest.param({'on': False}, id='sym-false'),
])
@pytest.mark.parametrize("tda", [
    pytest.param({'on': True}, id='tda-true'),
    pytest.param({'on': False}, id='tda-false'),
])
@pytest.mark.parametrize("basis", [
    pytest.param({'name': "sto-6g"}, id='sto6g'),
    pytest.param({'name': "cc-pvdz"}, id='dz'),
])
def test_uks_tdscf(inp, symmetry, tda, basis, gauxc_exec):
    h2o = psi4.geometry("""
        1 2
        O
        H 1 1.0
        H 1 1.0 2 101.5
    """)

    if not symmetry["on"]: h2o.reset_point_group("c1")

    psi4.set_options({
        "gauxc_dft_enable": False,
        "basis": basis["name"],
        "d_convergence": 10,
        "dft_radial_points": 100,
        "dft_spherical_points": 770,
        "reference": "uks",
        "tdscf_tda": tda["on"],
        "tdscf_states": [2, 1, 1, 2] if symmetry["on"] else [6],
        # LRC df integrals are non-hermitian due to the use of a
        # non-erf metric, so we need pk integrals for tight convergence.
        "scf_type": "pk" if "wb97x" == inp["name"] else "df"
    })

    psi_wfn = psi4.energy(inp["name"], return_wfn=True)[1]
    psi4.tdscf(psi_wfn)

    psi4.set_options({
        "gauxc_dft_enable": True,
        "gauxc_use_gpu": gauxc_exec["use_gpu"],
        "gauxc_radial_points": 80,
        "gauxc_spherical_points": 590,
        "dft_enable_psi": False,
    })

    gau_wfn = psi4.energy(inp["name"], return_wfn=True)[1]
    psi4.tdscf(gau_wfn)

    for n in range(6):
        ref_en = psi_wfn.variable(f"TD-{inp['name']} ROOT 0 -> ROOT {n+1} EXCITATION ENERGY")
        gau_en = gau_wfn.variable(f"TD-{inp['name']} ROOT 0 -> ROOT {n+1} EXCITATION ENERGY")
        assert psi4.compare_values(ref_en, gau_en, 4, f"TD-{inp['name']}/cc-pvdz ROOT 0 -> ROOT {n+1} EXCITATION ENERGY")


