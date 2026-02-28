import pytest

import psi4

from addons import uusing


@uusing("gauxc")
@pytest.mark.parametrize("inp", [
    pytest.param({'name': 'svwn'}, id='svwn'),
    pytest.param({'name': 'pbe'}, id='pbe'),
    pytest.param({'name': 'b3lyp'}, id='b3lyp'),
    pytest.param({'name': 'wb97x'}, id='wb97x'),
    pytest.param({'name': 'b2plyp'}, id='b2plyp')
])
@pytest.mark.parametrize("symmetry", [
    pytest.param({'on': True}, id='sym-true'),
    pytest.param({'off': False}, id='sym-false', marks=pytest.mark.xfail),
])
@pytest.mark.parametrize("basis", [
    pytest.param({'name': "sto-6g"}, id='sto6g'),
    pytest.param({'name': "cc-pvdz"}, id='dz'),
])
def test_dft(inp, symmetry, basis):
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
    """)

    if symmetry["on"]: h2o.reset_point_group("c1")

    psi4.set_options({
        "gauxc_integrate": False,
        "basis": basis["name"],
        "d_convergence": 10,
        "dft_radial_points": 80,
        "dft_spherical_points": 590
    })

    enPsi = psi4.energy(inp["name"])

    psi4.set_options({"gauxc_integrate": True, "gauxc_radial_points": 80, "gauxc_spherical_points": 590})

    enGau = psi4.energy(inp["name"])

    assert psi4.compare_values(enPsi, enGau, 6, f"{inp['name']} energies")
