import pytest

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def _dipole_polarizability(screening):
    """Static dipole polarizability of water with SCF_TYPE DIRECT under the given SCREENING."""
    psi4.core.clean()
    psi4.core.clean_variables()
    psi4.core.clean_options()
    mol = psi4.geometry("""
        O   0.000000000000     0.000000000000    -0.075367841241
        H   0.000000000000    -0.869758531118     0.598071166797
        H   0.000000000000     0.869758531118     0.598071166797
    """)
    psi4.set_options({"basis": "cc-pVDZ", "scf_type": "direct", "screening": screening})
    psi4.properties("scf", properties=["DIPOLE_POLARIZABILITIES"], molecule=mol)
    return [psi4.variable(f"DIPOLE POLARIZABILITY {c}") for c in ("XX", "YY", "ZZ")]


def test_density_screening_many_densities():
    """SCREENING=DENSITY must bound every density handed to the JK object.

    The CPHF solve passes one density per perturbation, three here. A bound
    built from the first two only screens away the quartets that carry the
    third (z) response, so alpha_zz comes out wrong."""
    ref = _dipole_polarizability("csam")
    dens = _dipole_polarizability("density")
    for r, d, c in zip(ref, dens, ("XX", "YY", "ZZ")):
        assert psi4.compare_values(r, d, 6, f"alpha_{c}: DENSITY vs CSAM screening")
