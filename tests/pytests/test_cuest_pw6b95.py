import numpy as np
import psi4
import pytest

from addons import uusing


pytestmark = [pytest.mark.psi, pytest.mark.api]


def _format_gradient(gradient):
    return np.array2string(gradient, precision=10, separator=", ", floatmode="fixed")


@uusing("cuest")
@uusing("cuda_cc8")
def test_cuest_pw6b95_mixed_rung_gradient():
    """Run with ``pytest -s`` to print values for migration into test_cuest.py."""
    psi4.core.set_num_threads(4)
    molecule = psi4.geometry(
        """
        units bohr
        N   -1.2443656662    1.5116296128    1.1401094834
        C    1.2072325506    0.2976885350    0.9288948097
        H    2.1975245576    0.3690634856    2.7411987245
        H    2.3661826021    1.2428916800   -0.4976593391
        H    0.9782544525   -1.6841238010    0.3901528276
        H   -2.1240520112    1.4799955998   -0.5728893307
        H   -1.0022350753    3.3701184308    1.5836092757
        """
    )

    options = {
        "basis": "def2-svp",
        "scf_type": "df",
        "dft_nuclear_scheme": "stratmann",
        "dft_spherical_points": 590,
        "dft_radial_points": 100,
        "df_basis_scf": "def2-universal-JKFIT",
        "maxiter": 300,
        "e_convergence": 10,
        "d_convergence": 9,
        "puream": True,
        "reference": "rhf",
        "cuest_mixed_precision": False,
    }

    psi4.set_options({**options, "use_cuest": False})
    psi4_gradient, psi4_wfn = psi4.gradient("pw6b95", return_wfn=True, molecule=molecule)
    psi4_energy = psi4_wfn.get_energies("Total Energy")
    psi4_gradient = psi4_gradient.np.copy()

    psi4.core.clean()
    psi4.set_options({**options, "use_cuest": True})
    cuest_gradient, cuest_wfn = psi4.gradient("pw6b95", return_wfn=True, molecule=molecule)
    cuest_energy = cuest_wfn.get_energies("Total Energy")
    cuest_gradient = cuest_gradient.np.copy()

    print("\nPW6B95 values for test_cuest.py")
    print(f'Psi4 "ref_energy": {psi4_energy:.15f},')
    print(f'Psi4 "ref_gradient": {_format_gradient(psi4_gradient)}')
    print(f'cuEST "energy": {cuest_energy:.15f},')
    print(f'cuEST "gradient": {_format_gradient(cuest_gradient)}')
    print(f"energy difference: {cuest_energy - psi4_energy:.15e}")
    print(f"maximum gradient difference: {np.max(np.abs(cuest_gradient - psi4_gradient)):.15e}")

    assert psi4.compare_values(
        psi4_energy,
        cuest_energy,
        atol=5.0e-6,
        label="PW6B95 Psi4 DF versus cuEST DF energy",
    )
    assert psi4.compare_values(
        psi4_gradient,
        cuest_gradient,
        atol=5.0e-5,
        label="PW6B95 Psi4 DF versus cuEST DF gradient",
    )
