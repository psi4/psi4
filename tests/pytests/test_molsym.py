import pytest
import psi4
import numpy as np

from addons import using

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.findif]

_fd_molecules = [
    # Water
    """
    O        1.2091536548      1.7664118189     -0.0171613972
    H        2.1984800075      1.7977100627      0.0121161719
    H        0.9197881882      2.4580185570      0.6297938832
    """,
    # Ammonia
    """
    N  0.          0.          0.07056746
    H  0.46649474 -0.80799259 -0.32682968
    H  0.46649474  0.80799259 -0.32682968
    H -0.93298949  0.         -0.32682968
    """,
    # Methane
    """
    C       0.0000000       0.0000000       0.0000000
    H       0.6268910       0.6268910       0.6268910
    H       -0.6268910      -0.6268910      0.6268910
    H       -0.6268910      0.6268910       -0.6268910
    H       0.6268910       -0.6268910      -0.6268910
    """
]

_fd_ref_frequencies = [np.array([2115.0425274, 4095.99314115, 4390.02370649]),
np.array([1234.75063586, 2042.0645856, 2042.0698831, 4016.05557387, 4304.82565468, 4304.83414845]),
    np.array([1681.91158886, 1681.91976851, 1681.91976851, 1905.658396,
 1905.6663143, 3502.50350807, 3761.53410271, 3761.53949881, 3761.53949881])
]
_fd_ref_intensities = [
np.array([8.68657198, 39.02873184, 23.97605778]),
np.array([107.59260692,   2.6209917,    2.62099019,   5.22323669,   2.32675017, 2.32677372]),
np.array([6.30809029e+00, 6.30807464e+00, 6.30807464e+00, 1.23772972e-23,
 2.74201817e-30, 1.02555259e-20, 5.59524983e-01, 5.59450528e-01, 5.59450528e-01])]


@pytest.mark.parametrize("salc_package", 
    [
        pytest.param("internal"),
        pytest.param("molsym", marks=using("molsym")),
    ]
)
@pytest.mark.parametrize("molsym_exploit_degeneracy", [False, True])
@pytest.mark.parametrize("dertype", ["energy","gradient"])
def test_molsym_findif(salc_package, molsym_exploit_degeneracy,dertype):
    for m, mol in enumerate(_fd_molecules):
        molecule = psi4.geometry(mol)
        psi4.set_options({
            'basis': 'sto-3g',
            'scf_type': 'pk',
            'salc_package': salc_package,
            'molsym_exploit_degeneracy': molsym_exploit_degeneracy,
            'disp_size': 0.005,
            'e_convergence': 1e-11,
            'd_convergence': 1e-11,
            'fd_project': True,
        })

        e, wfn = psi4.frequency(
            'scf',
            dertype=dertype,
            return_wfn=True
        )

        vib_deg_free = molecule.natom() * 3 - 6

        frequencies = np.real(
            wfn.frequency_analysis["omega"].data[-vib_deg_free:]
        )

        intensities = np.real(
            wfn.frequency_analysis["IR_intensity"].data[-vib_deg_free:]
        )

        assert np.allclose(
            frequencies,
            _fd_ref_frequencies[m],
            atol=1e-1
        )

        assert np.allclose(
            intensities,
            _fd_ref_intensities[m],
            atol=1e-1
        )
