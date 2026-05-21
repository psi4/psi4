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



_fd_ref_freqs_energy = [np.array([2115.02012832, 4096.00140028, 4390.01690838]),
    np.array([1234.75439971, 2041.99743682, 2042.07941254, 4016.0550554,
 4304.81615418,4304.85256983]),
    np.array([
        1681.87436846, 1681.90038394, 1681.92651734,
        1905.63629032, 1905.66591125, 3502.51332279,
        3761.53141463, 3761.53877669, 3761.55470204
    ])
]
_fd_ref_intensities_energy = [np.array([ 8.68661, 39.02869, 23.97606]),
np.array([107.59281,   2.62096,   2.62098,   5.22301,   2.32676,   2.32678]),
np.array([6.30806, 6.30806, 6.30805, 0.,      0.,      0.,      0.55949, 0.55949, 0.5595 ])]

_fd_ref_freqs_gradient = [np.array([2115.0425274, 4095.99314115, 4390.02370649]),
np.array([1234.75063586, 2042.0645856, 2042.0698831, 4016.05557387, 4304.82565468, 4304.83414845]),
    np.array([1681.91158886, 1681.91976851, 1681.91976851, 1905.658396,
 1905.6663143, 3502.50350807, 3761.53410271, 3761.53949881, 3761.53949881])
]
_fd_ref_intensities_gradient = [
np.array([8.68657198, 39.02873184, 23.97605778]),
np.array([107.59260692,   2.6209917,    2.62099019,   5.22323669,   2.32675017, 2.32677372]),
np.array([6.30809029e+00, 6.30807464e+00, 6.30807464e+00, 1.23772972e-23,
 2.74201817e-30, 1.02555259e-20, 5.59524983e-01, 5.59450528e-01, 5.59450528e-01])]


@pytest.mark.parametrize(
    "engine",
    [
        pytest.param(False),
        pytest.param(True, marks=using("molsym")),
    ]
)
@pytest.mark.parametrize("exploit_degeneracy", [False, True])
@pytest.mark.parametrize(
    "dertype, ref_freqs, ref_intensities",
    [
        ("energy", _fd_ref_freqs_energy, _fd_ref_intensities_energy),
        ("gradient", _fd_ref_freqs_gradient, _fd_ref_intensities_gradient),
    ]
)
def test_molsym_findif(engine, exploit_degeneracy,
                       dertype, ref_freqs, ref_intensities):

    for m, mol in enumerate(_fd_molecules):

        molecule = psi4.geometry(mol)

        psi4.set_options({
            'basis': 'sto-3g',
            'scf_type': 'pk',
            'molsym': engine,
            'disp_size': 0.005,
            'e_convergence': 1e-11,
            'd_convergence': 1e-11,
            'fd_project': True,
            'exploit_degeneracy': exploit_degeneracy,
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
            ref_freqs[m],
            atol=1e-1
        )

        assert np.allclose(
            intensities,
            ref_intensities[m],
            atol=1e-1
        )
