import pytest
import psi4
import numpy as np

from utils import compare_values
from addons import using, uusing

pytestmark = [pytest.mark.psi, pytest.mark.api]

# Define test molecules
_fd0_molecules = [
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

# Reference frequencies (in cm^-1)
ref_freqs = [np.array([2115.02012832, 4096.00140028, 4390.01690838]),
    np.array([1234.75439971, 2041.99743682, 2042.07941254, 4016.0550554,
 4304.81615418,4304.85256983]),
    np.array([
        1681.87436846, 1681.90038394, 1681.92651734,
        1905.63629032, 1905.66591125, 3502.51332279,
        3761.53141463, 3761.53877669, 3761.55470204
    ])
]
ref_intensities = [np.array([ 8.68661, 39.02869, 23.97606]),
np.array([107.59281,   2.62096,   2.62098,   5.22301,   2.32676,   2.32678]),
np.array([6.30806, 6.30806, 6.30805, 0.,      0.,      0.,      0.55949, 0.55949, 0.5595 ])]

@pytest.mark.findif
@pytest.mark.parametrize("engine", [pytest.param(False), pytest.param(True, marks=using("molsym"))])
def test_fd_energy_engine(engine):
    print(f"\nTesting with MolSym = {engine}")
    for m, mol in enumerate(_fd0_molecules):
        molecule = psi4.geometry(mol)

        psi4.set_options({
            'basis': 'sto-3g',
            'scf_type': 'pk',
            'molsym': engine,
            'disp_size': 0.005,
            'e_convergence': 1e-11,
            'd_convergence': 1e-11,
            'fd_project': True,
        })

        e, wfn = psi4.frequency('scf', dertype='energy', return_wfn=True)

        vib_deg_free = molecule.natom() * 3 - 6
        frequencies = np.real(wfn.frequency_analysis["omega"].data[-vib_deg_free:])
        intensities = np.real(wfn.frequency_analysis["IR_intensity"].data[-vib_deg_free:])

        print("Computed frequencies (cm^-1):", np.round(frequencies, 1))
        print("Reference frequencies (cm^-1):", np.round(ref_freqs[m], 1))

        # Compare to 0.1 cm^-1 tolerance
        assert np.allclose(frequencies, ref_freqs[m], atol=1e-1), (
            f"Frequencies differ for molecule {m} (MolSym={engine}).\n"
            f"Expected: {ref_freqs[m]}\nGot: {results}"
        )

        assert np.allclose(intensities, ref_intensities[m], atol=1e-1), (
            f"Frequencies differ for molecule {m} (MolSym={engine}).\n"
            f"Expected: {ref_freqs[m]}\nGot: {results}"
        )
