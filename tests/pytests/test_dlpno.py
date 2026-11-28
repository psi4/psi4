import pytest

from utils import compare_values

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.dlpno, pytest.mark.quick]


# Two water molecules ~10 Angstrom apart (O--O along z, second monomer rotated 90 deg).
# Every intermolecular LMO pair is far outside the differential-overlap cutoff, so
# whether it is computed explicitly or screened into the dipole correction is
# controlled entirely by T_CUT_PRE.
_water_dimer_10A = """
0 1
O    0.000000    0.000000    0.000000
H    0.757000    0.586000    0.000000
H   -0.757000    0.586000    0.000000
--
0 1
O    0.000000    0.000000   10.000000
H    0.000000    0.757000   10.586000
H    0.000000   -0.757000   10.586000
symmetry c1
"""


def test_dlpno_dipole_pair_correction():
    """Distant (LMO, LMO) pairs screened by T_CUT_PRE are replaced by a dipole
    (SC-MP2) estimate of their pair energy. That estimate must reproduce the
    explicit LMP2 energy of the very same pairs. Before psi4/psi4#3619 every
    screened off-diagonal pair was counted for both (i,j) and (j,i), doubling
    the dipole correction, so the two DLPNO-CCSD energies below disagreed by
    roughly the whole correction.
    """

    dimer = psi4.geometry(_water_dimer_10A)

    common_options = {
        "basis": "cc-pvdz",
        "freeze_core": True,
        "scf_type": "df",
        "d_convergence": 1.0e-10,
        "e_convergence": 1.0e-10,
        "r_convergence": 1.0e-9,
        # keep every surviving distant pair at the iterative LMP2 level
        "t_cut_pairs_mp2": 0.0,
    }

    # (a) every intermolecular pair is screened -> dipole correction only
    psi4.set_options({**common_options, "t_cut_pre": 1.0e-3})
    e_screened = psi4.energy("dlpno-ccsd", molecule=dimer)
    de_dipole = psi4.variable("DLPNO DIPOLE ENERGY")

    psi4.core.clean()

    # (b) no pair-energy prescreening -> the same pairs are computed explicitly
    #     (as LMP2 weak pairs), so no dipole correction is applied
    psi4.set_options({**common_options, "t_cut_pre": 0.0})
    e_explicit = psi4.energy("dlpno-ccsd", molecule=dimer)
    de_dipole_explicit = psi4.variable("DLPNO DIPOLE ENERGY")

    # sanity: (a) really did screen the intermolecular pairs, (b) really did not
    assert de_dipole < -1.0e-8, f"no LMO pairs were screened into the dipole correction ({de_dipole})"
    assert compare_values(0.0, de_dipole_explicit, atol=1.0e-12, label="no dipole pairs with T_CUT_PRE = 0")

    # Dipole estimate vs. explicit treatment of the same pairs. With the
    # double-counting bug the difference is ~50% of the dipole correction;
    # once fixed it is only the (small) multipole-truncation error.
    assert compare_values(e_explicit,
                          e_screened,
                          atol=0.25 * abs(de_dipole),
                          label="dipole-screened vs. explicit distant pairs")
