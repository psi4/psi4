import pytest
import numpy as np
from utils import *
from addons import using, uusing

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api]

@uusing("dftd3")
def test_dftd3_dft_grad_lr3():
    """modified VV10-less B97 functional gradient wB97X-V -> wB97X-D3BJ"""

    # stored data from finite differences
    FD_wb97x_d3 = psi4.core.Matrix.from_list([
       [  0.03637802044642,    0.06718963272193,    0.00000000000000],
       [  0.04955519892514,   -0.06340333481039,    0.00000000000000],
       [ -0.07009043821383,   -0.00834477190196,    0.00000000000000],
       [  0.02732425404378,   -0.05883094637658,    0.00000000000000],
       [ -0.02158351760075,    0.03169471018350,    0.05342791683461],
       [ -0.02158351760075,    0.03169471018350,   -0.05342791683461]])

    psi4.geometry("""
    0 1
    O         -1.65542       -0.12330        0.00000
    O          1.24621        0.10269        0.00000
    H         -0.70409        0.03193        0.00000
    H         -2.03867        0.75372        0.00000
    H          1.57599       -0.38252       -0.75856
    H          1.57599       -0.38252        0.75856
    """)

    psi4.set_options({
        'scf_type': 'pk',
        'basis': 'minix',
        'dft_radial_points': 99,
        'dft_spherical_points': 302,
        'e_convergence': 8,
    })

    analytic = psi4.gradient('wB97X-D3BJ', dertype=1)
    assert compare_matrices(analytic, FD_wb97x_d3, 5, "wB97X-D3BJ Analytic vs Store")


# Orca about 400,770 ! HF-3c NoRI TightSCF angs %METHOD AngularGrid 7 IntAcc 29.0 GridPruning Unpruned end
s16_hf3c_xc_orca = np.array([-153.79456839, -77.49902272, -76.29481930])
s16_hf3c_disp_orca = np.array([-0.024367755, -0.014320824, -0.008066062])
# grep 'gCP+bas correction'
s16_hf3c_gcp_orca = np.array([-0.016246604, -0.001871046, -0.014680228])
# grep 'Total correction to HF/MINIX'
s16_hf3c_gcp_orca = np.array([-0.040614359, -0.016191870, -0.022746290])
s16_hf3c_final_orca = np.array([-153.835182748371, -77.515214592777, -76.317565590919])

# > s-dftd3 s16di.tmol --bj-param 1.0 0.8777 0.4171 2.9149
s16_hf3c_sdftd3 = np.array([-2.4367755433608E-02, -1.4320824084920E-02, -8.0660621301533E-03])

# Orca about 400,770 ! PBEh-3c NoRI TightSCF angs %METHOD AngularGrid 7 IntAcc 29.0 GridPruning Unpruned end
# grep 'Total Energy'
s16_pbeh3c_xc_orca = np.array([-155.55536461, -78.40881605, -77.14361375])
# grep 'Dispersion correction'
s16_pbeh3c_disp_orca = np.array([-0.002717357, -0.001301388 -0.000698117])
# grep 'gCP correction'
s16_pbeh3c_gcp_orca = np.array([0.007221336, 0.004956005, 0.001998163])
# grep 'FINAL SINGLE POINT ENERGY'
s16_pbeh3c_final_orca = np.array([-155.550860629116, -78.405161432347, -77.142313701120])

# > s-dftd3 s16di.tmol --bj-param 1.0 0.0 0.4860 4.5000
s16_pbeh3c_sdftd3 = np.array([-2.7171815555163E-03, -1.3016040877044E-03, -6.9806156762005E-04])

# Orca about 400,770 ! B97-3c NoRI TightSCF angs %METHOD AngularGrid 7 IntAcc 29.0 GridPruning Unpruned end
s16_b973c_xc_orca = np.array([-155.79596752, -78.52480800, -77.26992220])
s16_b973c_disp_orca = np.array([-0.010633212, -0.005635464, -0.003247342])
s16_b973c_gcp_orca = np.array([-0.027189940, -0.014493561, -0.012690835])
s16_b973c_final_orca = np.array([-155.833790674968, -78.544937023820, -77.285860375503])

# > s-dftd3 s16di.tmol --bj-param 1.0 1.50 0.37 4.10
s16_b973c_sdftd3 = np.array([-1.0633036221941E-02, -5.6356801553795E-03, -3.2472858339059E-03])

# Orca about 400,770 ! R2SCAN NoRI def2-mTZVPP angs %METHOD AngularGrid 7 IntAcc 29.0 GridPruning Unpruned end (note not tight)
s16_r2scan_final_orca = np.array([-155.86137806, -78.55690019, -77.30233380])
# > dftd4 s16di.tmol --param 1.0 0.0 0.42 5.65 --mbdscale 2.0 --zeta 2.0 1.0
s16_r2scan3c_dftd4 = np.array([-1.4455457680017e-03, -6.2041477166000e-04, -3.3971225644558e-04])
# > mctc-gcp s16di.tmol -l r2scan3
s16_r2scan3c_mctcgcp = np.array([0.0020969558, 0.0016932666, 0.0003797674])

# Orca about 400,770 ! wB97X-D4 NoRI TightSCF angs %METHOD AngularGrid 7 IntAcc 29.0 GridPruning Unpruned
#   D4A1 0.2464 D4A2 4.737 D4S6 1.0 D4S8 0.0 D4S9 1.0 end (plus vDZP user-defined basis)
s16_wb97x3c_xc_orca = np.array([-26.27034106, -13.76904046, -12.49962806])
s16_wb97x3c_disp_orca = np.array([-0.004825743, -0.002498220, -0.001432302])
s16_wb97x3c_final_orca = np.array([-26.275166807430, -13.771538675222, -12.501060363754])

# wb97x-d3 from suppmat p35 of https://doi.org/10.1063/5.0133026
# > dftd4 s16di.tmol --param 1.0 0.0 0.2464 4.737 --mbdscale 1.0
s16_wb97x3c_dftd4 = np.array([-4.8257433806220e-03, -2.4982196402379e-03, -1.4323024694511e-03])

#d4.bj-eeq-atm = { s8=0.60187490, a1=0.51559235, a2=5.77342911, doi="10.1063/5.0041008" }
# > dftd4 s16ma.tmol --param 1.0 0.60187490 0.51559235 5.77342911 --mbdscale 1.0
s16_r2scan_dftd4 = np.array([-1.1296999851733E-03, -4.4567923427274E-04, -2.4728471580494E-04])

# > dftd4 s16di.tmol --param 1.0 0.8324 0.4944 5.9019 --mbdscale 1.0
#s16_r2scanh_dftd4 = np.array([

#[parameter.r2scan0]
# > dftd4 s16di.tmol --param 1.0 0.8992 0.4778 5.8779 --mbdscale 1.0
#s16_r2scan0_dftd4 = np.array([

# > dftd4 s16di.tmol --param 1.0 1.0471 0.4574 5.8969 --mbdscale 1.0
s16_r2scan50_dftd4 = np.array([-1.4008964960682E-03, -5.5551918498700E-04, -3.1465561625320E-04])

# refs from psi4, not external
s16_r2scan0_psi4 = np.array([-155.6718657, -78.4640942, -77.2052683])
s16_r2scanh_psi4 = np.array([-155.6812544, -78.4677956, -77.2110203])
s16_r2scan50_psi4 = np.array([-155.6575747, -78.4586347, -77.1962996])


@pytest.mark.nbody
@pytest.mark.parametrize("mode", [
    pytest.param("abs", marks=pytest.mark.long),
    pytest.param("ie", marks=pytest.mark.quick),
])
@pytest.mark.parametrize("mtdbas,ref", [
    pytest.param("hf-3c/", s16_hf3c_final_orca, marks=[*using("dftd3"), *using("gcp")]),  # MINIX
    pytest.param("pbeh-3c/", s16_pbeh3c_final_orca, marks=[*using("dftd3"), *using("gcp")]),  # def2-mSVP
    pytest.param("b97-3c/", s16_b973c_final_orca, marks=[*using("dftd3"), *using("mctc-gcp")]),  # def2-mTZVP
    pytest.param("r2scan0/def2-svp", s16_r2scan0_psi4),
    pytest.param("r2scanh/def2-svp", s16_r2scanh_psi4),
    pytest.param("r2scan50-d4/def2-svp", s16_r2scan50_psi4 + s16_r2scan50_dftd4, marks=using("dftd4")),
    pytest.param("r2scan50/def2-svp", s16_r2scan50_psi4),
    pytest.param("r2scan-d4/def2-mTZVPP", s16_r2scan_final_orca + s16_r2scan_dftd4, marks=using("dftd4")),
    pytest.param("r2scan/def2-mTZVPP", s16_r2scan_final_orca),
    pytest.param("r2scan-3c/", s16_r2scan_final_orca + s16_r2scan3c_dftd4 + s16_r2scan3c_mctcgcp, marks=[*using("dftd4_350"), *using("mctc-gcp")]),  # def2-mTZVPP
    pytest.param("wb97x-3c/", s16_wb97x3c_xc_orca + s16_wb97x3c_dftd4, marks=using("dftd4")),  # vDZP
])
def test_grimme_3c(mtdbas, ref, mode):

    s16di = psi4.geometry("""
    # ang C   0.000000  -0.667578  -2.124659
    # ang C   0.000000   0.667578  -2.124659
    # ang H   0.923621  -1.232253  -2.126185
    # ang H  -0.923621  -1.232253  -2.126185
    # ang H  -0.923621   1.232253  -2.126185
    # ang H   0.923621   1.232253  -2.126185
    # ang --
    # ang C   0.000000   0.000000   2.900503
    # ang C   0.000000   0.000000   1.693240
    # ang H   0.000000   0.000000   0.627352
    # ang H   0.000000   0.000000   3.963929
    units au
    C   0.00000000 -1.26153959 -4.01502362
    C   0.00000000  1.26153959 -4.01502362
    H   1.74539073 -2.32862069 -4.01790734
    H  -1.74539073 -2.32862069 -4.01790734
    H  -1.74539073  2.32862069 -4.01790734
    H   1.74539073  2.32862069 -4.01790734
    --
    C   0.00000000  0.00000000  5.48115630
    C   0.00000000  0.00000000  3.19975986
    H   0.00000000  0.00000000  1.18552346
    H   0.00000000  0.00000000  7.49074019
    symmetry c1
    """)
    kcal = psi4.driver.constants.hartree2kcalmol

    if mode == "abs":
        psi4.set_options({
            "scf_type": "direct",
            "dft_radial_points": 300,
            "dft_spherical_points": 770,
            "e_convergence": 8,
            "d_convergence": 8,
        })

        ene = psi4.energy(mtdbas)
        assert psi4.compare_values(ref[0], ene, 5, mtdbas)

    if mode == "ie":
        psi4.set_options({
            "scf_type": "pk",
            "dft_radial_points": 99,
            "dft_spherical_points": 302,
            "e_convergence": 8,
            "d_convergence": 8,
        })
        psi4.set_options({'basis': 'cc-pvdz'})  # try to confuse method
        ene = psi4.energy(mtdbas, bsse_type='nocp')
        assert psi4.compare_values(kcal * (ref[0] - ref[1] - ref[2]), kcal * ene, 1.1e-3, mtdbas)

