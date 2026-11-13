import numpy as np
import pytest

from addons import uusing
from utils import compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api]


@uusing("ecpint")
@pytest.mark.parametrize("puream", [False, True], ids=["cartesian", "spherical"])
@pytest.mark.parametrize("basis_names", [
    pytest.param(("def2-svp", "def2-svp"), id="svp-svp"),
    pytest.param(("def2-svp", "def2-tzvp"), id="svp-tzvp"),
    pytest.param(("def2-svp", "def2-qzvp"), id="svp-qzvp"),
    pytest.param(("def2-tzvp", "def2-qzvp"), id="tzvp-qzvp"),
])
def test_ecp_mixed_basis_transpose(basis_names, puream):
    """The same ECP operator must transpose when its bra and ket are exchanged."""
    import psi4

    mol = psi4.geometry("""
        -1 1
        I 0 0 0
        symmetry c1
        no_com
        no_reorient
    """)
    # These def2 bases all carry the same 28-electron ECP for iodine.
    bra, ket = [
        psi4.core.BasisSet.build(mol, "ORBITAL", name, puream=puream)
        for name in basis_names
    ]
    assert bra.n_ecp_core() == ket.n_ecp_core() == 28
    if basis_names[0] != basis_names[1]:
        assert bra.max_am() < ket.max_am()

    mints = psi4.core.MintsHelper(bra)
    square = mints.ao_ecp().np.copy()
    np.testing.assert_allclose(square, square.T, rtol=1e-10, atol=1e-10)
    forward = mints.ao_ecp(bra, ket).np.copy()
    reverse = mints.ao_ecp(ket, bra).np.copy()
    assert forward.shape == (bra.nbf(), ket.nbf())
    assert reverse.shape == (ket.nbf(), bra.nbf())
    assert np.isfinite(forward).all()
    assert np.isfinite(reverse).all()
    np.testing.assert_allclose(forward, reverse.T, rtol=1e-10, atol=1e-10)


@uusing("ecpint")
def test_ecp_ghost_sanity():
    import psi4

    psi4.set_memory("4 GB")
    psi4.geometry(
        """
        units Angstrom
        no_com
        no_reorient
        symmetry c1
        0 1
        Ar      0.000000   0.000000  -1.500000
        Gh(I)   0.000000   0.000000   1.500000
"""
    )

    ### ECP + valence
    psi4.basis_helper(
        """
        assign aug-cc-pvtz
        assign I myecp
        [myecp]
        I     0
        S   11   1.00
            5.546500D+03           1.560000D-04
            8.382140D+02           9.860000D-04
            1.821870D+02           2.792000D-03
            3.121230D+01          -4.325100D-02
            1.953140D+01           2.341340D-01
            8.240990D+00          -7.509430D-01
            2.194550D+00           8.829680D-01
            1.109110D+00           4.620610D-01
            3.746410D-01           2.228600D-02
            1.770800D-01          -4.353000D-03
            8.106100D-02           1.102000D-03
        S   11   1.00
            5.546500D+03          -7.300000D-05
            8.382140D+02          -5.080000D-04
            1.821870D+02          -1.158000D-03
            3.121230D+01           1.219300D-02
            1.953140D+01          -8.785400D-02
            8.240990D+00           3.382000D-01
            2.194550D+00          -5.765500D-01
            1.109110D+00          -4.092980D-01
            3.746410D-01           5.674590D-01
            1.770800D-01           6.124890D-01
            8.106100D-02           1.432310D-01
        S   1   1.00
            3.746410D-01           1.000000D+00
        S   1   1.00
            1.770800D-01           1.000000D+00
        S   1   1.00
            1.001000D-01           1.000000D+00
        S   1   1.00
            4.120000D-02           1.0000000
        P   9   1.00
            1.889880D+02           5.850000D-04
            2.128680D+01           3.692300D-02
            1.003960D+01          -2.353240D-01
            3.451800D+00           3.414830D-01
            1.974560D+00           5.347880D-01
            1.024200D+00           2.651410D-01
            4.494370D-01           2.578700D-02
            1.866480D-01           5.220000D-04
            7.348100D-02           6.060000D-04
        P   9   1.00
            1.889880D+02          -2.560000D-04
            2.128680D+01          -1.168200D-02
            1.003960D+01           8.319200D-02
            3.451800D+00          -1.569700D-01
            1.974560D+00          -2.245180D-01
            1.024200D+00          -1.144510D-01
            4.494370D-01           3.753560D-01
            1.866480D-01           5.751360D-01
            7.348100D-02           2.459170D-01
        P   1   1.00
            5.981000D-01           1.000000D+00
        P   1   1.00
            9.618000D-02           1.000000D+00
        P   1   1.00
            3.590000D-02           1.0000000
        D   9   1.00
            1.326620D+02           5.720000D-04
            3.760540D+01           4.402000D-03
            1.038910D+01          -4.092200D-02
            6.490170D+00           9.966100D-02
            3.454510D+00           3.226630D-01
            1.844130D+00           4.003430D-01
            9.624780D-01           2.683060D-01
            4.728530D-01           8.484700D-02
            1.932000D-01           7.632000D-03
        D   1   1.00
            4.728530D-01           1.000000D+00
        D   1   1.00
            1.932000D-01           1.000000D+00
        D   1   1.00
            7.820000D-02           1.0000000
        F   1   1.00
            4.064000D-01           1.000000D+00
        F   1   1.00
            1.848000D-01           1.0000000
        ****
        I     0
        I-ECP     4     28
        g potential
        1
        2      1.00000000             0.00000000
        s-g potential
        3
        2     40.03337600            49.98964900
        2     17.30057600           281.00655600
        2      8.85172000            61.41673900
        p-g potential
        4
        2     15.72014100            67.41623900
        2     15.20822200           134.80769600
        2      8.29418600            14.56654800
        2      7.75394900            28.96842200
        d-g potential
        4
        2     13.81775100            35.53875600
        2     13.58780500            53.33975900
        2      6.94763000             9.71646600
        2      6.96009900            14.97750000
        f-g potential
        4
        2     18.52295000           -20.17661800
        2     18.25103500           -26.08807700
        2      7.55790100            -0.22043400
        2      7.59740400            -0.22164600
        """,
        key="basis",
        name="ECP+valence",
    )

    psi4.set_options(
        {
            "e_convergence": 1e-8,
            "d_convergence": 1e-8,
        }
    )

    psi4.core.set_output_file("output_ecp_plus_valence.dat", False)

    psi4.energy("scf")
    e1 = psi4.variable("HF TOTAL ENERGY")

    ### valence only
    psi4.basis_helper(
        """
        assign aug-cc-pvtz
        assign I noecp
        [noecp]
        I     0
        S   11   1.00
            5.546500D+03           1.560000D-04
            8.382140D+02           9.860000D-04
            1.821870D+02           2.792000D-03
            3.121230D+01          -4.325100D-02
            1.953140D+01           2.341340D-01
            8.240990D+00          -7.509430D-01
            2.194550D+00           8.829680D-01
            1.109110D+00           4.620610D-01
            3.746410D-01           2.228600D-02
            1.770800D-01          -4.353000D-03
            8.106100D-02           1.102000D-03
        S   11   1.00
            5.546500D+03          -7.300000D-05
            8.382140D+02          -5.080000D-04
            1.821870D+02          -1.158000D-03
            3.121230D+01           1.219300D-02
            1.953140D+01          -8.785400D-02
            8.240990D+00           3.382000D-01
            2.194550D+00          -5.765500D-01
            1.109110D+00          -4.092980D-01
            3.746410D-01           5.674590D-01
            1.770800D-01           6.124890D-01
            8.106100D-02           1.432310D-01
        S   1   1.00
            3.746410D-01           1.000000D+00
        S   1   1.00
            1.770800D-01           1.000000D+00
        S   1   1.00
            1.001000D-01           1.000000D+00
        S   1   1.00
            4.120000D-02           1.0000000
        P   9   1.00
            1.889880D+02           5.850000D-04
            2.128680D+01           3.692300D-02
            1.003960D+01          -2.353240D-01
            3.451800D+00           3.414830D-01
            1.974560D+00           5.347880D-01
            1.024200D+00           2.651410D-01
            4.494370D-01           2.578700D-02
            1.866480D-01           5.220000D-04
            7.348100D-02           6.060000D-04
        P   9   1.00
            1.889880D+02          -2.560000D-04
            2.128680D+01          -1.168200D-02
            1.003960D+01           8.319200D-02
            3.451800D+00          -1.569700D-01
            1.974560D+00          -2.245180D-01
            1.024200D+00          -1.144510D-01
            4.494370D-01           3.753560D-01
            1.866480D-01           5.751360D-01
            7.348100D-02           2.459170D-01
        P   1   1.00
            5.981000D-01           1.000000D+00
        P   1   1.00
            9.618000D-02           1.000000D+00
        P   1   1.00
            3.590000D-02           1.0000000
        D   9   1.00
            1.326620D+02           5.720000D-04
            3.760540D+01           4.402000D-03
            1.038910D+01          -4.092200D-02
            6.490170D+00           9.966100D-02
            3.454510D+00           3.226630D-01
            1.844130D+00           4.003430D-01
            9.624780D-01           2.683060D-01
            4.728530D-01           8.484700D-02
            1.932000D-01           7.632000D-03
        D   1   1.00
            4.728530D-01           1.000000D+00
        D   1   1.00
            1.932000D-01           1.000000D+00
        D   1   1.00
            7.820000D-02           1.0000000
        F   1   1.00
            4.064000D-01           1.000000D+00
        F   1   1.00
            1.848000D-01           1.0000000
      ****
        """,
        key="basis",
        name="only_valence",
    )

    psi4.set_options(
        {
            "e_convergence": 1e-8,
            "d_convergence": 1e-8,
        }
    )

    psi4.core.set_output_file("output_only_valence.dat", False)

    psi4.energy("scf")
    e2 = psi4.variable("HF TOTAL ENERGY")

    ### Compare
    diff = e2 - e1

    assert compare_values(0.0, diff, 7, "ECP + valence vs. only valence")



@uusing("ecpint")
@pytest.mark.medlong
def test_ecp_la2o3_grad_findif():
    """Analytic HF gradient of La2O3 (two 46-electron La ECPs) must match 5-point finite differences of energies.

    Before libecpint's fixes for narrow type-2 radial integrands and the early-stopping one-point quadrature,
    ECP integrals jumped discontinuously with geometry: the energy was off by 7.5e-5 Eh and findif deviated from
    analytic by 1.6e-4 Eh/a0 (no amount of tightening convergence helped). Now agreement is ~1e-7.
    """
    import psi4

    psi4.geometry("""
        La  0.0      0.944863  -0.3286479
        La  0.0     -0.944863  -0.3286479
        O   0.0      1.88973    1.5610781
        O   0.0      0.0        1.5610781
        O   0.0     -1.88973    1.5610781
        units au
        no_com
        no_reorient
    """)
    psi4.set_options({
        "scf_type": "pk",
        "basis": "def2-svp",
        "puream": True,
        "e_convergence": 10,
        "d_convergence": 10,
        "points": 5,
    })

    analytic, wfn = psi4.gradient("hf", dertype=1, return_wfn=True)
    findif = psi4.gradient("hf", dertype=0)
    assert compare_values(analytic, findif, 6, "La2O3 HF gradient: analytic vs findif")
    assert compare_values(-261.405809790, wfn.energy(), 6, "La2O3 HF energy")


@uusing("ecpint")
def test_ecp_ghost_beside_real_same_element():
    """A ghost atom carries no ECP, even when a real atom with the same label does.

    CH2I radical in the CH2I2 basis (one iodine ghosted). Previously the ECP shells were keyed by label,
    so Gh(I) picked up the real I's ECP: E(UHF) was 2.2e-3 Eh too high. References from NWChem 7.2.3
    (I def2-svp + def2-ecp from its library, bqI with basis only) at the same geometry;
    outputs in test_ecp/nwchem_ch2i_gi_pure_{hf,mp2ae}.out.
    """
    import psi4

    mol = psi4.geometry("""
        0 2
        C
        I     1 2.15
        Gh(I) 1 2.15 2 112.0
        X     1 1.0  2 56.0  3 0.0
        H     1 1.09 4 56.5  2 90.0
        H     1 1.09 4 56.5  2 -90.0
        symmetry c1
        no_com
        no_reorient
    """)
    psi4.set_options({
        "basis": "def2-svp",
        "reference": "uhf",
        "scf_type": "pk",
        "mp2_type": "conv",
        "e_convergence": 10,
        "d_convergence": 10,
    })

    nwchem_scf = -335.495596142434
    nwchem_scf_grad = np.array([
        [ 0.026200,  0.046304,  0.000000],
        [-0.030250,  0.081196,  0.000000],
        [-0.000487,  0.000317,  0.000000],
        [ 0.002268, -0.063909,  0.013732],
        [ 0.002268, -0.063909, -0.013732],
    ])
    nwchem_mp2 = -335.762780255744
    nwchem_mp2_grad = np.array([
        [ 0.011159,  0.058749,  0.000000],
        [-0.026022,  0.042160,  0.000000],
        [-0.002514,  0.001507,  0.000000],
        [ 0.008688, -0.051208,  0.020881],
        [ 0.008688, -0.051208, -0.020881],
    ])

    grad, wfn = psi4.gradient("hf", return_wfn=True)
    bs = wfn.basisset()
    assert bs.n_ecp_core() == 28
    assert [bs.n_ecp_shell_on_center(a) for a in range(mol.natom())] == [0, 4, 0, 0, 0]
    assert compare_values(nwchem_scf, wfn.energy(), 7, "CH2I + Gh(I) UHF energy")
    assert compare_values(nwchem_scf_grad, grad, 5, "CH2I + Gh(I) UHF gradient")

    grad = psi4.gradient("mp2")
    assert compare_values(nwchem_mp2, psi4.variable("CURRENT ENERGY"), 7, "CH2I + Gh(I) MP2 energy")
    assert compare_values(nwchem_mp2_grad, grad, 5, "CH2I + Gh(I) MP2 gradient")


@uusing("ecpint")
def test_ecp_ghost_cp_same_element():
    """Counterpoise through nbody ghosts by fragment; each monomer-in-dimer-basis job has a real I beside a ghost I.

    References from NWChem 7.2.3 (I def2-svp + def2-ecp, bq ghosts); outputs in test_ecp/nwchem_hi2_cp_*.out.
    CP IE = AB - A(AB) - B(AB). The CP gradient psi4 returns is the CP-corrected total (return_total_data=True default
    for gradients), which also adds the monomers in their own basis: AB - A(AB) - B(AB) + A(A) + B(B).
    """
    import psi4

    psi4.geometry("""
        0 1
        H  0.0  0.0  0.0
        I  0.0  0.0  1.609
        --
        0 1
        H  0.0  2.9  2.4
        I  0.0  4.2  3.3
        symmetry c1
        no_com
        no_reorient
    """)
    psi4.set_options({
        "basis": "def2-svp",
        "scf_type": "pk",
        "e_convergence": 10,
        "d_convergence": 10,
    })

    nwchem_cp_ie = 0.001110341095
    nwchem_cp_grad = np.array([  # CP-corrected total gradient; NWChem prints 6 decimals, five jobs summed
        [0.0, -0.000011, -0.002939],
        [0.0,  0.003179,  0.004034],
        [0.0,  0.006250,  0.004895],
        [0.0, -0.009416, -0.005990],
    ])

    psi4.energy("hf", bsse_type="cp")
    assert compare_values(nwchem_cp_ie, psi4.variable("CP-CORRECTED INTERACTION ENERGY"), 8, "HI dimer CP IE")
    grad = psi4.gradient("hf", bsse_type="cp")
    assert compare_values(nwchem_cp_grad, grad, atol=3.e-6, label="HI dimer CP gradient")
