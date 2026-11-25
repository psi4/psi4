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


# ---- ghost/ECP regression set -------------------------------------------------------------------------------------
# Which centers carry an ECP must follow each center's real/ghost status, never its element or label. A ghost keeps
# its orbital basis and contributes no ECP; a real ECP atom's potential acts on every basis function, ghost-centered
# ones included. Checks are internal (no external refs): ECP-matrix identities in a fixed basis (exact), then
# analytic derivatives vs findif, translational invariance summed over all centers (ghosts included), and a
# nonzero ghost-row gradient, so ghost-center (basis) derivatives of the real ECP integrals are exercised.
# HF/def2-svp, PK; I carries the 28-electron def2 ECP.

_ghost_ecp_systems = {
    # name: (charge, [(label, x, y, z) Angstrom], ECP shells per center)
    "F+GhI": (-1, [("Gh(I)", 0.0, 0.0, 0.0), ("F", 0.0, 0.0, 1.91)], [0, 0]),
    "I+GhF": (-1, [("I", 0.0, 0.0, 0.0), ("Gh(F)", 0.0, 0.0, 1.91)], [4, 0]),
    "I+GhI": (-1, [("I", 0.0, 0.0, 0.0), ("Gh(I)", 0.0, 0.0, 2.67)], [4, 0]),
    "GhI+I": (-1, [("Gh(I)", 0.0, 0.0, 0.0), ("I", 0.0, 0.0, 2.67)], [0, 4]),
    "I+I": (0, [("I", 0.0, 0.0, 0.0), ("I", 0.0, 0.0, 2.67)], [4, 4]),
    "IF": (0, [("I", 0.0, 0.0, 0.0), ("F", 0.0, 0.0, 1.91)], [4, 0]),
    # bent, so ECP derivative integrals with three distinct centers (ghost bra, real ket, ECP on a third atom) enter
    "IF+GhI": (0, [("I", 0.0, 0.0, 0.0), ("F", 0.0, 0.0, 1.91), ("Gh(I)", 2.630, 0.0, -0.464)], [4, 0, 0]),
}


def _ghost_ecp_molecule(name, fixed_frame=True):
    import psi4

    charge, atoms, _ = _ghost_ecp_systems[name]
    lines = "\n".join(f"{lbl} {x} {y} {z}" for lbl, x, y, z in atoms)
    frame = "no_com\nno_reorient\n" if fixed_frame else ""
    return psi4.geometry(f"{charge} 1\n{lines}\nunits angstrom\nsymmetry c1\n{frame}", name=name)


def _ecp_nuclear_repulsion(mol, name):
    """NRE from ECP-reduced charges (I: 53 - 28), ghosts zero, independent of psi4's Molecule.Z bookkeeping."""
    zeff = {"I": 25.0, "F": 9.0}
    z = [0.0 if lbl.startswith("Gh(") else zeff[lbl] for lbl, *_ in _ghost_ecp_systems[name][1]]
    geom = np.asarray(mol.geometry())
    return sum(z[i] * z[j] / np.linalg.norm(geom[i] - geom[j])
               for i in range(len(z)) for j in range(i + 1, len(z)))


@uusing("ecpint")
def test_ecp_ghost_integral_identities():
    """ECP matrices in fixed bases: ghost I gives none; a ghost neighbor doesn't change the real I's ECP; two real
    iodines' ECP is the sum of each alone with the other ghosted (fails if ECP centers are picked by element/label)."""
    import psi4

    V, nbf = {}, {}
    for name, (_, atoms, shells) in _ghost_ecp_systems.items():
        mol = _ghost_ecp_molecule(name)
        bs = psi4.core.BasisSet.build(mol, "ORBITAL", "def2-svp")
        assert [bs.n_ecp_shell_on_center(a) for a in range(mol.natom())] == shells, name
        assert bs.n_ecp_core() == 28 * sum(n > 0 for n in shells), name
        V[name], nbf[name] = psi4.core.MintsHelper(bs).ao_ecp().np, bs.nbf()

    # ghosts keep their orbital basis
    assert nbf["F+GhI"] == nbf["I+GhF"] == nbf["IF"]
    assert nbf["I+GhI"] == nbf["GhI+I"] == nbf["I+I"]

    assert compare_values(np.zeros_like(V["F+GhI"]), V["F+GhI"], atol=1.e-14, label="F + Gh(I): no ECP")
    assert compare_values(V["IF"], V["I+GhF"], atol=1.e-12, label="I + Gh(F) == I + F")
    assert compare_values(V["I+I"], V["I+GhI"] + V["GhI+I"], atol=1.e-12, label="I + I == (I + Gh(I)) + (Gh(I) + I)")


@uusing("ecpint")
@pytest.mark.d2ints
@pytest.mark.parametrize("name", [
    pytest.param("F+GhI"),
    pytest.param("I+GhF"),
    pytest.param("I+GhI"),
    pytest.param("I+I", marks=pytest.mark.medlong),
    pytest.param("I+I-free"),
    pytest.param("IF+GhI", marks=pytest.mark.medlong),
])
def test_ecp_ghost_derivatives(name):
    """Analytic HF gradient vs 5-point findif of energies (Cartesian steps) and analytic Hessian vs findif of gradients,
    over all centers including ghosts. Prototype agreement: gradients <= 1e-10, Hessians <= 3e-8.

    Also, analytic and findif paths must agree on the side results: NUCLEAR REPULSION ENERGY from ECP-reduced charges
    (findif's own molecule never builds the basis, so it once reported full-Z NRE), and the caller's fix_com/
    fix_orientation left as given (findif once froze them on the input molecule). "-free" runs in psi4's frame."""
    import psi4

    fixed_frame = not name.endswith("-free")
    name = name.removesuffix("-free")
    mol = _ghost_ecp_molecule(name, fixed_frame=fixed_frame)
    mol.update_geometry()
    natom = mol.natom()
    nre = _ecp_nuclear_repulsion(mol, name)

    def assert_side_results(label, wfn):
        for store, val in [("core", psi4.variable("NUCLEAR REPULSION ENERGY")),
                           ("wfn", wfn.variable("NUCLEAR REPULSION ENERGY"))]:
            assert compare_values(nre, val, atol=1.e-10, label=f"{name} {label} {store} NUCLEAR REPULSION ENERGY")
        assert (mol.com_fixed(), mol.orientation_fixed()) == (fixed_frame, fixed_frame), \
            f"{name} {label} changed input molecule frame flags to {(mol.com_fixed(), mol.orientation_fixed())}"

    ghosts = [a for a, (lbl, *_) in enumerate(_ghost_ecp_systems[name][1]) if lbl.startswith("Gh(")]
    psi4.set_options({
        "basis": "def2-svp",
        "scf_type": "pk",
        "e_convergence": 11,
        "d_convergence": 11,
        "points": 5,
        "fd_project": False,
        "disp_mass_weight": False,
    })

    # findif first: building a basis with ECPs lowers Z on the molecule it's built on, so an analytic calc beforehand
    #   would mask a findif NRE taken from its own (never basis-built) molecule
    g0, wfn = psi4.gradient("hf", molecule=mol, dertype=0, return_wfn=True)
    assert_side_results("grd0", wfn)
    g1, wfn = psi4.gradient("hf", molecule=mol, return_wfn=True)
    assert_side_results("grd1", wfn)
    g1, g0 = g1.np, g0.np
    assert compare_values(g0, g1, atol=1.e-8, label=f"{name} analytic vs findif gradient")
    assert compare_values(np.zeros(3), g1.sum(axis=0), atol=1.e-9, label=f"{name} gradient translational sum")
    for a in ghosts:
        assert np.abs(g1[a]).max() > 1.e-5, f"{name} ghost center {a} gradient vanishes: {g1[a]}"

    h2, wfn = psi4.hessian("hf", molecule=mol, return_wfn=True)
    assert_side_results("hes2", wfn)
    h1, wfn = psi4.hessian("hf", molecule=mol, dertype=1, return_wfn=True)
    assert_side_results("hes1", wfn)
    h2, h1 = h2.np, h1.np
    assert compare_values(h1, h2, atol=1.e-7, label=f"{name} analytic vs findif Hessian")
    assert compare_values(np.zeros((natom, 3, 3)), h2.reshape(natom, 3, natom, 3).sum(axis=2), atol=1.e-9,
                          label=f"{name} Hessian translational sum")

    if not fixed_frame:
        # QCSchema results: properties.nuclear_repulsion_energy came from the top-level (full-Z) molecule for findif
        for dertype in [1, 0]:
            atin = {
                "molecule": mol.to_schema(dtype=3),
                "specification": {
                    "driver": "gradient",
                    "model": {"method": "hf", "basis": "def2-svp"},
                    "keywords": {"scf_type": "pk", "e_convergence": 11, "d_convergence": 11, "points": 5,
                                 "fd_project": False, "disp_mass_weight": False, "function_kwargs": {"dertype": dertype}}},
            }
            atres = psi4.schema_wrapper.run_qcschema(atin)
            assert atres.success, atres.error.error_message
            for store, val in [("properties", atres.properties.nuclear_repulsion_energy),
                               ("qcvars", atres.extras["qcvars"]["NUCLEAR REPULSION ENERGY"])]:
                assert compare_values(nre, float(val), atol=1.e-10,
                                      label=f"{name} schema grd{dertype} {store} NUCLEAR REPULSION ENERGY")
            assert compare_values(g1, np.asarray(atres.return_result).reshape(-1, 3), atol=1.e-8,
                                  label=f"{name} schema grd{dertype} vs analytic gradient")


@uusing("ecpint")
@pytest.mark.d2ints
def test_ecp_vibanal_isotopologue():
    """vibanal_wfn with a separately defined molecule (isotope analysis) must accept a matching molecule though building
    an ECP basis lowered Z on wfn's molecule (I: 53 -> 25), and still reject a different geometry. Its NRE check once
    compared wfn's ECP-reduced NRE (8.217) with the fresh molecule's full-Z one (17.420) and raised for any ECP system."""
    import psi4

    hi = psi4.geometry("H 0 0 0\nI 0 0 1.61\nsymmetry c1\nno_com\nno_reorient\n", name="hi")
    psi4.set_options({"basis": "def2-svp", "scf_type": "pk", "d_convergence": 10})
    _, wfn = psi4.frequency("hf", molecule=hi, return_wfn=True)
    omega_hi = wfn.frequencies().np[-1]

    di = psi4.geometry("H@2.014101778 0 0 0\nI 0 0 1.61\nsymmetry c1\nno_com\nno_reorient\n", name="di")
    omega_di = psi4.vibanal_wfn(wfn, molecule=di)["omega"].data[-1].real
    m_H, m_D, m_I = 1.00782503223, 2.014101778, hi.mass(1)
    ratio = np.sqrt((m_H * m_I / (m_H + m_I)) / (m_D * m_I / (m_D + m_I)))
    assert compare_values(omega_hi * ratio, omega_di, atol=1.e-3, label="DI harmonic frequency from HI Hessian")

    stretched = psi4.geometry("H 0 0 0\nI 0 0 1.71\nsymmetry c1\nno_com\nno_reorient\n", name="hi_stretched")
    with pytest.raises(psi4.ValidationError, match="Impostor molecule"):
        psi4.vibanal_wfn(wfn, molecule=stretched)


@uusing("ecpint")
@pytest.mark.parametrize("method,kw", [
    pytest.param("sapt0", {}, id="sapt0"),
    pytest.param("sapt(dft)", {"sapt_dft_functional": "hf", "sapt_dft_grac_shift_a": 0.0,
                               "sapt_dft_grac_shift_b": 0.0, "sapt_dft_do_dhf": False}, id="sapt(dft)"),
    pytest.param("sf-sapt", {"reference": "rohf"}, id="sf-sapt"),
])
def test_ecp_sapt_qcschema_nre(method, kw):
    """SAPT through QCSchema must report the dimer's ECP-reduced NRE (I: 53 -> 25) in properties. The top-level molecule
    is split into dimer/monomers before any basis is built, so it once gave full Z (here 379.58 vs 96.79). Also checks
    SF-SAPT's returned dimer wfn, once built on monomer A's basis (B ghosted), so its molecule kept full Z."""
    import psi4

    dimer = psi4.geometry("""
0 1
H 0.0 0.0 -3.0
I 0.0 0.0 -1.39
--
0 1
H 0.0 0.0 1.5
I 0.0 0.0 3.11
units angstrom
symmetry c1
no_com
no_reorient
""", name="hi2")
    zeff = [1.0, 25.0, 1.0, 25.0]
    geom = np.asarray(dimer.geometry())
    nre = sum(zeff[i] * zeff[j] / np.linalg.norm(geom[i] - geom[j]) for i in range(4) for j in range(i + 1, 4))

    atin = {
        "molecule": dimer.to_schema(dtype=3),
        "specification": {
            "driver": "energy",
            "model": {"method": method, "basis": "def2-svp"},
            "keywords": {"scf_type": "df", "e_convergence": 8, "d_convergence": 8, **kw}},
    }
    atres = psi4.schema_wrapper.run_qcschema(atin)
    assert atres.success, atres.error.error_message
    assert compare_values(nre, atres.properties.nuclear_repulsion_energy, atol=1.e-8,
                          label=f"{method} schema properties NUCLEAR REPULSION ENERGY")

    if method == "sf-sapt":
        psi4.set_options({"basis": "def2-svp", "scf_type": "df", **kw})
        _, wfn = psi4.energy(method, molecule=dimer, return_wfn=True)
        assert compare_values(nre, wfn.molecule().nuclear_repulsion_energy(), atol=1.e-8,
                              label="sf-sapt dimer wfn molecule NRE")
        assert wfn.basisset().n_ecp_core() == 56, "sf-sapt dimer wfn basis has both iodines' ECPs"
