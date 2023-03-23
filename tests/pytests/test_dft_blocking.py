import numpy as np
import pytest
import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def test_dft_atomic_blocking():
    """calculate XC energy per atom with psi4numpy"""

    mol = psi4.geometry(
        """
    0 1
    O 0.000000000000  0.000000000000 -0.068516219310
    H 0.000000000000 -0.790689573744  0.543701060724
    H 0.000000000000  0.790689573744  0.543701060724
    no_com
    no_reorient
    symmetry c1
    """
    )
    psi4.set_options(
        {
            "BASIS": "def2-SVP",
            "DFT_BLOCK_SCHEME": "ATOMIC",
            "DFT_SPHERICAL_POINTS": 590,
            "DFT_RADIAL_POINTS": 85,
            "DFT_REMOVE_DISTANT_POINTS": 0,
            "D_convergence": 1e-8,
        }
    )

    svwn_w, wfn = psi4.energy("SVWN", return_wfn=True)
    xc_ref = -8.9396369264084168  # Default Octree blocking
    xc_psi4 = wfn.variable("DFT XC ENERGY")
    Vpot = wfn.V_potential()
    D = np.array(wfn.Da())
    V = np.zeros_like(D)

    xc_e = 0.0

    rho = []
    points_func = Vpot.properties()[0]
    superfunc = Vpot.functional()
    grid = Vpot.grid()
    # Loop over atoms and their blocks
    for atom in range(mol.nallatom()):
        xc_atom = 0
        for block in grid.atomic_blocks()[atom]:

            # Obtain block information
            points_func.compute_points(block)
            npoints = block.npoints()
            lpos = np.array(block.functions_local_to_global())
            psi4.core.print_out(f"This block belongs to atom {block.parent_atom()}\n")

            # Obtain the grid weight
            w = np.array(block.w())

            # Compute phi!
            phi = np.array(points_func.basis_values()["PHI"])[:npoints, : lpos.shape[0]]

            # Build a local slice of D
            lD = D[(lpos[:, None], lpos)]

            # Copmute rho
            rho = 2.0 * np.einsum("pm,mn,pn->p", phi, lD, phi, optimize=True)

            inp = {}
            inp["RHO_A"] = psi4.core.Vector.from_array(rho)

            # Compute the kernel
            ret = superfunc.compute_functional(inp, -1, True)

            # Compute the XC energy
            vk = np.array(ret["V"])[:npoints]
            xc_atom += np.einsum("a,a->", w, vk, optimize=True)
        xc_e += xc_atom
        psi4.core.print_out(f"XC Energy per Atom {xc_atom}\n")

    psi4.core.print_out(f"Total XC Energy {xc_e}\n")

    assert psi4.compare_values(xc_ref, xc_psi4, 7, "psi4 XC energy")
    assert psi4.compare_values(xc_ref, xc_e, 7, "psi4numpy XC energy")


@pytest.mark.parametrize("scheme", ["OCTREE", "NAIVE", "ATOMIC"])
def test_dft_block_schemes(scheme):
    """all DFT_BLOCK_SCHEME should give same results and number
    of grid points. Water dimer with ghost atoms"""

    mol = psi4.geometry(
        """
    0 1
    O           -1.490196515110    -0.043256842172     0.000000000000
    H           -1.845932568294     0.844902886698     0.000000000000
    H           -0.533258283804     0.073267064698     0.000000000000
    @O            1.416663724802     0.038738966977     0.000000000000
    @H            1.773104797767    -0.423233996755     0.760023878024
    @H            1.773104797767    -0.423233996755    -0.760023878024
    no_com
    no_reorient
    symmetry c1
    """
    )
    psi4.set_options(
        {
            "BASIS": "def2-SVPD",
            "DFT_SPHERICAL_POINTS": 590,
            "DFT_RADIAL_POINTS": 85,
            "D_convergence": 1e-8,
            "DFT_WEIGHTS_TOLERANCE": -1.0,
        }
    )
    ref = {"XC GRID TOTAL POINTS": 300900, "DFT XC ENERGY": -9.218561399189895}
    psi4.set_options({"DFT_BLOCK_SCHEME": scheme})
    e, wfn = psi4.energy("pbe/def2-SVPD", return_wfn=True)
    P = psi4.variable("XC GRID TOTAL POINTS")
    XC = wfn.variable("DFT XC ENERGY")
    assert psi4.compare_integers(ref["XC GRID TOTAL POINTS"], P, f" {scheme} GRID POINTS:")
    assert psi4.compare_values(ref["DFT XC ENERGY"], XC, f" {scheme} XC ENERGY:")


def test_dft_block_scheme_distantpoints():
    """Test removal of distant grid points. all DFT_BLOCK_SCHEME should give same results and number
    of grid points."""

    mol = psi4.geometry(
        """
    0 1
    O  -1.551007  -0.114520   0.000000
    H  -1.934259   0.762503   0.000000
    H  -0.599677   0.040712   0.000000
    no_com
    no_reorient
    symmetry c1
    """
    )
    psi4.set_options(
        {
            "BASIS": "sto-3g",
            "maxiter": 1,
            "FAIL_ON_MAXITER": False,
            "DFT_PRUNING_SCHEME": "ROBUST",
            "DFT_WEIGHTS_TOLERANCE": -1.0,
        }
    )
    ref = {"True": 45929, "False": 46890}
    YESNO = [True,False]
    SCHEMES = ["OCTREE", "NAIVE", "ATOMIC"]
    for YN in YESNO:
        psi4.set_options({"DFT_REMOVE_DISTANT_POINTS":YN})
        for S in SCHEMES:
            psi4.set_options({"DFT_BLOCK_SCHEME": S})
            e, wfn = psi4.energy("pbe", return_wfn=True)
            P = psi4.variable("XC GRID TOTAL POINTS")
            XC = wfn.variable("DFT XC ENERGY")
            assert psi4.compare_integers(ref[f"{YN}"], P, f" scheme={S}; distant points={YN} ")
