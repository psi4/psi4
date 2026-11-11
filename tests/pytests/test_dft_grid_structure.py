"""Structural checks on single-atom DFT grids.

The invariant that matters is *which angular grid sits on which radial shell*, judged by the shell's
radius. Radial schemes are free to emit their points in either order (TREUTLER and MULTIEXP run
outermost-first; BECKE, EM, and MURA run nucleus-first), so anything that assigns angular grids by
array index rather than by radius is a bug. These tests observe only the finished grid.
"""

import math

import numpy as np
import pytest
import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]

RADIAL_SCHEMES = ["TREUTLER", "BECKE", "MULTIEXP", "EM", "MURA"]
REGION_PRUNING = ["ROBUST", "TREUTLER"]
FUNCTION_PRUNING = ["FLAT", "P_SLATER", "D_SLATER", "LOG_SLATER", "P_GAUSSIAN", "D_GAUSSIAN", "LOG_GAUSSIAN"]

NRAD = 75
NSPH = 302
ELEMENTS = {"H": (1, 2), "C": (6, 3), "Cl": (17, 2)}  # symbol: (Z, multiplicity)

# Lebedev (order, npoints) pairs, as in LebedevGridMgr (standard published values)
# yapf: disable
LEBEDEV = [(0, 1), (3, 6), (5, 14), (7, 26), (9, 38), (11, 50), (13, 74), (15, 86), (17, 110), (19, 146),
           (21, 170), (23, 194), (25, 230), (27, 266), (29, 302), (31, 350), (35, 434), (41, 590), (47, 770),
           (53, 974), (59, 1202), (65, 1454), (71, 1730), (77, 2030), (83, 2354), (89, 2702), (95, 3074),
           (101, 3470), (107, 3890), (113, 4334), (119, 4802), (125, 5294), (131, 5810)]
# yapf: enable
# Bragg-Slater radii [bohr] for the elements tested, copied from psi4's GetBSRadius (NOT independently verified)
BS_RADIUS = {1: 0.661, 6: 1.323, 17: 1.890}


def _order_of(npts):
    return dict((n, o) for o, n in LEBEDEV)[npts]


def _npts_round_up(order):
    return next(n for o, n in LEBEDEV if o >= order)


def _rle(seq):
    """Expand [(value, repeat), ...] into a flat list."""
    return [v for v, n in seq for _ in range(n)]


# ==> Named grids: angular point counts from the nucleus outward <==

# symbol: (Z, multiplicity) for elements with SG-0/SG-1 definitions in psi4
NAMED_ELEMENTS = {
    "H": (1, 2), "Li": (3, 2), "Be": (4, 1), "B": (5, 2), "C": (6, 3), "N": (7, 4), "O": (8, 3), "F": (9, 2),
    "Na": (11, 2), "Mg": (12, 1), "Al": (13, 2), "Si": (14, 3), "P": (15, 4), "S": (16, 3), "Cl": (17, 2),
}  # yapf: disable

# Chien & Gill, J. Comput. Chem. 27 (2006) 730, Table 1, transcribed from the paper (not from psi4).
# symbol: (SG-0 Ntot, SG-1 Ntot, SG-0 Lebedev partition [(npoints, nshells), ...] from the nucleus outward)
# The table does not state the partition direction; nucleus-outward is established by SG-0 accuracy
# (see test_dft_named_grids.py) and by the 6-point rule belonging at the nucleus as in SG-1.
# yapf: disable
CHIEN_GILL_TABLE1 = {
    "H": (1406, 3752, [(6, 6), (18, 3), (26, 1), (38, 1), (74, 1), (110, 1), (146, 6), (86, 1), (50, 1), (38, 1), (18, 1)]),
    "Li": (1406, 3816, [(6, 6), (18, 3), (26, 1), (38, 1), (74, 1), (110, 1), (146, 6), (86, 1), (50, 1), (38, 1), (18, 1)]),
    "Be": (1390, 3816, [(6, 4), (18, 2), (26, 1), (38, 2), (74, 1), (86, 1), (110, 2), (146, 5), (50, 1), (38, 1), (18, 1), (6, 2)]),
    "B": (1426, 3816, [(6, 4), (26, 4), (38, 3), (86, 3), (146, 6), (38, 1), (6, 2)]),
    "C": (1390, 3816, [(6, 6), (18, 2), (26, 1), (38, 2), (50, 2), (86, 1), (110, 1), (146, 1), (170, 2), (146, 2), (86, 1), (38, 1), (18, 1)]),
    "N": (1414, 3816, [(6, 6), (18, 3), (26, 1), (38, 2), (74, 2), (110, 1), (170, 2), (146, 3), (86, 1), (50, 2)]),
    "O": (1154, 3816, [(6, 5), (18, 1), (26, 2), (38, 1), (50, 4), (86, 1), (110, 5), (86, 1), (50, 1), (38, 1), (6, 1)]),
    "F": (1494, 3816, [(6, 4), (38, 2), (50, 4), (74, 2), (110, 2), (146, 2), (110, 2), (86, 3), (50, 1), (6, 1)]),
    "Na": (1328, 3760, [(6, 6), (18, 2), (26, 3), (38, 1), (50, 2), (110, 8), (74, 2), (6, 2)]),
    # The paper lists Ntot = 1492 for Mg, but its own partition sums to 1468; psi4 follows the partition.
    "Mg": (1468, 3760, [(6, 5), (18, 2), (26, 2), (38, 2), (50, 2), (74, 1), (110, 2), (146, 4), (110, 1), (86, 1), (38, 2), (18, 1), (6, 1)]),
    "Al": (1496, 3760, [(6, 6), (18, 2), (26, 1), (38, 2), (50, 2), (74, 1), (86, 1), (146, 2), (170, 2), (110, 2), (86, 1), (74, 1), (26, 1), (18, 1), (6, 1)]),
    "Si": (1496, 3760, [(6, 5), (18, 4), (38, 4), (50, 3), (74, 1), (110, 2), (146, 1), (170, 3), (86, 1), (50, 1), (6, 1)]),
    "P": (1496, 3760, [(6, 5), (18, 4), (38, 4), (50, 3), (74, 1), (110, 2), (146, 1), (170, 3), (86, 1), (50, 1), (6, 1)]),
    "S": (1456, 3760, [(6, 4), (18, 1), (26, 8), (38, 2), (50, 1), (74, 2), (110, 1), (170, 3), (146, 1), (110, 1), (50, 1), (6, 1)]),
    # Cl: NOT independently verified -- taken from psi4's Initialize_SG0; SG-1 Ntot follows the Na-Ar row.
    "Cl": (1480, 3760, [(6, 4), (18, 7), (26, 2), (38, 2), (50, 1), (74, 1), (110, 2), (170, 3), (146, 1), (110, 1), (86, 1), (6, 1)]),
}
# yapf: enable
SG0 = {el: _rle(part) for el, (_, _, part) in CHIEN_GILL_TABLE1.items()}
SG0_NPTS = {el: n0 for el, (n0, _, _) in CHIEN_GILL_TABLE1.items()}

# SG-1: Gill, Johnson & Pople, Chem. Phys. Lett. 209 (1993) 506. Totals are from Chien & Gill Table 1;
# the per-row partitions are psi4's Q-Chem-compatible choice and are NOT independently verified.
_SG1_ROW = {
    1: _rle([(6, 16), (38, 5), (86, 4), (194, 9), (86, 16)]),
    2: _rle([(6, 14), (38, 7), (86, 3), (194, 9), (86, 17)]),
    3: _rle([(6, 12), (38, 7), (86, 5), (194, 7), (86, 19)]),
}
SG1 = {el: _SG1_ROW[1 if Z <= 2 else 2 if Z <= 10 else 3] for el, (Z, _) in NAMED_ELEMENTS.items()}
SG1_NPTS = {el: n1 for el, (_, n1, _) in CHIEN_GILL_TABLE1.items()}


# ==> Expected angular counts for psi4's own pruning schemes <==


def _region_pruned_npts(rank_from_outside, Z, nrad, nominal_npts, scheme):
    """Mirror of RadialPruneMgr::{Shell,TreutlerShell}Pruning, written against the shell's rank
    counted from the *outermost* shell (rank 0), which is what those routines document they assume."""
    ri = rank_from_outside
    nominal = _order_of(nominal_npts)
    if Z >= 36:
        return _npts_round_up(nominal)
    order = nominal - 1 if Z <= 2 else nominal
    if scheme == "TREUTLER":
        if ri >= 2 * nrad / 3.0:
            order = 7
        elif ri >= nrad / 2.0:
            order = 11
    elif scheme == "ROBUST":
        nr = nrad / 4.0
        if ri <= nr:
            order = nominal
        elif ri <= 2 * nr:
            order = nominal
        elif ri <= 3 * nr:
            order = nominal - 6
        else:
            order = 7
    return _npts_round_up(order)


_PRUNE_FN = {
    "FLAT": lambda rho, a: 1.0,
    "P_SLATER": lambda rho, a: rho * math.exp(1 - rho),
    "D_SLATER": lambda rho, a: (rho * math.exp(1 - rho))**2,
    "LOG_SLATER": lambda rho, a: math.exp(-a * abs(math.log(rho))),
    "P_GAUSSIAN": lambda rho, a: rho * math.exp((1 - rho * rho) / 2.0),
    "D_GAUSSIAN": lambda rho, a: rho * rho * math.exp(1 - rho * rho),
    "LOG_GAUSSIAN": lambda rho, a: math.exp(-a * math.log(rho)**2),
}


def _function_pruned_npts(r, Z, nominal_npts, fn, bs_alpha=1.0, pruning_alpha=1.0):
    """Mirror of RadialPruneMgr::GetPrunedNumAngPts."""
    rho = r / (BS_RADIUS[Z] * bs_alpha)
    order = math.ceil(_order_of(nominal_npts) * _PRUNE_FN[fn](rho, pruning_alpha) - 1.0e-10)
    return _npts_round_up(order)


# ==> Grid construction and analysis <==


def _atom_shells(element, str_opts, int_opts=None):
    """Build a one-atom grid with no points discarded; return shell radii (ascending) and, for each
    shell, its angular point count and summed weight. Also returns the raw arrays."""
    Z, mult = NAMED_ELEMENTS[element]
    # Zero basis tolerance gives infinite basis-function extents, so no point is removed as distant
    # and no block is discarded for lacking basis functions; zero weights tolerance keeps every point.
    psi4.set_options({"dft_basis_tolerance": 0.0, "dft_weights_tolerance": 0.0})
    mol = psi4.geometry(f"0 {mult}\n{element} 0 0 0\nsymmetry c1\nno_com\nno_reorient")
    basis = psi4.core.BasisSet.build(mol, "ORBITAL", "STO-3G", quiet=True)
    grid = psi4.core.DFTGrid.build(mol, basis, int_opts or {}, str_opts)

    x, y, z, w = (np.concatenate([getattr(b, c)().np for b in grid.blocks()]) for c in "xyzw")
    r = np.sqrt(x * x + y * y + z * z)

    # Group points into shells by radius, with a relative tolerance (outer EM/BECKE shells reach ~1e4 bohr)
    order = np.argsort(r)
    rs, ws = r[order], w[order]
    breaks = np.nonzero(np.diff(rs) > 1.0e-8 * rs[1:])[0] + 1
    shells_r = np.split(rs, breaks)
    shells_w = np.split(ws, breaks)
    return {
        "Z": Z,
        "npoints": grid.npoints(),
        "radii": np.array([s.mean() for s in shells_r]),
        "angular": [len(s) for s in shells_r],
        "shell_w": np.array([s.sum() for s in shells_w]),
        "r": r,
        "w": w,
    }


def _check_radial_sanity(g, nshell):
    assert len(g["radii"]) == nshell, "wrong number of radial shells"
    assert g["npoints"] == len(g["r"]) == sum(g["angular"]), "point totals inconsistent"
    assert np.all(np.isfinite(g["radii"])) and g["radii"][0] > 0.0
    assert np.all(np.diff(g["radii"]) > 0.0)
    # Individual points may carry negative weight (e.g., the 74-point Lebedev rule), but each shell's
    # total (its radial weight times 4 pi) must be positive.
    assert np.all(g["shell_w"] > 0.0), "non-positive radial weight"
    # Integrate the normalized hydrogen 1s density; catches scrambled r <-> weight pairing
    norm = np.sum(g["w"] * np.exp(-2.0 * g["r"]) / np.pi)
    assert norm == pytest.approx(1.0, abs=1.0e-3)


# ==> Tests <==


@pytest.mark.parametrize("element", list(NAMED_ELEMENTS))
@pytest.mark.parametrize("name,table,npts", [("SG0", SG0, SG0_NPTS), ("SG1", SG1, SG1_NPTS)])
def test_named_grid(name, table, npts, element):
    g = _atom_shells(element, {"DFT_GRID_NAME": name})
    _check_radial_sanity(g, len(table[element]))
    assert g["npoints"] == npts[element]
    assert g["angular"] == table[element], f"{name} angular grids not applied nucleus-outward"


def test_sg0_energy_accuracy():
    """SG-0 should reproduce a fine-grid B3LYP energy to well under 1 mEh for water."""
    mol = psi4.geometry(
        """
    O 0.000  0.000  0.117
    H 0.000  0.757 -0.467
    H 0.000 -0.757 -0.467
    symmetry c1
    """
    )
    psi4.set_options(
        {
            "basis": "cc-pvdz",
            "scf_type": "direct",
            "e_convergence": 1e-10,
            "d_convergence": 1e-8,
            "dft_grid_name": "SG0",
        }
    )
    e_sg0 = psi4.energy("b3lyp", molecule=mol)

    # Fine-grid reference: psi4 (150, 974) Treutler grid, DIRECT (same as PK to 1e-10); (200, 1454) gives
    # -76.4202850183. Cross-checked against NWChem: -76.420285024 at both (150, 3074) and (200, 5810).
    e_ref = -76.4202850142
    # SG-0 error now about -0.32 mEh. With the 2013-2026 backwards SG-0 grid it was about -14.68 mEh.
    assert abs(e_sg0 - e_ref) < 1.0e-3
    assert psi4.compare_values(-76.4206068189, e_sg0, 6, "B3LYP/cc-pVDZ water, SG-0 grid")
    # 2013-2026 value with the backwards SG-0 grid (psi4archive 1dcefd12 until this fix); it used to pass:
    # assert psi4.compare_values(-76.4349631854, e_sg0, 6, "B3LYP/cc-pVDZ water, SG-0 grid (pre-fix)")


@pytest.mark.parametrize("element", list(ELEMENTS))
@pytest.mark.parametrize("pruning", ["NONE"] + REGION_PRUNING + FUNCTION_PRUNING)
@pytest.mark.parametrize("radial", RADIAL_SCHEMES)
def test_radial_pruning_grid(radial, pruning, element):
    g = _atom_shells(
        element,
        {"DFT_GRID_NAME": "", "DFT_RADIAL_SCHEME": radial, "DFT_PRUNING_SCHEME": pruning},
        {"DFT_RADIAL_POINTS": NRAD, "DFT_SPHERICAL_POINTS": NSPH},
    )
    Z = g["Z"]
    _check_radial_sanity(g, NRAD)

    if pruning == "NONE":
        expected = [NSPH] * NRAD
    elif pruning in REGION_PRUNING:
        # ascending-radius shell k is rank NRAD-1-k from the outside
        expected = [_region_pruned_npts(NRAD - 1 - k, Z, NRAD, NSPH, pruning) for k in range(NRAD)]
        # sanity of the expectation itself: never coarser going outward
        assert expected == sorted(expected)
    else:
        expected = [_function_pruned_npts(r, Z, NSPH, pruning) for r in g["radii"]]

    assert g["angular"] == expected, f"angular grids misassigned for {radial} radial + {pruning} pruning"
