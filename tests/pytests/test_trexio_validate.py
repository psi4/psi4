"""Validate Psi4's TrexIO files against the TrexIO specification.

trexio-validate (https://github.com/TREX-CoE/trexio-validate) recomputes the
contents of a TrexIO file from the basis set stored in it, with libcint and the
conventions of the specification rather than those of any one program: AO
order and solid-harmonic phases, normalization, MO orthonormality, the one- and
two-electron integrals in both bases. Psi4's own round-trip tests cannot catch a
convention error, because they undo it on the way back in; this can.

Two things trexio-validate cannot see are checked here as well: that each MO
carries its own occupation (it checks only their sum), and -- by rebuilding the
SCF energy from the file alone -- that the stored core Hamiltonian is the one
the calculation used. The second is what covers effective core potentials,
whose integrals the validator does not recompute.

Each case isolates one thing a writer can get wrong. The validator is found as
``$TREXIO_VALIDATE`` or ``trexio-validate`` on PATH. Without it these tests
skip -- unless ``PSI4_TREXIO_VALIDATE_REQUIRED`` is set, as CI does, so that a
missing validator can never turn the job green by checking nothing.
"""

import os
import shutil
import subprocess

import numpy as np
import pytest

import psi4

trexio = pytest.importorskip("trexio")

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.trexio]

_REQUIRED = bool(os.environ.get("PSI4_TREXIO_VALIDATE_REQUIRED"))

# Psi4 writes no dipole integrals; everything else is required to be present.
_NO_DIPOLES = [f"{b}_1e_int_dipole_{c}" for b in ("ao", "mo") for c in "xyz"]
_ERI_CHECKS = ["ao_2e_int_eri", "mo_2e_int_eri"]

_WATER = "O\nH 1 0.96\nH 1 0.96 2 104.5"

# (id, geometry, options, method, write ERIs, what it isolates)
CASES = [
    ("h2_sto3g", "0 1\nH\nH 1 0.74\nsymmetry c1", {"basis": "sto-3g"}, "hf", True,
     "s functions only"),
    ("h2o_ccpvdz_c1", "0 1\n" + _WATER + "\nsymmetry c1", {"basis": "cc-pvdz"}, "hf", True,
     "baseline s, p, d"),
    ("h2o_ccpvdz_c2v", "0 1\n" + _WATER, {"basis": "cc-pvdz"}, "hf", True,
     "point-group symmetry: SO -> AO back-transformation, MOs across irreps"),
    ("n2_ccpvtz_d2h", "0 1\nN\nN 1 1.10", {"basis": "cc-pvtz"}, "hf", True,
     "f functions; D2h subgroup of a linear molecule"),
    ("ne_ccpvqz", "0 1\nNe", {"basis": "cc-pvqz"}, "hf", True,
     "g functions: solid-harmonic phase at l=4"),
    ("ne_ccpv5z", "0 1\nNe", {"basis": "cc-pv5z"}, "hf", False,
     "h functions: solid-harmonic phase at l=5 (91 AOs, too many for the ERI checks)"),
    ("h2o_631gs_cart", "0 1\n" + _WATER + "\nsymmetry c1", {"basis": "6-31g*"}, "hf", True,
     "Cartesian d, as the basis file specifies"),
    ("hf_ccpvtz_cart", "0 1\nH\nF 1 0.92\nsymmetry c1", {"basis": "cc-pvtz", "puream": False}, "hf", True,
     "Cartesian f: per-AO normalization at l=3"),
    ("h2o_631g_sp", "0 1\n" + _WATER + "\nsymmetry c1", {"basis": "6-31g"}, "hf", True,
     "SP shells, split into separate s and p shells"),
    ("ch2_triplet_uhf", "0 3\nC\nH 1 1.08\nH 1 1.08 2 134.0\nsymmetry c1", {"basis": "cc-pvdz", "reference": "uhf"}, "hf", True,
     "UHF: alpha and beta MOs, mo.spin, all spin blocks of the MO integrals"),
    ("nh2_doublet_rohf", "0 2\nN\nH 1 1.02\nH 1 1.02 2 103.0\nsymmetry c1", {"basis": "cc-pvdz", "reference": "rohf"}, "hf", True,
     "ROHF: one MO set with doubly and singly occupied orbitals"),
    ("oh_uks_b3lyp", "0 2\nO\nH 1 0.97\nsymmetry c1", {"basis": "cc-pvdz", "reference": "uks"}, "b3lyp", True,
     "Kohn-Sham orbitals, unrestricted"),
    ("oh_anion_augdz", "-1 1\nO\nH 1 0.97\nsymmetry c1", {"basis": "aug-cc-pvdz"}, "hf", True,
     "anion with diffuse functions"),
    ("h3o_cation", "1 1\nO\nH 1 0.98\nH 1 0.98 2 110.0\nH 1 0.98 2 110.0 3 120.0\nsymmetry c1", {"basis": "cc-pvdz"}, "hf", True,
     "cation"),
    ("h2o_ghost_he", "0 1\n" + _WATER + "\n@He 1 3.0 2 90.0 3 180.0\nsymmetry c1", {"basis": "cc-pvdz"}, "hf", True,
     "ghost atom: basis functions on a centre with zero charge"),
    ("h2o_lindep", "0 1\n" + _WATER + "\nsymmetry c1", {"basis": "aug-cc-pvdz", "s_tolerance": 1.0e-2}, "hf", True,
     "linear dependencies removed, so mo.num < ao.num"),
    ("h2o_mixed_basis", "0 1\n" + _WATER + "\nsymmetry c1", {"basis": "mixed_h2o"}, "hf", True,
     "a different basis set per element"),
    ("hi_def2svp_ecp", "0 1\nH\nI 1 1.61\nsymmetry c1", {"basis": "def2-svp"}, "hf", True,
     "ECP on one atom but not the other: ecp group, Z_eff charges, V_ECP in the core Hamiltonian"),
    ("xe_def2svp_ecp", "0 1\nXe", {"basis": "def2-svp"}, "hf", True,
     "single ECP atom"),
    ("hi_def2svp_ecp_cart", "0 1\nH\nI 1 1.61\nsymmetry c1", {"basis": "def2-svp", "puream": False}, "hf", True,
     "ECP integrals in a Cartesian basis: the per-AO normalization applies to them too"),
    ("ag_doublet_uhf_ecp", "0 2\nAg\nsymmetry c1", {"basis": "def2-svp", "reference": "uhf"}, "hf", True,
     "ECP with an unrestricted wavefunction"),
    ("h2o_ghost_i_ecp", "0 1\n" + _WATER + "\n@I 1 3.5 2 90.0 3 180.0\nsymmetry c1", {"basis": "def2-svp"}, "hf", True,
     "ghost atom with an ECP basis: its basis functions, but no core removed and no ECP applied"),
]


def _validator():
    exe = os.environ.get("TREXIO_VALIDATE") or shutil.which("trexio-validate")
    if exe is None:
        msg = "trexio-validate not found (set TREXIO_VALIDATE or put it on PATH)"
        if _REQUIRED:
            pytest.fail(msg + ", and PSI4_TREXIO_VALIDATE_REQUIRED is set")
        pytest.skip(msg)
    return exe


def _require_hdf5():
    # trexio.pytr is the SWIG module in both the old (top-level pytrexio) and
    # the new (trexio.pytrexio) package layouts
    if not trexio.pytr.trexio_has_backend(trexio.TREXIO_HDF5):
        msg = "the trexio Python module was built without the HDF5 back end"
        if _REQUIRED:
            pytest.fail(msg + ", and PSI4_TREXIO_VALIDATE_REQUIRED is set")
        pytest.skip(msg)


def _check_occupations(path, wfn):
    """What the validator cannot see: that each MO carries its own occupation.

    trexio-validate checks only that the occupations sum to the electron count,
    so marking a virtual orbital occupied and an occupied one empty passes it.
    For these ground states the occupied MOs of each spin must be the lowest.
    """
    with trexio.File(path, mode="r", back_end=trexio.TREXIO_HDF5) as f:
        energy = np.asarray(trexio.read_mo_energy(f))
        occ = np.asarray(trexio.read_mo_occupation(f))
        spin = np.asarray(trexio.read_mo_spin(f))

    restricted = wfn.same_a_b_orbs()
    for s, nel in ((0, wfn.nalpha()), (1, wfn.nbeta())):
        if restricted and s == 1:
            break
        e_s, o_s = energy[spin == s], occ[spin == s]
        assert np.all(np.diff(e_s) >= -1e-10), f"spin {s}: MOs not sorted by energy"
        nocc = int(np.count_nonzero(o_s))
        assert np.all(o_s[:nocc] > 0) and np.all(o_s[nocc:] == 0), \
            f"spin {s}: occupied MOs are not the lowest in energy: {o_s.tolist()}"
    if restricted:
        expected = [2.0] * wfn.nbeta() + [1.0] * (wfn.nalpha() - wfn.nbeta())
        assert occ[: wfn.nalpha()].tolist() == expected
    else:
        assert np.count_nonzero(occ[spin == 0]) == wfn.nalpha()
        assert np.count_nonzero(occ[spin == 1]) == wfn.nbeta()


def _chemists_eri(block, n):
    """Dense chemists' (ij|kl) from TrexIO's sparse physicists' <pq|rs> = (pr|qs)."""
    eri = np.zeros((n, n, n, n))
    for (p, q, r, s_), v in zip(block["indices"], block["values"]):
        i, j, k, l = p, r, q, s_
        for a, b in ((i, j), (j, i)):
            for c, d in ((k, l), (l, k)):
                eri[a, b, c, d] = v
                eri[c, d, a, b] = v
    return eri


def _check_energy(path, energy):
    """Rebuild the SCF energy from the file alone and compare with Psi4's.

    Uses only nucleus.repulsion, the MO core Hamiltonian, the MO ERIs and the
    occupations and spins, so it fails if any of them disagrees with what the
    calculation used -- in particular a core Hamiltonian missing the ECP term,
    which trexio-validate cannot check.
    """
    data = psi4.trexio_from_file(path)
    h = np.asarray(data["mo_1e"]["core_hamiltonian"])
    n = h.shape[0]
    g = _chemists_eri(data["mo_2e"], n)
    occ = np.asarray(data["mo"]["occupation"])
    spin = np.asarray(data["mo"]["spin"])
    if np.any(spin == 1):  # unrestricted: alpha and beta are separate MOs
        sets = [np.flatnonzero((spin == 0) & (occ > 0)), np.flatnonzero((spin == 1) & (occ > 0))]
    else:                  # restricted (RHF, ROHF): one MO set, 2 = both spins
        sets = [np.flatnonzero(occ >= 1), np.flatnonzero(occ >= 2)]

    e = float(data["nucleus"]["repulsion"])
    for occ_s in sets:
        e += np.trace(h[np.ix_(occ_s, occ_s)])
    everyone = np.concatenate(sets)
    e += 0.5 * sum(g[i, i, j, j] for i in everyone for j in everyone)  # Coulomb
    for occ_s in sets:                                                  # exchange, same spin
        e -= 0.5 * sum(g[i, j, j, i] for i in occ_s for j in occ_s)
    assert abs(e - energy) < 1.0e-8, f"energy from the file {e:.12f} != Psi4's {energy:.12f}"


@pytest.mark.parametrize("cid,geometry,options,method,eri,what", CASES, ids=[c[0] for c in CASES])
def test_trexio_validate(tmp_path, cid, geometry, options, method, eri, what):
    exe = _validator()
    _require_hdf5()

    psi4.core.clean()
    psi4.core.clean_options()
    psi4.core.clean_variables()
    psi4.geometry(geometry)
    options = dict(options)
    if options.get("basis") == "mixed_h2o":
        psi4.basis_helper("assign cc-pvdz\nassign H sto-3g", name="mixed_h2o")
        options.pop("basis")
    psi4.set_options({"scf_type": "pk", "e_convergence": 10, "d_convergence": 10, **options})

    try:
        energy, wfn = psi4.energy(method, return_wfn=True)
    except RuntimeError as err:
        if "libecpint addon not enabled" not in str(err):
            raise
        if _REQUIRED:
            pytest.fail("this Psi4 was built without libecpint, and PSI4_TREXIO_VALIDATE_REQUIRED is set")
        pytest.skip("this Psi4 was built without libecpint")

    path = str(tmp_path / f"{cid}.h5")
    psi4.trexio(wfn, path, back_end=trexio.TREXIO_HDF5,
                save_ao_integrals=True, save_mo_integrals=True, save_eri=eri, save_mo_eri=eri)

    skip = _NO_DIPOLES + ([] if eri else _ERI_CHECKS)
    if wfn.basisset().has_ECP():
        # trexio-validate does not recompute ECP integrals, and reports the core
        # Hamiltonian as unsupported once an ECP is present. _check_energy below
        # is what verifies it.
        skip += ["ao_1e_int_core_hamiltonian", "mo_1e_int_core_hamiltonian"]
    proc = subprocess.run([exe, "--require", "all", "--skip", ",".join(skip), path],
                          capture_output=True, text=True)
    report = proc.stdout + proc.stderr
    # 77 means nothing could be checked: for these files that is a failure too.
    assert proc.returncode == 0, f"{cid} ({what}): trexio-validate exit {proc.returncode}\n{report}"

    _check_occupations(path, wfn)
    if eri and method == "hf":
        _check_energy(path, energy)


def _ecp_wfn():
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.geometry("0 1\nH\nI 1 1.61\nsymmetry c1")
    psi4.set_options({"basis": "def2-svp", "scf_type": "pk", "e_convergence": 10, "d_convergence": 10})
    try:
        return psi4.energy("hf", return_wfn=True)
    except RuntimeError as err:
        if "libecpint addon not enabled" in str(err) and not _REQUIRED:
            pytest.skip("this Psi4 was built without libecpint")
        raise


def test_trexio_ecp_roundtrip(tmp_path):
    """An ECP file rebuilds the same ECP in Psi4: charges, core, and integrals."""
    _require_hdf5()
    energy, wfn = _ecp_wfn()
    path = str(tmp_path / "hi_ecp.h5")
    psi4.trexio(wfn, path, back_end=trexio.TREXIO_HDF5, save_ao_integrals=True)

    data = psi4.trexio_from_file(path)
    rebuilt = psi4.trexio_to_wavefunction(path)
    bs, bs0 = rebuilt.basisset(), wfn.basisset()
    mol, mol0 = rebuilt.molecule(), wfn.molecule()

    assert bs.has_ECP()
    assert bs.n_ecp_core() == bs0.n_ecp_core()
    np.testing.assert_allclose([mol.Z(i) for i in range(mol.natom())],
                               [mol0.Z(i) for i in range(mol0.natom())])
    assert abs(mol.nuclear_repulsion_energy() - mol0.nuclear_repulsion_energy()) < 1e-10

    # The ECP integrals recomputed from the rebuilt basis match the stored ones,
    # which were computed from the original: the ecp group lost nothing.
    from psi4.driver.p4util.trexio import _canonical_ao_normalization
    ao_norm = _canonical_ao_normalization(bs)
    recomputed = np.asarray(psi4.core.MintsHelper(bs).ao_ecp()) * np.outer(ao_norm, ao_norm)
    np.testing.assert_allclose(recomputed, np.asarray(data["ao_1e"]["ecp"]), atol=1e-12)


@pytest.mark.parametrize("cid", ["h2o_ccpvdz_c1", "hi_def2svp_ecp"])
def test_trexio_pyscf_reads_psi4(tmp_path, cid):
    """An independent reader reproduces Psi4's SCF energy from Psi4's file.

    The ECP round trip above cannot catch a convention error made symmetrically
    on write and read -- a radial power off by two, say. pyscf-forge's reader
    interprets the file on its own, and pyscf computes the ECP integrals with its
    own code. Water is the control: if it agrees and HI does not, the ECP is at
    fault. Runs only where pyscf-forge is installed.
    """
    pf_trexio = pytest.importorskip("pyscf.tools.trexio")
    scf = pytest.importorskip("pyscf.scf")
    _require_hdf5()
    geometry, options = {c[0]: (c[1], c[2]) for c in CASES}[cid]
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.geometry(geometry)
    psi4.set_options({"scf_type": "pk", "e_convergence": 10, "d_convergence": 10, **options})
    energy, wfn = psi4.energy("hf", return_wfn=True)
    path = str(tmp_path / f"{cid}.h5")
    psi4.trexio(wfn, path, back_end=trexio.TREXIO_HDF5)

    mol = pf_trexio.from_trexio(path)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-10
    assert abs(mf.kernel() - energy) < 1e-6
