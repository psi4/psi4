"""DFT functional energies across the internal (CPU) and cuEST (GPU) paths.

This is the PsiAPI companion to the ctest ``dft1`` case, not a replacement for
it.  ``tests/dft1/input.dat`` is a longstanding regression pinned to
``SCF_TYPE DIRECT``, and it is left exactly as it is; an earlier attempt to
parametrize that input over cuEST was reverted because the two engines cannot
share its settings (see below), which forced the input to carry two disjoint
reference sets and a conditional functional list.

Here both legs instead run *identical* settings, chosen to give the two engines
the best chance of producing the same number:

  * ``SCF_TYPE DF``.  cuESTJK is density-fitted only -- ``jk.cc`` builds it for
    MEM_DF/DISK_DF -- so the DIRECT of the ctest case is not available.
  * ``DFT_NUCLEAR_SCHEME STRATMANN``.  cuEST reproduces internal Psi4 to ~1e-9
    only under Stratmann nuclear partitioning.  With the default scheme the
    meta-GGA drifts by ~5e-5, which is grid partitioning rather than a defect
    in either engine.
  * ``E_CONVERGENCE``/``D_CONVERGENCE`` of 1e-9, tighter than the ctest case's
    1e-8, so that the SCF is converged well below the tolerance at which the
    two engines are compared against each other.

Every case is checked twice: against the stored reference for its own engine at
a tight tolerance, and against the *internal* reference at a looser one.  The
first is the regression check proper; the second is the cross-engine check, and
is the reason the settings above are shared rather than tuned per engine.  On
the internal leg the two comparisons coincide, which is harmless and keeps the
test body uniform.
"""

import pytest

import psi4

from addons import using
from utils import compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.scf, pytest.mark.dft, pytest.mark.df]


# Geometries as in tests/dft1/input.dat.
_GEOMS = {
    "water": """
0 1
O 0.000000000000  0.000000000000 -0.068516219310
H 0.000000000000 -0.790689573744  0.543701060724
H 0.000000000000  0.790689573744  0.543701060724
no_com
no_reorient
""",
    "water_cation": """
1 2
O
H 1 1.0
H 1 1.0 2 104.5
""",
    "water_anion": """
-1 2
O
H 1 1.0
H 1 1.0 2 104.5
""",
}

# tag -> (label, geometry, reference, guess)
_BLOCKS = {
    "rks_neutral": ("RKS  0 1", "water", "rks", None),
    "uks_neutral": ("UKS  0 1", "water", "uks", None),
    "uks_cation": ("UKS  1 2", "water_cation", "uks", None),
    "uks_anion": ("UKS -1 2", "water_anion", "uks", "sad"),
}

_FUNCTIONALS = ["svwn", "b3lyp", "wB97", "wB97X", "b86bpbe", "pw86pbe", "m05", "spw92"]

# cuESTJK::compute_JK() computes J and K but never wK, and nothing consults
# do_wK_, so a range-separated functional does not raise on the cuEST leg -- it
# silently converges 1-2 Eh off (wB97 on the neutral case: -73.264 rather than
# -75.319).  The gradient code already sidesteps this (scf_grad.cc requires
# !functional_->is_x_lrc() before taking the cuEST JK-gradient path); the energy
# path has no such check.
#
# These are xfail(strict=True) rather than skipped so the gap stays visible in
# the test report and, once cuEST grows wK support, the resulting XPASS fails
# the suite and says so.  At that point drop this tuple and fill in the two
# missing cuEST references.
_NEEDS_WK = ("wB97", "wB97X")

_OPTIONS = {
    "basis": "sto-3g",
    "scf_type": "df",
    "df_basis_scf": "def2-universal-JKFIT",
    "dft_nuclear_scheme": "stratmann",
    "dft_spherical_points": 302,
    "dft_radial_points": 99,
    "e_convergence": 9,
    "d_convergence": 9,
}

# Reference energies, keyed (block tag, functional), each set produced by the
# engine it is filed under with the settings above.
_REFS = {
    "internal": {
        # RKS  0 1
        ("rks_neutral", "svwn"):         -74.9361457674,
        ("rks_neutral", "b3lyp"):        -75.3200834739,
        ("rks_neutral", "wB97"):         -75.3189727239,
        ("rks_neutral", "wB97X"):        -75.3132187890,
        ("rks_neutral", "b86bpbe"):      -75.2981206541,
        ("rks_neutral", "pw86pbe"):      -75.3503184102,
        ("rks_neutral", "m05"):          -75.3214214214,
        ("rks_neutral", "spw92"):        -74.7371639178,
        # UKS  0 1
        ("uks_neutral", "svwn"):         -74.9361457674,
        ("uks_neutral", "b3lyp"):        -75.3200834739,
        ("uks_neutral", "wB97"):         -75.3189727239,
        ("uks_neutral", "wB97X"):        -75.3132187890,
        ("uks_neutral", "b86bpbe"):      -75.2981206541,
        ("uks_neutral", "pw86pbe"):      -75.3503184102,
        ("uks_neutral", "m05"):          -75.3214214214,
        ("uks_neutral", "spw92"):        -74.7371639178,
        # UKS  1 2
        ("uks_cation", "svwn"):          -74.5658273814,
        ("uks_cation", "b3lyp"):         -74.9678744373,
        ("uks_cation", "wB97"):          -74.9713472832,
        ("uks_cation", "wB97X"):         -74.9636555640,
        ("uks_cation", "b86bpbe"):       -74.9477314548,
        ("uks_cation", "pw86pbe"):       -74.9997520073,
        ("uks_cation", "m05"):           -74.9703623086,
        ("uks_cation", "spw92"):         -74.3874998062,
        # UKS -1 2
        ("uks_anion", "svwn"):           -74.4319601928,
        ("uks_anion", "b3lyp"):          -74.8158925335,
        ("uks_anion", "wB97"):           -74.8016635922,
        ("uks_anion", "wB97X"):          -74.7982781550,
        ("uks_anion", "b86bpbe"):        -74.7908482778,
        ("uks_anion", "pw86pbe"):        -74.8450791750,
        ("uks_anion", "m05"):            -74.8016980089,
        ("uks_anion", "spw92"):          -74.2143762266,
    },
    # No wB97/wB97X entries: see _NEEDS_WK above.
    "cuest": {
        # RKS  0 1
        ("rks_neutral", "svwn"):         -74.9361457674,
        ("rks_neutral", "b3lyp"):        -75.3200834739,
        ("rks_neutral", "b86bpbe"):      -75.2981206541,
        ("rks_neutral", "pw86pbe"):      -75.3503184102,
        ("rks_neutral", "m05"):          -75.3214214214,
        ("rks_neutral", "spw92"):        -74.7371639179,
        # UKS  0 1
        ("uks_neutral", "svwn"):         -74.9361457674,
        ("uks_neutral", "b3lyp"):        -75.3200834739,
        ("uks_neutral", "b86bpbe"):      -75.2981206541,
        ("uks_neutral", "pw86pbe"):      -75.3503184102,
        ("uks_neutral", "m05"):          -75.3214214214,
        ("uks_neutral", "spw92"):        -74.7371639179,
        # UKS  1 2
        ("uks_cation", "svwn"):          -74.5658273814,
        ("uks_cation", "b3lyp"):         -74.9678744373,
        ("uks_cation", "b86bpbe"):       -74.9477314548,
        ("uks_cation", "pw86pbe"):       -74.9997520074,
        ("uks_cation", "m05"):           -74.9703623086,
        ("uks_cation", "spw92"):         -74.3874998062,
        # UKS -1 2
        ("uks_anion", "svwn"):           -74.4319601928,
        ("uks_anion", "b3lyp"):          -74.8158925336,
        ("uks_anion", "b86bpbe"):        -74.7908482778,
        ("uks_anion", "pw86pbe"):        -74.8450791750,
        ("uks_anion", "m05"):            -74.8016980089,
        ("uks_anion", "spw92"):          -74.2143762266,
    },
}

# Tolerances, from measurement over this 4-block x 8-functional matrix on an
# L40S (sm_89):
#
#   cuEST run-to-run, same GPU        <= 9.9e-12
#   internal vs cuEST                 <= 3.7e-11
#   rounding, refs stored to 10 dp    <= 4.7e-11
#
# so a comparison against a stored reference carries ~6e-11 of unavoidable
# noise and a cross-engine one ~1e-10.  Note the 3.7e-11 is what the tightened
# convergence above buys: at E/D_CONVERGENCE 1e-8 the same spread was 8.6e-10.
_ATOL_OWN = 1.0e-9

# Looser, as the module docstring explains: this one crosses engines, and the
# spread above was measured on a single GPU architecture.
_ATOL_CROSS = 1.0e-8


def _run_case(tag, func, cuest):
    """Run one (block, functional) case on the requested engine."""
    label, geom, reference, guess = _BLOCKS[tag]

    mol = psi4.geometry(_GEOMS[geom])

    psi4.set_options(_OPTIONS)
    psi4.set_options({"reference": reference, "use_cuest": cuest})
    if guess:
        psi4.set_options({"guess": guess})

    return psi4.energy("scf", dft_functional=func, molecule=mol)


def _missing_wk(mode, func):
    """True where the engine cannot do this functional's range separation."""
    return mode == "cuest" and func in _NEEDS_WK


def _cases():
    for mode in ("internal", "cuest"):
        mode_marks = [] if mode == "internal" else [*using("cuest"), *using("cuda_cc8")]

        for tag in _BLOCKS:
            for func in _FUNCTIONALS:
                marks = list(mode_marks)
                if _missing_wk(mode, func):
                    marks.append(
                        pytest.mark.xfail(
                            strict=True,
                            reason="cuESTJK computes J and K but not wK, so range-separated "
                            "functionals converge to a wrong answer instead of failing",
                        )
                    )

                yield pytest.param(mode, tag, func, marks=marks, id=f"{mode}-{tag}-{func.lower()}")


@pytest.mark.parametrize("mode,tag,func", list(_cases()))
def test_dft1(mode, tag, func):
    label = _BLOCKS[tag][0]
    val = _run_case(tag, func, cuest=(mode == "cuest"))

    name = f"{label} {func.upper():>7}"

    # The wK cases deliberately have no cuEST reference -- storing one would be
    # filing a known-wrong number as expected behaviour. Their xfail therefore
    # has to come from the cross-engine comparison below, which is the check
    # that will start passing (and so XPASS, and so fail) once wK lands. Every
    # other case is pinned to its own engine's reference as well.
    if not _missing_wk(mode, func):
        assert compare_values(_REFS[mode][tag, func], val, _ATOL_OWN, f"{name} vs {mode} ref")

    assert compare_values(_REFS["internal"][tag, func], val, _ATOL_CROSS, f"{name} vs internal")
