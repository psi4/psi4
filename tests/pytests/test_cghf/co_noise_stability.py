#!/usr/bin/env python3
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2026 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE

"""Sample CGHF energies for CO at several bond lengths with machine-epsilon noise.

This script tests the numerical stability of the CGHF SCF by adding noise at the
level of machine epsilon to the MO coefficients immediately after the initial
guess.  For each bond length many independent noise realizations are run.  The
converged energies are summarized; a bond length is flagged as stable if the
standard deviation of the converged samples is below twice the energy convergence
threshold."""

import argparse
import sys
import numpy as np
import psi4
from psi4 import core
from psi4.driver.p4util.exceptions import SCFConvergenceError


# Options that make CO/cc-pVDZ converge with CGHF.  This mirrors the converging
# DIIS case from the CGHF test suite; the only difference is that we intercept
# the guess to perturb the MO coefficients.
DEFAULT_OPTIONS = {
    "basis": "cc-pVDZ",
    "reference": "cghf",
    "guess": "core",
    "scf_type": "direct",
    "df_scf_guess": False,
    "diis": False,
    "scf_initial_accelerator": "none",
    "orbital_optimizer_package": "internal",
    "maxiter": 20,
    "fail_on_maxiter": True,
    "e_convergence": 1.0e-6,
    "d_convergence": 1.0e-4,
}

def _co_molecule(r_co):
    """Build a neutral singlet CO molecule with the given C-O bond length (A)."""
    return f"""
        0 1
        C
        O 1 {r_co}
        symmetry c1
        """


def _add_noise_to_coefficients(C, rng, noise_scale):
    """Add machine-epsilon-scale complex noise to the MO coefficients in place.

    The perturbation is relative to the magnitude of each coefficient:

        C_{ij} <- C_{ij} * (1 + eps * (a + i b))

    where ``eps`` is machine epsilon scaled by ``noise_scale`` and ``a``,
    ``b`` are uniform random numbers in ``[-1, 1]``.  This is equivalent to
    flipping the least significant bits of the stored double-precision values.
    """
    eps = np.finfo(np.float64).eps * noise_scale
    if eps == 0.0:
        return

    for view in C.array_interface():
        arr = np.asarray(view)
        if arr.size == 0:
            continue

        shape = arr.shape
        re_pert = (rng.random(shape, dtype=np.float64) - 0.5) * 2.0
        im_pert = (rng.random(shape, dtype=np.float64) - 0.5) * 2.0
        arr *= (1.0 + eps * (re_pert + 1.0j * im_pert))


def _make_noisy_guess(original_guess, rng, noise_scale):
    """Return a CGHF.guess wrapper that perturbs the coefficients after the guess."""
    def noisy_guess(self):
        original_guess(self)
        _add_noise_to_coefficients(self.C(), rng, noise_scale)
        # Rebuild the density from the noisy orbitals so the first SCF step uses
        # the perturbed guess rather than the density from the unperturbed orbitals.
        self.form_D()
    return noisy_guess


def sample_bond_length(r_co, n_samples, noise_scale, rng, options):
    """Run ``n_samples`` CGHF calculations for a single CO bond length.

    Returns
    -------
    energies : np.ndarray
        Converged total energies (Eh).
    n_converged : int
        Number of samples that converged.
    n_failed : int
        Number of samples that failed to converge.
    """
    mol = psi4.geometry(_co_molecule(r_co))
    energies = []
    n_converged = 0
    n_failed = 0

    original_guess = core.CGHF.guess
    core.CGHF.guess = _make_noisy_guess(original_guess, rng, noise_scale)

    try:
        for _ in range(n_samples):
            psi4.core.clean()
            psi4.core.clean_options()
            psi4.set_options(options)
            try:
                e = psi4.energy("scf", molecule=mol)
                energies.append(e)
                n_converged += 1
            except SCFConvergenceError:
                n_failed += 1
            except Exception as exc:
                n_failed += 1
                print(f"    warning: unexpected failure at R={r_co}: {exc}", file=sys.stderr)
    finally:
        core.CGHF.guess = original_guess

    return np.array(energies, dtype=float), n_converged, n_failed


def main():
    parser = argparse.ArgumentParser(
        description="Sample CGHF energies for CO with machine-epsilon noise injected after the guess.",
    )
    parser.add_argument(
        "--n-samples", type=int, default=30,
        help="Number of independent noise realizations per bond length (default: 30).",
    )
    parser.add_argument(
        "--noise-scale", type=float, default=1.0,
        help="Scale factor for the machine-epsilon noise (default: 1.0).",
    )
    parser.add_argument(
        "--bond-lengths", type=float, nargs="+",
        default=[0.95, 1.0, 1.05, 1.10, 1.110, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40],
        help="CO bond lengths in Angstrom (default: a range around 1.110).",
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility (default: 42).",
    )
    parser.add_argument(
        "--threads", type=int, default=1,
        help="Number of threads to use for the SCF calculations (default: 1).",
    )
    parser.add_argument(
        "--output", default="cghf_co_noise_stability.log",
        help="Psi4 output file (default: cghf_co_noise_stability.log).",
    )
    args = parser.parse_args()

    if not psi4.extras.addons("einsums"):
        print("CGHF requires Einsums. Rebuild with -DENABLE_Einsums=ON.", file=sys.stderr)
        sys.exit(1)

    psi4.set_output_file(args.output, append=False)
    psi4.set_memory("2 GB")
    psi4.set_num_threads(args.threads, quiet=True)

    rng = np.random.default_rng(args.seed)
    options = dict(DEFAULT_OPTIONS)
    e_conv = options["e_convergence"]

    print("CGHF CO noise-stability sampling")
    print(f"  noise scale = {args.noise_scale} * eps")
    print(f"  samples per bond length = {args.n_samples}")
    print(f"  energy convergence threshold = {e_conv}")
    print(f"  stability criterion: std < {2 * e_conv}")
    print(f"  bond lengths (A) = {args.bond_lengths}")
    print(f"  threads = {args.threads}")
    print()

    header = f"{'R_CO':>8}  {'mean':>18}  {'std':>12}  {'min':>18}  {'max':>18}  {'conv':>5}  {'fail':>5}  {'stable':>7}"
    print(header)
    print("-" * len(header))

    for r in args.bond_lengths:
        energies, n_converged, n_failed = sample_bond_length(
            r, args.n_samples, args.noise_scale, rng, options
        )

        if n_converged == 0:
            print(f"{r:8.4f}  no converged samples ({n_failed} failed)")
            continue

        mean_e = energies.mean()
        min_e = energies.min()
        max_e = energies.max()
        if n_converged >= 2:
            std_e = energies.std(ddof=1)
            stable = std_e < 2 * e_conv
            std_str = f"{std_e:12.6e}"
            stable_str = "yes" if stable else "no"
        else:
            std_str = "         n/a"
            stable_str = "n/a"

        print(
            f"{r:8.4f}  {mean_e:18.10f}  {std_str}  {min_e:18.10f}  {max_e:18.10f}  "
            f"{n_converged:5d}  {n_failed:5d}  {stable_str:>7}"
        )


if __name__ == "__main__":
    main()
