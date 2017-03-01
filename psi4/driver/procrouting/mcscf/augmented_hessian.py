#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import numpy as np
import os

from psi4 import core
from psi4.driver import p4util
from psi4.driver.p4util.exceptions import *

# np.set_printoptions(precision=5, linewidth=200, threshold=2000, suppress=True)

def ah_iteration(mcscf_obj, tol=1e-3, max_iter=15, lindep=1e-14, print_micro=True):
    """
    Solve the generalized eigenvalue problem:
    | 0,  g.T | | 1/l | =   | 1/l |
    | g,  H/l | | X   | = e | X   |

    Where g is the gradient, H is the orbital Hessian, X is our orbital update step,
    and l is the eigenvalue.

    In some ways this is the subspace reduction of the full MCSCF Hessian where the
    CC part has been solved exactly. When this occurs the OC and CO elements collapse
    to the above and the CC Hessian becomes diagonally dominant.

    We can solve this through Davidson iterations where we condition the edges. It's the
    Pulay equations all over again, just iterative.

    Watch out for lambdas that are zero. Looking for the lambda that is ~1.

    """

    # Unpack information
    orb_grad = mcscf_obj.gradient()
    precon = mcscf_obj.H_approx_diag()
    approx_step = mcscf_obj.approx_solve()
    orb_grad_ssq = orb_grad.sum_of_squares()

    # Gears
    min_lambda = 0.3
    converged = False
    warning_neg = False
    warning_mult = False

    fullG = np.zeros((max_iter + 2, max_iter + 2))
    fullS = np.zeros((max_iter + 2, max_iter + 2))
    fullS[np.diag_indices_from(fullS)] = 1

    guesses = []
    sigma_list = []
    guesses.append(approx_step)
    sigma_list.append(mcscf_obj.compute_Hk(approx_step))

    if print_micro:
        core.print_out("\n                             Eigenvalue          Rel dE          dX \n")

    # Run Davidson look for lambda ~ 1
    old_val = 0
    for microi in range(1, max_iter + 1):

        # Gradient
        fullG[0, microi] = guesses[-1].vector_dot(orb_grad)
        for i in range(microi):
            fullG[i + 1, microi] = guesses[-1].vector_dot(sigma_list[i])
            fullS[i + 1, microi] = guesses[-1].vector_dot(guesses[i])

        fullG[microi] = fullG[:, microi]
        fullS[microi] = fullS[:, microi]

        wlast = old_val

        # Slice out relevant S and G
        S = fullS[:microi + 1, :microi + 1]
        G = fullG[:microi + 1, :microi + 1]

        # Solve Gv = lSv
        v, L = np.linalg.eigh(S)
        mask = v > (np.min(np.abs(v)) * 1.e-10)
        invL = L[:, mask] * (v[mask] ** -0.5)

        # Solve in S basis, rotate back
        evals, evecs = np.linalg.eigh(np.dot(invL.T, G).dot(invL))
        vectors = np.dot(invL, evecs)

        # Figure out the right root to follow
        if np.sum(np.abs(vectors[0]) > min_lambda) == 0:
            raise PsiException("Augmented Hessian: Could not find the correct root!\n"\
                               "Try starting AH when the MCSCF wavefunction is more converged.")

        if np.sum(np.abs(vectors[0]) > min_lambda) > 1 and not warning_mult:
            core.print_out("   Warning! Multiple eigenvectors found to follow. Following closest to \lambda = 1.\n")
            warning_mult = True

        idx = (np.abs(1 - np.abs(vectors[0]))).argmin()
        lam = abs(vectors[0, idx])
        subspace_vec = vectors[1:, idx]

        # Negative roots should go away?
        if idx > 0 and evals[idx] < -5.0e-6 and not warning_neg:
            core.print_out('   Warning! AH might follow negative eigenvalues!\n')
            warning_neg = True

        diff_val = evals[idx] - old_val
        old_val = evals[idx]

        new_guess = guesses[0].clone()
        new_guess.zero()
        for num, c in enumerate(subspace_vec / lam):
            new_guess.axpy(c, guesses[num])

        # Build estimated sigma vector
        new_dx = sigma_list[0].clone()
        new_dx.zero()
        for num, c in enumerate(subspace_vec):
            new_dx.axpy(c, sigma_list[num])

        # Consider restraints
        new_dx.axpy(lam, orb_grad)
        new_dx.axpy(old_val * lam, new_guess)

        norm_dx = (new_dx.sum_of_squares() / orb_grad_ssq) ** 0.5

        if print_micro:
            core.print_out("      AH microiter %2d   % 18.12e   % 6.4e   % 6.4e\n" % (microi, evals[idx],
                                    diff_val / evals[idx], norm_dx))

        if abs(old_val - wlast) < tol and norm_dx < (tol ** 0.5):
            converged = True
            break

        # Apply preconditioner
        tmp = precon.clone()
        val = tmp.clone()
        val.set(evals[idx])
        tmp.subtract(val)
        new_dx.apply_denominator(tmp)

        guesses.append(new_dx)
        sigma_list.append(mcscf_obj.compute_Hk(new_dx))

    if print_micro and converged:
        core.print_out("\n")
        #    core.print_out("      AH converged!       \n\n")

    #if not converged:
    #    core.print_out("      !Warning. Augmented Hessian did not converge.\n")

    new_guess.scale(-1.0)

    return converged, microi, new_guess

