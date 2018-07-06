#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2018 The Psi4 Developers.
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
#

from __future__ import print_function
"""
Generalized iterative solvers for Psi4.

"""
import time

import numpy as np

from psi4 import core
from .exceptions import *


def cg_solver(rhs_vec, hx_function, preconditioner, guess=None, printer=None, printlvl=1, maxiter=20, rcond=1.e-6):
    """
    Solves the Ax = b linear equations via Conjugate Gradient. The `A` matrix must be a hermitian, positive definite matrix.

    Parameters
    ----------
    rhs_vec : list of :py:class:`~psi4.core.Matrix`
        The RHS vector in the Ax=b equation.
    hx_function : function
        Takes in a list of :py:class:`~psi4.core.Matrix` objects and a mask of active indices. Returns the Hessian-vector product.
    preconditioner : function
        Takes in a list of :py:class:`~psi4.core.Matrix` objects and a mask of active indices. Returns the preconditioned value.
    guess : list of :py:class:`~psi4.core.Matrix`, optional
        Starting vectors, if None use a preconditioner(rhs) guess
    printer : function, optional
        Takes in a list of current x and residual vectors and provides a print function. This function can also
        return a value that represents the current residual.
    printlvl : int, optional
        The level of printing provided by this function.
    maxiter : int, optional
        The maximum number of iterations this function will take.
    rcond : float, optional
        The residual norm for convergence.

    Returns
    -------
    ret : tuple, list of :py:class:`~psi4.core.Matrix`
        Returns the solved `x` vectors and `r` vectors.

    Notes
    -----
    This is a generalized cg solver that can also take advantage of solving multiple RHS's simultaneously when
    it is advantageous to do so.

    Examples
    --------



    """

    tstart = time.time()
    if printlvl:
        core.print_out("\n   -----------------------------------------------------\n")
        core.print_out("   " + "Generalized CG Solver".center(52) + "\n")
        core.print_out("   " + "by Daniel. G. A. Smith".center(52) + "\n")
        core.print_out("   -----------------------------------------------------\n")
        core.print_out("    Maxiter             = %11d\n" % maxiter)
        core.print_out("    Convergence         = %11.3E\n" % rcond)
        core.print_out("    Number of equations = %11ld\n\n" % len(rhs_vec))
        core.print_out("     %4s %14s %12s  %6s  %6s\n" %
                       ("Iter", "Residual RMS", "Max RMS", "Remain", "Time [s]"))
        core.print_out("   -----------------------------------------------------\n")

    nrhs = len(rhs_vec)
    active_mask = [True for x in range(nrhs)]

    # Start function
    if guess is None:
        x_vec = preconditioner(rhs_vec, active_mask)
    else:
        if len(guess) != len(rhs_vec):
            raise ValidationError("CG Solver: Guess vector length does not match RHS vector length.")
        x_vec = [x.clone() for x in guess]

    Ax_vec = hx_function(x_vec, active_mask)

    # Set it up
    r_vec = []  # Residual vectors
    for x in range(nrhs):
        tmp_r = rhs_vec[x].clone()
        tmp_r.axpy(-1.0, Ax_vec[x])
        r_vec.append(tmp_r)

    z_vec = preconditioner(r_vec, active_mask)
    p_vec = [x.clone() for x in z_vec]

    # First RMS
    grad_dot = [x.sum_of_squares() for x in rhs_vec]

    resid = [(r_vec[x].sum_of_squares() / grad_dot[x]) ** 0.5 for x in range(nrhs)]

    if printer:
        resid = printer(0, x_vec, r_vec)
    elif printlvl:
        # core.print_out('         CG Iteration Guess:    Rel. RMS = %1.5e\n' %  np.mean(resid))
        core.print_out("    %5s %14.3e %12.3e %7d %9d\n" % (
            "Guess", np.mean(resid), np.max(resid), len(z_vec), time.time() - tstart))

    rms = np.mean(resid)
    rz_old = [0.0 for x in range(nrhs)]
    alpha = [0.0 for x in range(nrhs)]
    active = np.where(active_mask)[0]

    # CG iterations
    for rot_iter in range(maxiter):

        # Build old RZ so we can discard vectors
        for x in active:
            rz_old[x] = r_vec[x].vector_dot(z_vec[x])

        # Build Hx product
        Ap_vec = hx_function(p_vec, active_mask)

        # Update x and r
        for x in active:
            alpha[x] = rz_old[x] / Ap_vec[x].vector_dot(p_vec[x])
            if np.isnan(alpha)[0]:
                core.print_out("CG: Alpha is NaN for vector %d. Stopping vector." % x)
                active_mask[x] = False
                continue

            x_vec[x].axpy(alpha[x], p_vec[x])
            r_vec[x].axpy(-alpha[x], Ap_vec[x])
            resid[x] = (r_vec[x].sum_of_squares() / grad_dot[x]) ** 0.5


        # Print out or compute the resid function
        if printer:
            resid = printer(rot_iter + 1, x_vec, r_vec)

        # Figure out active updated active mask
        for x in active:
            if (resid[x] < rcond):
                active_mask[x] = False

        # Print out if requested
        if printlvl:
            core.print_out("    %5d %14.3e %12.3e %7d %9d\n" % (
                rot_iter + 1, np.mean(resid), np.max(resid), sum(active_mask), time.time() - tstart))


        active = np.where(active_mask)[0]

        if sum(active_mask) == 0:
            break

        # Update p
        z_vec = preconditioner(r_vec, active_mask)
        for x in active:
            beta = r_vec[x].vector_dot(z_vec[x]) / rz_old[x]
            p_vec[x].scale(beta)
            p_vec[x].axpy(1.0, z_vec[x])

    if printlvl:
        core.print_out("   -----------------------------------------------------\n")

    return x_vec, r_vec


class DIIS(object):

    """
    An object to assist in the DIIS extrpolation procedure.
    """

    def __init__(self, max_vec=6, removal_policy="OLDEST"):
        """
        An object to assist in the DIIS extrpolation procedure.

        Parameters
        ----------
        max_vect : int, optional
            The maximum number of error and state vectors to hold. These are pruned based off the removal policy.
        removal_policy : {"OLDEST", "LARGEST"}, optional
            How the state and error vectors are removed once at the maximum. OLDEST will remove the oldest vector while
            largest will remove the residual with the largest RMS value.

        """
        self.error = []
        self.state = []
        self.max_vec = max_vec
        self.removal_policy = removal_policy.upper()

        if self.removal_policy not in ["LARGEST", "OLDEST"]:
            raise ValidationError("DIIS: removal_policy must either be oldest or largest.")

    def add(self, state, error):
        """
        Adds a DIIS state and error vector to the DIIS object.

        state : :py:class:`~psi4.core.Matrix`
            The current state vector.
        error : :py:class:`~psi4.core.Matrix`
            The current error vector.

        """
        self.error.append(error.clone())
        self.state.append(state.clone())

    def extrapolate(self, out=None):
        """
        Extrapolates next state vector from the current set of state and error vectors.

        Parameters
        ----------
        out : :py:class:`~psi4.core.Matrix`, optional
            A array in which to place the next state vector.

        Returns
        -------
        ret : :py:class:`~psi4.core.Matrix`
            Returns the next state vector.

        """

        # Limit size of DIIS vector
        diis_count = len(self.state)

        if diis_count == 0:
            raise ValidationError("DIIS: No previous vectors.")
        if diis_count == 1:
            return self.state[0]

        if diis_count > self.max_vec:

            if self.removal_policy == "OLDEST":
                pos = 0
            else:
                pos = np.argmax([x.rms() for x in self.error])

            del self.state[pos]
            del self.error[pos]
            diis_count -= 1

        # Build error matrix B
        B = np.empty((diis_count + 1, diis_count + 1))
        B[-1, :] = 1
        B[:, -1] = 1
        B[-1, -1] = 0
        for num1, e1 in enumerate(self.error):
            B[num1, num1] = e1.vector_dot(e1)
            for num2, e2 in enumerate(self.error):
                if num2 >= num1:
                    continue
                val = e1.vector_dot(e2)
                B[num1, num2] = B[num2, num1] = val

        # Build residual vector
        resid = np.zeros(diis_count + 1)
        resid[-1] = 1

        # Solve pulay equations

        # Yea, yea this is unstable make it stable
        iszero = np.any(np.diag(B)[:-1] <= 0.0)
        if iszero:
            S = np.ones((diis_count + 1))
        else:
            S = np.diag(B).copy()
            S[:-1] **= -0.5
            S[-1] = 1

        # Then we gotta do a custom inverse
        B *= S[:, None] * S
        invB = core.Matrix.from_array(B)
        invB.power(-1.0, 1.e-12)

        ci = np.dot(invB, resid)
        ci *= S

        # combination of previous fock matrices
        if out is None:
            out = core.Matrix("DIIS result", self.state[0].rowdim(), self.state[1].coldim())
        else:
            out.zero()

        for num, c in enumerate(ci[:-1]):
            out.axpy(c, self.state[num])

        return out


# Generalized Davidson solver

def davidson_solver(ax_function, preconditioner, guess, e_conv=1.0E-8,
        r_conv=None, no_eigs=1, max_vecs_per_root=20, maxiter=100):
    """
    Solves for the lowest few eigenvalues and eigenvectors of the given real symmetric matrix

    Parameters
    -----------
    ax_function : function
        Takes in a guess vector and returns the Matrix-vector product.
    preconditioner : function
        Takes in a list of :py:class:`~psi4.core.Matrix` objects and a mask of active indices. Returns the preconditioned value.
    guess : list of :py:class:`~psi4.core.Matrix`
        Starting vectors, required
    e_conv : float
        Convergence tolerance for eigenvalues
    r_conv : float
        Convergence tolerance for residual vectors
    no_eigs : int
        Number of eigenvalues needed
    maxiter : int
        The maximum number of iterations this function will take.

    Returns
    -----------

    Notes
    -----------

    Examples
    -----------

    """

    print("Starting Davidson algo:")
    print("options:")
    for k,v in dict(e_conv=e_conv, r_conv=r_conv, no_eigs=no_eigs,
            max_vecs_per_root = max_vecs_per_root, maxiter=maxiter).items():
        if v == None:
            v = 'None'
        print("{:<25}{:.>100}".format(k,v))
    if r_conv == None:
        r_conv = e_conv * 100
    d_tol = 1.0E-8

    # using the shape of the guess vectors to set the dimension of the matrix
    N = guess.shape[0]

    #sanity check, guess subspace must be at least equal to number of eigenvalues
    nli = guess.shape[1]
    if nli < no_eigs:
        raise ValueError("Not enough guess vectors provided!")

    nl = nli
    converged=False
    count = 0
    A_w_old = np.ones(no_eigs)
    max_ss_size = no_eigs * max_vecs_per_root
    B = guess

    ### begin loop
    conv_roots = [False] * no_eigs
    while count < maxiter:
        # active_mask = [True for x in range(nl)] # never used later
        # Apply QR decomposition on B to orthogonalize the new vectors wrto all other subspace vectors
        ## orthogonalize preconditioned residuals against all other vectors in the search subspace
        B, r = np.linalg.qr(B)
        nl = B.shape[1]

        print("Davidson: Iter={:<3} nl = {:<4}".format(count, nl))
        # compute sigma vectors corresponding to the new vectors sigma_i = A B_i
        sigma = np.zeros_like(B)
        for i in range(nl):
            sigma[:,i] = ax_function(B[:,i])

        # compute subspace matrix A_b = Btranspose sigma
        A_b = np.dot(B.T, sigma)

        # solve eigenvalue problem for subspace matrix; choose n lowest eigenvalue eigpairs
        A_w, A_v = np.linalg.eig(A_b)

        # sorting eigenvalues and corresponding eigenvectors
        idx = A_w.argsort()[:no_eigs]
        A_v_imag = A_v[:, idx].imag
        A_w_imag = A_w[idx].imag
        A_v = A_v[:, idx].real
        A_w = A_w[idx].real

        # Print warning if complex parts are too large
        for i,w_i in enumerate(A_w_imag):
            if abs(w_i) > 1.0e-10:
                print("*** WARNING***")
                print("    Root {:>5}: |Imag[A_w[{}]]| > 1.0E-10!".format(i))
                print("              : |Imag[A_v[{}]]| = {:.2E}".format(
                    i, np.linalg.norm(A_v_imag[:, i])))


        # here, check if no residuals > max no residuals, if so, collapse subspace
        if nl >= max_ss_size:
            print("Subspace too big. Collapsing.\n")
            B = np.dot(B, A_v)
            continue
        # else, build residual vectors
        ## residual_i = sum(j) sigma_j* eigvec_i - eigval_i * B_j * eigvec_i
        norm = np.zeros(no_eigs)
        new_Bs = []
        # Only need a residual for each desired root, not one for each guess
        force_collapse = False
        for i in range(no_eigs):
            if conv_roots[i]:
                print("    ROOT {:<3}: CONVERGED!".format(i))
                continue
            residual = np.dot(sigma, A_v[:, i]) - A_w[i] * np.dot(B, A_v[:,i])

            # check for convergence by norm of residuals
            norm[i] = np.linalg.norm(residual)
            # apply the preconditioner
            precon_resid = preconditioner(residual, i, A_w[i])
            de = abs(A_w_old[i] - A_w[i])
            print("    ROOT {:<3} de = {:<10.6f} |r| = {:<10.6f}".format(i, de,
                norm[i]))
            conv_roots[i] = (de < e_conv) and (norm[i] < r_conv)
            if conv_roots[i]:
                force_collapse = True
            else:
                new_Bs.append(precon_resid)

        # check for convergence by diff of eigvals and residual norms
        # r_norm = np.linalg.norm(norm)
        # eig_norm = np.linalg.norm(A_w - A_w_old)
        A_w_old = A_w
        #if( r_norm < r_conv) and (eig_norm < e_conv):
        if all(conv_roots):
            converged = True
            print("Davidson converged at iteration number {}".format(count))
            print("{:<3}|{:^20}|".format('count','eigenvalues'))
            print("{:<3}|{:^20}|".format('='*3, '='*20))
            retvecs = np.dot(B,A_v)
            for i, val in enumerate(A_w):
                print("{:<3}|{:<20.12f}".format(i, val))
            return A_w, retvecs
        else:
            if force_collapse:
                B = np.dot(B, A_v)
            n_left_to_converge = np.count_nonzero(np.logical_not(conv_roots))
            n_converged = np.count_nonzero(conv_roots)
            max_ss_size = n_converged + (n_left_to_converge * max_vecs_per_root)
            B = np.column_stack(tuple(B[:, i] for i in range(B.shape[-1])) + tuple(new_Bs))
        count += 1

    if not converged:
        print("Davidson did not converge. Max iterations exceeded.")
        return None, None
