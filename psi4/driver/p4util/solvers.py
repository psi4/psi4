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
import math
import time

import numpy as np

from psi4 import core
from .exceptions import ValidationError


"""
Generalized iterative solvers for Psi4.

"""


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
        core.print_out("     %4s %14s %12s  %6s  %6s\n" % ("Iter", "Residual RMS", "Max RMS", "Remain", "Time [s]"))
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

    resid = [(r_vec[x].sum_of_squares() / grad_dot[x])**0.5 for x in range(nrhs)]

    if printer:
        resid = printer(0, x_vec, r_vec)
    elif printlvl:
        # core.print_out('         CG Iteration Guess:    Rel. RMS = %1.5e\n' %  np.mean(resid))
        core.print_out("    %5s %14.3e %12.3e %7d %9d\n" % ("Guess", np.mean(resid), np.max(resid), len(z_vec),
                                                            time.time() - tstart))

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
            resid[x] = (r_vec[x].sum_of_squares() / grad_dot[x])**0.5

        # Print out or compute the resid function
        if printer:
            resid = printer(rot_iter + 1, x_vec, r_vec)

        # Figure out active updated active mask
        for x in active:
            if (resid[x] < rcond):
                active_mask[x] = False

        # Print out if requested
        if printlvl:
            core.print_out("    %5d %14.3e %12.3e %7d %9d\n" % (rot_iter + 1, np.mean(resid), np.max(resid),
                                                                sum(active_mask), time.time() - tstart))

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


def _diag_print_heading(solver_name, max_ss_size, nroot, e_tol, r_tol, maxiter, verbose=1):
    if verbose < 1:
        # no printing
        return
    else:
        # summary of options
        core.print_out("\n\n")
        core.print_out("   " + "{} options".format(solver_name).center(77) + "\n")
        core.print_out("\n  -----------------------------------------------------\n")
        core.print_out("    Maxiter                         = {:<5d}\n".format(maxiter))
        core.print_out("    Eigenvalue tolerance            = {:11.5e}\n".format(e_tol))
        core.print_out("    Eigenvector tolerance           = {:11.5e}\n".format(r_tol))
        core.print_out("    Max number of expansion vectors = {:<5d}\n".format(max_ss_size))
        core.print_out("\n")
        core.print_out("  => Iterations <=\n")
    if verbose == 1:
        core.print_out("  {}           {}      {}\n".format(" " * len(solver_name), "Max[delta val]", "Max[|R|]"))


def _diag_print_info(solver_name, info, verbose=1):
    if verbose < 1:
        # no printing
        return
    elif verbose == 1:
        # print iter  maxde max|R| conv/restart
        flags = []
        if info['collapse']:
            flags.append("Restart")
        if info['done']:
            flags.append("Converged")

        core.print_out("  {name} iter {ni:3d}:   {m_de:-11.5e} {m_r:12.5e} {flgs}\n".format(
            name=solver_name,
            ni=info['count'],
            m_de=np.max(info['delta_val']),
            m_r=np.amx(info['res_norm']),
            flgs="/".join(flags)))
    else:
        # print iter / ssdim folowed by de/|R| for each root
        core.print_out("  {name} iter {ni:3d}: {nv:4d} guess vectors\n".format(
            name=solver_name, ni=info['count'], nv=info['nvec']))
        for i, (de, rn) in enumerate(zip(info['delta_val'], info['res_norm'])):
            core.print_out("     {nr:2d}: {de:-11.5e} {rn:12.5e}\n".format(nr=i + 1, de=de, rn=rn))
            if info['done']:
                core.print_out("  Solver Converged! all roots\n\n")
            else:
                if info['collapse']:
                    core.print_out("  Subspace limits exceeded restarting\n\n")


def _diag_print_converged(solver_name, stats, vals, verbose=1, **kwargs):
    if verbose < 1:
        # no printing
        return
    if verbose == 1:
        # print values summary + number of iterations + # of "big" product evals
        core.print_out(" {} converged in {} iterations".format(solver_name, stats[-1]['count']))
        core.print_out("  Root #    eigenvalue\n")
        for (i, vi) in enumerate(vals):
            core.print_out("  {:^6}    {:20.12f}".format(i + 1, vi))
        max_nvec = max(istat['nvec'] for istat in stats)
        core.print_out("  Computed a total of {} Large products\n\n".format(stats[-1]['product_count']))


def davidson_solver(engine, guess, e_tol=1.0E-6, r_tol=1.0E-8, nroot=1, max_vecs_per_root=20, maxiter=100, verbose=1):
    """
    Solves for the lowest few eigenvalues and eigenvectors of a large problem emulated through an engine.


    If the large matrix `A` has dimension `{NxN}` and N is very large, and only a small number of roots, `k`
    are desired this algorithm is preferable to standard methods as uses on the order of `N * k` memory.
    One only needs to have the ability to compute the product of a times a vector.

    For non-hermitan `A` the basis of the algorithm breaks down. However in practice, for strongly diagonally-dominant `A`
    such as the similarity transformed hamiltonian in EOM-CC this algorithm commonly still used.

    Parameters
    -----------
    engine : object
       The engine drive all operations involving data structures that have one "large" dimension. The methods it is required
       to provide and their signatures are detailed below.
    guess : list {engine dependent}
       At least `nroot` initial expansion vectors
    e_tol : float
        Convergence tolerance for eigenvalues
    r_tol : float
        Convergence tolerance for residual vectors
    nroot : int
        Number of roots desired
    maxiter : int
        The maximum number of iterations
    verbose : int
        The amount of logging info to print (0 -> none, 1 -> some, >1 -> everything)

    Requirements of engine
    ----------------------
    orthogonalize_subspace(old_X, new_X) --> X, nx:
       old_X : list of `vectors`
           This set may be empty, if not, the vectors are mutually orthogonal.
       new_x : list of `trial vectors`
           The normalization or orthogonality of these vectors can't be assumed
       X : list of `vectors`
           A set of ortho-normal vectors constructed from old_X and new_X, need not have dimension = dim(old_X) + dim(new_X)
       nl : The number of vectors in the set.

    compute_products(X) -> AX :
       X : list of `vectors`
           As returned from `orthogonalize_subspace`
       Ax : list of `vectors`
           The product :math:`AX_{i}` for each `X_{i}` in `X`, in that order.

    subspace_project:(X, Ax) --> A
       X : list of `vectors` length d
           As returned from `orthogonalize_subspace`
       Ax : list of `vectors` length d
           As returned from `compute_products`
       A : :py:class:`numpy.ndarray` float {d, d}
           Each element of :math:`A_{ij}` is the dot product :math:`X_i(AX)_{j}`

    approximate_eigenvector(Ax, Vss_k) -> V
       Ax : list of `vectors` length d
           As returned from compute_products
       Vss_k : `numpy.ndarray` float {d, }
           The k'th eigenvector of the subspace projected problem
       V : single `vector`
           The approximate "true" eigenvector :math:`V_{k} = \sum_j (AX)_j\tilde{V}_k`

    residual(V, X, Vss_k, w_k) -> R
       V : single `vector`
           The approximate "true" eigenvector returned for `approximate_eigenvector`
       X : list of `vector` length d
           As returned from `orthogonalize_subspace`
       Vss_k : `numpy.ndarray` float {d,}
           As in `approximate_eigenvector`
       w_k : float
           The eigenvalue of the subspace problem associated with Vss_k
       R : single `vector`
           The residual vector :math:`R_{k} = V_k - \sum_j X_j\tilde{V}_k\tilde{\omega}_k`

    precondition_residual(R_k, w_k) -> new_X_k
       R_k : single `vector`
           As returned from `residual`
       w_k : float
           As in `residual`
       new_X_k : single `vector`
           This vector will be in the set passed to `orthogonalize_subspace` at the start of the next iteration.

     ..note:: The `vector` referred to here is intentionally vague, the solver does not care what it is and only
              holds individual or sets of them. In fact an individual `vector` could be split across two elements in a list,
              such as for different. This allowance is also why the `orthogonalize_subspace` function should return the dimension
              of the subspace. It is not required to be equal to `len(X)`.
    """
    # hard set for now --> make these options eventually
    imag_vec_tol = 1.0e-3
    imag_val_tol = 1.0e-3
    nl = len(guess)
    nk = nroot
    iter_info = dict(
        count=0,
        res_norm=np.zeros((nk)),
        val=np.zeros((nk)),
        delta_val=np.zeros((nk)),
        # conv defaults to true, and will be flipped when a non-conv root is hit
        done=True,
        nvec=nl,
        collapse=False,
        product_count=0,
    )

    print_name = "DavidsonSolver"
    if verbose != 0:
        core.print_out("\n  " + "Generalized Davidson Solver".center(53) + "\n")
        core.print_out("   " + "By Ruhee Dcunha".center(53) + "\n")

    max_ss_size = max_vecs_per_root * nk
    X_new = guess
    X = []

    stats = []

    _diag_print_heading(print_name, max_ss_size, nroot, r_tol, e_tol, maxiter, verbose)

    while iter_info['count'] < maxiter:
        # increment iteration/ save old vals
        iter_info['count'] += 1
        old_w = iter_info['val'].copy()
        # reset flags
        iter_info['collapse'] = False
        iter_info['done'] = True

        # expand subspace
        X, nl = engine.orthogonalize_subspace(X, X_new)
        X_new = []

        if nl >= max_ss_size:
            iter_info['collapse'] = True

        iter_info['nvec'] = nl

        Ax = engine.compute_products(X)
        iter_info['product_count'] += len(Ax)
        Ass = engine.subspace_project(Ax, X)

        # solve eigenvalue problem for subspace matrix
        w, Vss = np.linalg.eig(Ass)

        # sorting eigenvalues and corresponding eigenvectors & choose n lowest eigenvalue eigpairs
        idx = w.argsort()[:nk]
        Vss = Vss[:, idx]
        w = w[idx]
        imag_V = Vss.imag
        imag_w = w.imag
        Vss = Vss.real
        w = w.real

        if any(imag_w > imag_val_tol):
            raise Exception("DavidsonSolver: Im[eigenvalue] exceeded tolerance")
        if any(np.linalg.norm(imag_V, axis=1) > imag_vec_tol):
            raise Exception("DavidsonSolver: |Im[eigenvector]| exceeded tolerance")

        V = []
        for k in range(nk):
            # computes $\tilde{V}_k = sum_j Ax_j V^{ss}_k$
            V.append(engine.approximate_eigvector(Ax, Vss[:, k]))
            # $R^{v}_k = \tilde{V}_k - \sum_j \omega_k X_j V^{ss}_k$
            Rv = engine.residual(V[k], X, Vss[:, k], w[k])
            iter_info['res_norm'][k] = engine.vector_norm(Rv)
            iter_info['delta_val'][k] = abs(old_w[k] - w[k])
            iter_info['val'][k] = w[k]
            if (iter_info['res_norm'][k] < r_tol) and (iter_info['delta_val'][k] < e_tol):
                continue
            else:
                iter_info['done'] = False
                if iter_info['collapse']:
                    continue
                else:
                    X_new.append(engine.precondition_residual(Rv, shift=w[k]))

        _diag_print_info(print_name, iter_info, verbose)
        stats.append(iter_info.copy())
        if iter_info['done']:
            _diag_print_converged(print_name, stats, w, rvec=V, verbose=verbose)
            return w, V, stats
        if iter_info['collapse']:
            X = []
            X_new = V

    # If we get down here  we have exceeded max iterations without convergence, return none + stats
    return None, None, stats


def hamiltonian_solver(engine,
                       guess,
                       e_tol=1.0E-6,
                       r_tol=1.0E-8,
                       nroot=1,
                       max_vecs_per_root=20,
                       maxiter=100,
                       verbose=1):
    """
    Finds the smallest eigenvalues and associated right and left hand eigenvectors of a large real Hamiltonian eigenvalue problem
    emulated through an engine.

    Similar to the davidson algorithm but with special considerations that preserve the structure of the problem and lead to faster
    convergence.

    Parameters
    -----------
    engine : object
       The engine drive all operations involving data structures that have one "large" dimension. The methods it is required
       to provide and their signatures are detailed below.
    guess : list {engine dependent}
       At least `nroot` initial expansion vectors
    e_tol : float
        Convergence tolerance for eigenvalues
    r_tol : float
        Convergence tolerance for residual vectors
    nroot : int
        Number of roots desired
    maxiter : int
        The maximum number of iterations
    verbose : int
        The amount of logging info to print (0 -> none, 1 -> some, >1 -> everything)

    Requirements of engine
    ----------------------
    orthogonalize_subspace(old_X, new_X) --> X, nx:
       old_X : list of `vectors`
           This set may be empty, if not, the vectors are mutually orthogonal.
       new_x : list of `trial vectors`
           The normalization or orthogonality of these vectors can't be assumed
       X : list of `vectors`
           A set of ortho-normal vectors constructed from old_X and new_X, need not have dimension = dim(old_X) + dim(new_X)
       nl : The number of vectors in the set.

    compute_products(X) -> Px, Mx:
       X : list of `vectors`
           As returned from `orthogonalize_subspace`
       Px : list of `vectors`
           The product :math:`(A+B)X_{i}` for each `X_{i}` in `X`, in that order.
       Mx : list of `vectors`
           The product :math:`(A-B)X_{i}` for each `X_{i}` in `X`, in that order.

    subspace_project:(X, Ax) --> A
       X : list of `vectors` length d
           As returned from `orthogonalize_subspace`
       Ax : list of `vectors` length d
           As either of the returns from `compute_products`
       A : :py:class:`numpy.ndarray` float {d, d}
           Each element of :math:`A_{ij}` is the dot product :math:`X_i(AX)_{j}`

    approximate_eigenvector(Ax, Vss_k) -> V
       Ax : list of `vectors` length d
           As returned from `compute`_products
       Vss_k : `numpy.ndarray` float {d, }
           The k'th eigenvector of the subspace projected problem
       V : single `vector`
           The approximate "true" eigenvector :math:`V_{k} = \sum_j (AX)_j\tilde{V}_k`

    residual(V, X, Vss_k, w_k) -> R
       V : single `vector`
           The approximate "true" eigenvector returned for `approximate_eigenvector`
       X : list of `vector` length d
           As returned from `orthogonalize_subspace`
       Vss_k : `numpy.ndarray` float {d,}
           As in `approximate_eigenvector`
       w_k : float
           The eigenvalue of the subspace problem associated with Vss_k
       R : single `vector`
           The residual vector :math:`R_{k} = V_k - \sum_j X_j\tilde{V}_k\tilde{\omega}_k`

    precondition_residual(R_k, w_k) -> new_X_k
       R_k : single `vector`
           As returned from `residual`
       w_k : float
           As in `residual`
       new_X_k : single `vector`
           This vector will be in the set passed to `orthogonalize_subspace` at the start of the next iteration.

     ..note:: The `vector` referred to here is intentionally vague, the solver does not care what it is and only
              holds individual or sets of them. In fact an individual `vector` could be split across two elements in a list,
              such as for different. This allowance is also why the `orthogonalize_subspace` function should return the dimension
              of the subspace. It is not required to be equal to `len(X)`.

    """

    nk = nroot

    iter_info = dict(
        count=0,
        res_norm=np.zeros((nk)),
        val=np.zeros((nk)),
        delta_val=np.zeros((nk)),
        # conv defaults to true, and will be flipped when a non-conv root is hit
        conv=True,
        nvec=0,
        product_count=0,
    )
    ss_max = max_vecs_per_root * nk
    X_new = guess
    X = []

    stats = []
    print_name = "HamiltonianSolver"

    while iter_info['count'] < maxiter:
        # increment iteration/ save old vals
        iter_info['count'] += 1
        old_w = iter_info['val'].copy()
        # reset flags
        iter_info['collapse'] = False
        iter_info['done'] = True

        # expand subspace
        X, nl = engine.orthogonalize_subspace(X, X_new)
        X_new = []
        if nl >= ss_max:
            iter_info['collapse'] = True

        iter_info['nvec'] = nl
        # Step 2: compute (P)*bi and (M)*bi
        Px, Mx = engine.compute_products(X)
        iter_info['product_count'] += len(Px) * 2
        # Step 3: form Pss, Mss. The P/M matrices in the subspace
        Pss = engine.subspace_project(Px, X)
        Mss = engine.subspace_project(Mx, X)

        # Step 4: Hermitian Product (Subspace analog of M^{1/2} P M^{1/2})
        Mss_val, Mss_vec = np.linalg.eigh(Mss)
        if any(Mss_val < 0.0):
            raise Exception("A-B is not Positive Definite")
        Mss_half = np.einsum('ij,j,kj->ik', Mss_vec, np.sqrt(Mss_val), Mss_vec)
        Hss = np.einsum('ij,jk,km->im', Mss_half, Pss, Mss_half)

        # Step 5: diagonalize Hss -> w^2, Tss
        w, Tss = np.linalg.eigh(Hss)
        w = np.sqrt(w)
        idx = np.argsort(w)
        w = w[idx]
        Tss = Tss[:, idx]

        #Step 6a: extract Rss = M^{1/2}Tss
        Rss = np.dot(Mss_half, Tss)
        #Step 6b: extract Lss = Pss Rss * w^-1
        Lss = np.dot(Pss, Rss)
        Lss = np.einsum('ij,j->ij', Lss, np.divide(1.0, w))

        # store best R/L eigvals
        L = []
        R = []
        for k in range(nk):
            L.append(engine.approximate_eigvector(Px, Rss[:, k]))
            R.append(engine.approximate_eigvector(Mx, Lss[:, k]))
            WL = engine.residual(L[k], X, Lss[:, k], w[k])
            WR = engine.residual(R[k], X, Rss[:, k], w[k])
            norm_r = engine.vector_norm(WR)
            norm_l = engine.vector_norm(WL)
            iter_info['res_norm'][k] = norm_r + norm_l
            iter_info['delta_val'][k] = np.abs(old_w[k] - w[k])
            iter_info['val'][k] = w[k]
            if (iter_info['res_norm'][k] < r_tol) and (iter_info['delta_val'][k] < e_tol):
                continue
            else:
                iter_info['done'] = False
                if iter_info['collapse']:
                    # don't bother preconditioning, we are not going to use them
                    continue
                else:
                    # precondition residuals
                    X_new.append(engine.precondition_residual(WL, shift=w[k]))
                    X_new.append(engine.precondition_residual(WR, shift=w[k]))

        _diag_print_info(print_name, iter_info, verbose)
        stats.append(iter_info.copy())
        if iter_info['done']:
            _diag_print_converged(print_name, stats, w, rvec=R, lvec=L, verbose=verbose)
            return w[:nk], R, L, stats
        if iter_info['collapse']:
            # list add, not vector add, collapsed guess space will be 2*nroots in size(nroot L + nroot R)
            X = []
            X_new = R + L

    # if we get down here we have exceeded max iterations without convergence, return none + stats
    return None, None, None, stats
