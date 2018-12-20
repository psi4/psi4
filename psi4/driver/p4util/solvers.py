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


def _diag_print_heading(title_lines, solver_name, max_ss_size, nroot, e_tol, r_tol, maxiter, verbose=1):
    """Print a message to the output file when the solver has processed all options and is ready to begin"""
    if verbose < 1:
        # no printing
        return
    # show title if not silent
    core.print_out("\n\n")
    core.print_out("\n".join([x.center(77) for x in title_lines]))
    core.print_out("\n")
    if verbose > 1:
        # summarize options for verbose
        core.print_out("   " + "{} options".format(solver_name).center(77) + "\n")
        core.print_out("\n  -----------------------------------------------------\n")
        core.print_out("    Maxiter                         = {:<5d}\n".format(maxiter))
        core.print_out("    Eigenvalue tolerance            = {:11.5e}\n".format(e_tol))
        core.print_out("    Eigenvector tolerance           = {:11.5e}\n".format(r_tol))
        core.print_out("    Max number of expansion vectors = {:<5d}\n".format(max_ss_size))
        core.print_out("\n")
    # show iteration info headings if not silent
    core.print_out("  => Iterations <=\n")
    if verbose == 1:
        # default printing one line per iter max delta value and max residual norm
        core.print_out("  {}           {}      {}\n".format(" " * len(solver_name), "Max[D[value]]", "Max[|R|]"))
    else:
        # verbose printing, value, delta, and |R| for each root
        core.print_out("    {}       {}      {}      {}\n".format(" " * len(solver_name), "value", "D[value]", "|R|"))


def _diag_print_info(solver_name, info, verbose=1):
    """Print a message to the output file at each iteration"""
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
            m_r=np.max(info['res_norm']),
            flgs="/".join(flags)))
    else:
        # print iter / ssdim folowed by de/|R| for each root
        core.print_out("  {name} iter {ni:3d}: {nv:4d} guess vectors\n".format(
            name=solver_name, ni=info['count'], nv=info['nvec']))
        for i, (e, de, rn) in enumerate(zip(info['val'], info['delta_val'], info['res_norm'])):
            core.print_out("     {nr:2d}: {s:} {e:-11.5f} {de:-11.5e} {rn:12.5e}\n".format(
                nr=i + 1, s=" " * (len(solver_name) - 8), e=e, de=de, rn=rn))
        if info['done']:
            core.print_out("  Solver Converged! all roots\n\n")
        elif info['collapse']:
            core.print_out("  Subspace limits exceeded restarting\n\n")


def _diag_print_converged(solver_name, stats, vals, verbose=1, **kwargs):
    """Print a message to the output file when the solver is converged."""
    if verbose < 1:
        # no printing
        return
    if verbose >= 1:
        # print values summary + number of iterations + # of "big" product evals
        core.print_out(" {} converged in {} iterations\n".format(solver_name, stats[-1]['count']))
        core.print_out("  Root #    eigenvalue\n")
        for (i, vi) in enumerate(vals):
            core.print_out("  {:^6}    {:20.12f}\n".format(i + 1, vi))
        max_nvec = max(istat['nvec'] for istat in stats)
        core.print_out("  Computed a total of {} Large products\n\n".format(stats[-1]['product_count']))


def _gs_orth(engine, old_vecs, new_vec, thresh):
    """Perform GS orthonormalization of new_vec against some set old_vecs
    Parameters
    ----------
    engine : object
       The engine passed to the solver, required to define vector algebraic operations needed
    old_vecs : list of `vector`
       A set of orthonormal vectors, can be empty so that a set of vectors can be build via successive calls
    new_vec : single `vector`
       The vector to orthonormalize against old_vecs
    thresh : float
       If the orthogonalized vector has a norm smaller than this value it is considered LD to the set

    Returns
    -------
    new : single `vector` (if norm > thresh)
       new_vec after all components from old_vec have been projected out and the result normalized
    None: if (norm < thresh)
       When the final vector has norm smaller than thresh and should be discarded, none is returned.

    ..note:: Since discarded vectors are returned as None, when appending to a set with successive calls to the this function
    the caller should check the return before appending.
    """
    for old in old_vecs:
        dot = engine.vector_dot(old, new_vec)
        new_vec = engine.vector_axpy(-1.0 * dot, old, new_vec)
    norm = engine.vector_dot(new_vec, new_vec)
    norm = np.sqrt(norm)
    if norm >= thresh:
        new_vec = engine.vector_scale(1.0 / norm, new_vec)
        return new_vec
    else:
        return None


def _best_vectors(engine, ss_vectors, basis_vectors):
    """Best approximation of the true eigenvectors:

    ..math:: V^{k} = \Sum_{i} \tilde{V}^{k}_{i}X_{i}

    Where :math:`\tilde{V}^{k} is the `k`th eigenvector of the subspace problem and :math:`X_{i}` is an expansion vector

    Parameters
    ----------
    engine : object
       The engine passed to the solver, required to define vector algebraic operations needed
    ss_vectors : :py:class:`np.ndarray` {l, k}
       The k eigenvectors of the subspace problem, l = dimension of the subspace basis, and k is the number of roots
    basis_vectors : list of `vector` {l}
       The current basis vectors

    Returns
    -------
    new_vecs : list of `vector` {k}
       The approximations of the k true eigenvectors.
    """
    l, n = ss_vectors.shape
    new_vecs = []
    for i in range(n):
        cv_i = engine.new_vector()
        for j in range(l):
            cv_i = engine.vector_axpy(ss_vectors[j, i], basis_vectors[j], cv_i)
        new_vecs.append(cv_i)
    return new_vecs


def davidson_solver(engine,
                    guess,
                    e_tol=1.0E-6,
                    r_tol=1.0E-8,
                    nroot=1,
                    max_vecs_per_root=20,
                    maxiter=100,
                    verbose=1,
                    schmidt_add_tol=1.0e-4):
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
    schmidt_add_tol : float
        Correction vectors must have norm larger than this value to be added to the guess space
    verbose : int
        The amount of logging info to print (0 -> none, 1 -> some, >1 -> everything)

    Requirements of engine
    ----------------------
    compute_products(X) -> AX :
       X : list of `vectors`
       AX : list of `vectors`
           The product :math:`A x X_{i}` for each `X_{i}` in `X`, in that order. Where a is the hermitian matrix to be diagonalized.

    precondition(R_k, w_k) -> new_X_k
       R_k : single `vector`, the residual vector
       w_k : float, the eigenvalue associated with this vector
       new_X_k : single `vector`
           This vector will be added to the guess space after orthogonalization

    vector_dot(X, Y) -> a:
       X : single `vector`
       Y : single `vector`
       a : float, the dot product  (X x Y)

    vector_axpy(a, X, Y) - > Y':
       a : float
       X : single `vector`
       Y : single `vector`
       Y': single `vector`
         Y' = Y + a*X (it is assumed that the vector passed as Y is modified. This may not be true)

    vector_scale(a, X) -- > X':
       a : float
       X : single `vector`
       X': single `vector`
          X' = a * X (it is assumed that the vector passed as X is modified. This may not be true)

    vector_copy(X) -- > X':
       X : single `vector`
       X': single `vector`
         a copy of X

     new_vector() -- > X:
       X : single `vector`
         A new object that has the same shape as `vector` like quantities.


     ..note:: The `vector` referred to here is intentionally vague, the solver does not care what it is and only
              holds individual or sets of them. In fact an individual `vector` could be split across two elements in a list,
              such as for different spin. Whatever data type is used and individual vector should be a single element in a list such that
              len(list) returns the number of vector-like objects.
    """
    nk = nroot

    iter_info = {
        "count": 0,
        "res_norm": np.zeros((nk)),
        "val": np.zeros((nk)),
        "delta_val": np.zeros((nk)),
        # conv defaults to true, and will be flipped when a non-conv root is hit
        "done": True,
        "nvec": None,
        "collapse": False,
        "product_count": 0,
    }

    print_name = "DavidsonSolver"
    title_lines = ["Generalized Davidson Solver", "By Ruhee Dcunha"]
    max_ss_size = max_vecs_per_root * nk

    _diag_print_heading(title_lines, print_name, max_ss_size, nroot, r_tol, e_tol, maxiter, verbose)

    # NOTE: ignore passed guess always generate max_ss / vectors and force collapse on first iter
    vecs = engine.generate_guess(max_ss_size)
    # nguess_v_passed = len(guess)
    # vecs = []
    # # make guess set orthonormal
    # for v in guess:
    #     new = _gs_orth(engine, vecs, v, schmidt_add_tol)
    #     if new is not None:
    #         vecs.append(new)

    # nguess_v_after_orth = len(vecs)
    # # print warning if LD in passed guesses had to be removed.
    # if (nguess_v_after_orth < nguess_v_passed) and (verbose > 0):
    #     ndiscard = nguess_v_passed - nguess_v_after_orth
    #     core.print_out(
    #         "\n***Warning: Linear dependencies detected in initial guess vector, {} passed vectors discarded\n".format(
    #             ndiscard))

    # # raise exception if we don't have at least nroot guesses
    # if nguess_v_after_orth < nk:
    #     raise Exception(
    #         "At least nroot ({}) ortho-normal guess vectors must be provided. After orthonormalization only {} vectors remain".
    #         format(nk, nguess_v_after_orth))

    stats = []

    while iter_info['count'] < maxiter:
        # increment iteration/ save old vals
        iter_info['count'] += 1
        old_vals = iter_info['val'].copy()
        # reset flags
        iter_info['collapse'] = False
        iter_info['done'] = True

        l = len(vecs)
        iter_info['nvec'] = l
        if l >= max_ss_size:
            iter_info['collapse'] = True

        Ax = engine.compute_products(vecs)
        G = np.zeros((l, l))
        for i in range(l):
            for j in range(l):
                G[i, j] = engine.vector_dot(vecs[i], Ax[j])

        lam, alpha = np.linalg.eigh(G)
        idx = np.argsort(lam)
        lam = lam[idx]
        alpha = alpha[:, idx]

        # update best_solutions, failure return
        best_eigvecs = _best_vectors(engine, alpha[:, :nk], vecs)

        for k in range(nk):
            Rk = engine.new_vector()
            lam_k = lam[k]
            for i in range(l):
                Axi = Ax[i]
                Rk = engine.vector_axpy(alpha[i, k], Axi, Rk)

            Rk = engine.vector_axpy(-1.0 * lam_k, best_eigvecs[k], Rk)
            norm = engine.vector_dot(Rk, Rk)
            norm = np.sqrt(norm)

            iter_info['val'][k] = lam_k
            iter_info['delta_val'][k] = abs(old_vals[k] - lam_k)
            iter_info['res_norm'][k] = norm
            converged = (norm < r_tol) and (abs(old_vals[k] - lam_k) < e_tol)
            # We want to expand trial basis when, solution is not converged or when trial space is too small.
            if (not converged):
                Qk = engine.precondition(Rk, lam_k)
                iter_info['done'] = False
                Qk = engine.vector_scale(1.0 / norm, Rk)
                Qk = _gs_orth(engine, vecs, Qk, schmidt_add_tol)
                if Qk is not None:
                    vecs.append(Qk)

        _diag_print_info(print_name, iter_info, verbose)
        stats.append(iter_info.copy())
        if iter_info['done']:
            lam = lam[:nk]
            _diag_print_converged(print_name, stats, lam, rvec=best_eigvecs, verbose=verbose)
            return lam, best_eigvecs, stats
        if iter_info['collapse']:
            vecs = best_eigvecs
            best_eigvecs = []

    # If we get down here  we have exceeded max iterations without convergence, return None + stats, (let caller raise the error)
    return None, None, stats


def hamiltonian_solver(engine,
                       guess,
                       e_tol=1.0E-6,
                       r_tol=1.0E-8,
                       nroot=1,
                       max_vecs_per_root=20,
                       maxiter=100,
                       verbose=1,
                       schmidt_add_tol=1.0e-8):
    """
    Finds the smallest eigenvalues and associated right and left hand eigenvectors of a large real Hamiltonian eigenvalue problem
    emulated through an engine.

    A hamiltonian EVP has the structure with A, B of some large dimension N the problem is 2Nx2N:
    [A  B][X]  = [1   0](w)[X]
    [B  A][Y]    [0  -1](w)[Y]

    Which can be written as the NxN, non-hermitan EVP:
    (A+B)(A-B)(X+Y) = w^2(X+Y)


    With left-hand eigenvectors:
    (X-Y)(A-B)(A+B) = w^2(X-Y)

    if (A-B) is positive definite, we can transform the problem to arrive at the hermitian NxN EVP:
    (A-B)^1/2(A+B)(A-B)^1/2 = w^2 T

    Where T = (A-B)^-1/2(X+Y).

    We use a Davidson like where we transform (A+B) (H1) and (A-B)(H2) in to the subspace defined by the trial vectors. Transform the subspace quantities and solve the
    subspace analog of the NxN hermitian EVP, and back transform to extract X+Y, X-Y. The subspace is augmented to correct for both until convergence.

    Using the converged Left/Right eigenvectors of the NxN non-hermitian EVP, the components of the 2n*2N Hamiltonian problem (X, Y) are extracted.

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
    schmidt_add_tol : float
        Correction vectors must have norm larger than this value to be added to the guess space
    verbose : int
        The amount of logging info to print (0 -> none, 1 -> some, >1 -> everything)


    Requirements of engine
    ----------------------
    compute_products(X) -> H1X, H2X :
       X : list of `vectors`
       H1X : list of `vectors`
           The product :math:`H1 x X_{i}` for each `X_{i}` in `X`, in that order. Where H1 is the sum (A+B) of the blocks described above.
       H2X : list of `vectors`
           The product :math:H2 x X_{i} for each `X_{i}` in  X, in that order. Where H2 is the difference (A-B) of the blocks described above.

    precondition(R_k, w_k) -> new_X_k
       R_k : single `vector`, the residual vector
       w_k : float, the eigenvalue associated with this vector
       new_X_k : single `vector`
           This vector will be added to the guess space after orthogonalization

    vector_dot(X, Y) -> a:
       X : single `vector`
       Y : single `vector`
       a : float, the dot product  (X x Y)

    vector_axpy(a, X, Y) - > Y':
       a : float
       X : single `vector`
       Y : single `vector`
       Y': single `vector`
         Y' = Y + a*X (it is assumed that the vector passed as Y is modified. This may not be true)

    vector_scale(a, X) -- > X':
       a : float
       X : single `vector`
       X': single `vector`
          X' = a * X (it is assumed that the vector passed as X is modified. This may not be true)

    vector_copy(X) -- > X':
       X : single `vector`
       X': single `vector`
         a copy of X

     new_vector() -- > X:
       X : single `vector`
         A new object that has the same shape as `vector` like quantities.


     ..note:: The `vector` referred to here is intentionally vague, the solver does not care what it is and only
              holds individual or sets of them. In fact an individual `vector` could be split across two elements in a list,
              such as for different spin. Whatever data type is used and individual vector should be a single element in a list such that
              len(list) returns the number of vector-like objects.
    """

    nk = nroot

    iter_info = {
        "count": 0,
        "res_norm": np.zeros((nk)),
        "val": np.zeros((nk)),
        "delta_val": np.zeros((nk)),
        # conv defaults to true, and will be flipped when a non-conv root is hit
        "conv": True,
        "nvec": 0,
        "product_count": 0,
    }
    print_name = "HamiltonianSolver"
    title_lines = ["Generalized Hamiltonian Solver", "By Andrew M. James"]
    ss_max = max_vecs_per_root * nk

    _diag_print_heading(title_lines, print_name, ss_max, nroot, r_tol, e_tol, maxiter, verbose)

    # NOTE: Ignore passed guess, always generate 1/2 * max ss size vectors and allow collapse on first iter
    vecs = engine.generate_guess(ss_max)

    # NOTE: Re-enable guess ortho-normalization if user guess allowed.
    # nguess_v_passed = len(guess)
    # vecs = []
    # # make guess set orthonormal
    # for v in guess:
    #     new = _gs_orth(engine, vecs, v, schmidt_add_tol)
    #     if new is not None:
    #         vecs.append(new)

    # nguess_v_after_orth = len(vecs)
    # # print warning if LD in passed guesses had to be removed.
    # if (nguess_v_after_orth < nguess_v_passed) and (verbose > 0):
    #     ndiscard = nguess_v_passed - nguess_v_after_orth
    #     core.print_out(
    #         "\n***Warning: Linear dependencies detected in initial guess vector, {} passed vectors discarded\n".format(
    #             ndiscard))

    # raise exception if we don't have at least nroot guesses
    # if nguess_v_after_orth < nk:
    #     raise Exception(
    #         "At least nroot ({}) ortho-normal guess vectors must be provided. After orthonormalization only {} vectors remain".
    #         format(nk, nguess_v_after_orth))

    stats = []
    while iter_info['count'] < maxiter:
        # increment iteration/ save old vals
        iter_info['count'] += 1
        old_w = iter_info['val'].copy()
        # reset flags
        iter_info['collapse'] = False
        iter_info['done'] = True

        # expand subspace
        l = len(vecs)
        if l >= ss_max:
            iter_info['collapse'] = True

        iter_info['nvec'] = l
        # Step 2: compute (P)*bi and (M)*bi
        H1x, H2x = engine.compute_products(vecs)
        # Step 3: form Pss, Mss. The P/M matrices in the subspace
        H1_ss = np.zeros((l, l))
        H2_ss = np.zeros((l, l))
        for i in range(l):
            for j in range(l):
                H1_ss[i, j] = engine.vector_dot(vecs[i], H1x[j])
                H2_ss[i, j] = engine.vector_dot(vecs[i], H2x[j])

        # Step 4: Hermitian Product (Subspace analog of M^{1/2} P M^{1/2})
        H2_ss_val, H2_ss_vec = np.linalg.eigh(H2_ss)
        if any(H2_ss_val < 0.0):
            raise Exception("H2 is not Positive Definite")
        H2_ss_half = np.einsum('ij,j,kj->ik', H2_ss_vec, np.sqrt(H2_ss_val), H2_ss_vec)
        Hss = np.einsum('ij,jk,km->im', H2_ss_half, H1_ss, H2_ss_half)

        # Step 5: diagonalize Hss -> w^2, Tss
        w, Tss = np.linalg.eigh(Hss)
        w = np.sqrt(w)
        idx = np.argsort(w)
        w = w[idx]
        Tss = Tss[:, idx]

        #Step 6a: extract Rss = H2_ss^{1/2}Tss
        Rss = np.dot(H2_ss_half, Tss)
        #Step 6b: extract Lss = Pss Rss * w^-1
        Lss = np.dot(H1_ss, Rss)
        Lss = np.einsum('ij,j->ij', Lss, np.divide(1.0, w))

        # #step 6c: R/L in the subspace are already biorthogonal, we need to normalize
        for k in range(l):
            dot = np.sqrt(Rss[:, k].dot(Lss[:, k]))
            Rss[:, k] /= dot
            Lss[:, k] /= dot

        best_R = _best_vectors(engine, Rss[:, :nk], vecs)
        best_L = _best_vectors(engine, Lss[:, :nk], vecs)

        # compute residuals
        for k in range(nk):
            WR_k = engine.new_vector()
            WL_k = engine.new_vector()
            wk = w[k]
            for i in range(l):
                H1x_i = H1x[i]
                H2x_i = H2x[i]
                WL_k = engine.vector_axpy(Rss[i, k], H1x_i, WL_k)
                WR_k = engine.vector_axpy(Lss[i, k], H2x_i, WR_k)

            WL_k = engine.vector_axpy(-1.0 * wk, best_L[k], WL_k)
            WR_k = engine.vector_axpy(-1.0 * wk, best_R[k], WR_k)

            norm_R = engine.vector_dot(WR_k, WR_k)
            norm_R = np.sqrt(norm_R)

            norm_L = engine.vector_dot(WL_k, WL_k)
            norm_L = np.sqrt(norm_L)

            norm = norm_R + norm_L

            iter_info['res_norm'][k] = norm
            iter_info['delta_val'][k] = np.abs(old_w[k] - w[k])
            iter_info['val'][k] = w[k]
            if (iter_info['res_norm'][k] > r_tol) or (iter_info['delta_val'][k] > e_tol):
                iter_info['done'] = False

                QR_k = engine.precondition(WR_k, w[k])
                QR_k = _gs_orth(engine, vecs, QR_k, schmidt_add_tol)
                if QR_k is not None:
                    vecs.append(QR_k)

                QL_k = engine.precondition(WL_k, w[k])
                QL_k = _gs_orth(engine, vecs, QL_k, schmidt_add_tol)
                if QL_k is not None:
                    vecs.append(QL_k)

        _diag_print_info(print_name, iter_info, verbose)
        stats.append(iter_info.copy())
        if iter_info['done']:
            _diag_print_converged(print_name, stats, w[:nk], rvec=best_R, lvec=best_L, verbose=verbose)
            return w[:nk], best_R, best_L, stats
        # number of vectors hasn't changed (nothing new added), but we are not converged. Force collapse
        if len(vecs) == iter_info['nvec']:
            iter_info['collapse'] = True
        if iter_info['collapse']:
            vecs = []
            for r in best_R:
                new = _gs_orth(engine, vecs, r, schmidt_add_tol)
                if new is not None:
                    vecs.append(new)
            for l in best_L:
                new = _gs_orth(engine, vecs, l, schmidt_add_tol)
                if new is not None:
                    vecs.append(new)

    # if we get down here we have exceeded max iterations without convergence, return none + stats
    return None, None, None, stats
