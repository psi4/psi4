#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2022 The Psi4 Developers.
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

import numpy as np

from psi4 import core
from psi4.driver.p4util.exceptions import *
"""
This module provides ``engine`` objects that can be used by the :func:`~psi4.driver.p4util.solvers.davidson_solver` and
:func:`~psi4.driver.p4util.solvers.hamiltonian_solver`

Spin Orbital Expressions for the Relevant product components
------------------------------------------------------------
Aia,jb  = (E_a - E_i) +     J    -    K
        = (E_a - E_i) +  (ia|jb) - (ij|ab)

Bia,jb  =     J    -   K^T
        =  (ai|bj) - (aj|bi)

H2ia,jb =                A                   -           B
        = [(E_a - E_i) +     J    -    K   ]  - [   J     -    K^T ]
        = [(E_a - E_i) +  (ia|jb) - (ij|ab)]  - [  (ai|bj)  - (aj|bi)]
        = (E_a - E_i) -    K    +    K^T
        = (E_a - E_i) - (ij|ab) + (aj|bi)

H1ia,jb =              A                    +       B
        = [(E_a - E_i) +     J    -    K   ] + [  J        -   K^T  ]
        = [(E_a - E_i) +  (ia|jb) - (ij|ab)] + [  (ai|bj)  - (aj|bi)]
        = [(E_a - E_i) +     J    -    K     -   K^T]
        = [(E_a - E_i) +  (ia|jb) - (ij|ab)  - (aj|bi)]
"""


class SingleMatPerVector:
    """Operations for RHF-like systems where the `vector` is a single :py:class:`psi4.core.Matrix`
    """

    @staticmethod
    def vector_dot(X, Y):
        return X.vector_dot(Y)

    @staticmethod
    def vector_scale(a, X):
        X.scale(a)
        return X

    @staticmethod
    def vector_axpy(a, X, Y):
        Y.axpy(a, X)
        return Y

    @staticmethod
    def vector_copy(X):
        return X.clone()

    @staticmethod
    def vector_transpose(X):
        return X.transpose()


class PairedMatPerVector:
    """Operations for UHF-like systems where the vector is a pair of :py:class:`psi4.core.Matrix` objects holding
    (Alpha, Beta) components.
    """

    @staticmethod
    def vector_dot(X, Y):
        dot = X[0].vector_dot(Y[0])
        dot += X[1].vector_dot(Y[1])
        return dot

    @staticmethod
    def vector_scale(a, X):
        X[0].scale(a)
        X[1].scale(a)
        return X

    @staticmethod
    def vector_axpy(a, X, Y):
        Y[0].axpy(a, X[0])
        Y[1].axpy(a, X[1])
        return Y

    @staticmethod
    def vector_copy(X):
        return [X[0].clone(), X[1].clone()]

    @staticmethod
    def vector_transpose(X):
        return [X[0].transpose(), X[1].transpose()]


class ProductCache:
    """Caches product vectors
    """

    def __init__(self, *product_types):
        """Creates a new Product Cache

        Parameters
        ----------
        *product_types, list of str
            A list of product labels
        """
        self._products = {p: [] for p in product_types}

    def add(self, pkey, new_elements):
        """Adds new elements to a given key

        Parameters
        ----------
        pkey : str
            Product label
        new_elements : list of arrays
            New products to add to the cache

        """
        if pkey not in self._products.keys():
            raise AttributeError("No such product {}".format(pkey))

        for new in new_elements:
            self._products[pkey].append(new)

        return self._products[pkey].copy()

    def reset(self):
        """Resets the ProductCache by clearing all data.
        """
        for pkey in self._products.keys():
            self._products[pkey].clear()

    def count(self):
        """Return the number of cached products
        """
        lens = [len(self._products[pkey]) for pkey in self._products.keys()]
        all_same = all(lens[0] == x for x in lens)

        if all_same:
            return lens[0]
        else:
            raise ValueError("Cache lengths are not the same, invalid cache error. Call a developer.")


class TDRSCFEngine(SingleMatPerVector):
    """Engine for R(HF/KS) products

    Fulfills the API required by :class:`~psi4.driver.p4util.solvers.SolverEngine`


    Parameters
    ----------
    wfn : :py:class:`psi4.core.Wavefunction`
        The converged SCF wfn
    ptype : {'rpa', 'tda'}
        The product type to be evaluated. When ``ptype == 'rpa'``. The return of `compute_products` will be as
        expected by :func:`~psi4.driver.p4util.solvers.hamiltonian_solver`, when ``ptype == 'tda'`` the return of
        compute_products will be as expected by :func:`~psi4.driver.p4util.solvers.davidson_solver`.
    triplet : bool , optional
        Are products spin-adapted for triplet excitations?
    """

    def __init__(self, wfn, *, ptype, triplet=False):

        # primary data
        self.wfn = wfn
        self.ptype = ptype.lower()
        self.needs_K_like = self.wfn.functional().is_x_hybrid() or self.wfn.functional().is_x_lrc()

        if self.ptype not in ["rpa", "tda"]:
            raise KeyError(f"Product type {self.ptype} not understood")

        # product type
        self.singlet = not triplet
        if ptype == 'rpa':
            self.product_cache = ProductCache("H1", "H2")
        else:
            self.product_cache = ProductCache("A")

        # orbitals and eigenvalues
        self.Co = self.wfn.Ca_subset("SO", "OCC")
        self.Cv = self.wfn.Ca_subset("SO", "VIR")
        self.E_occ = self.wfn.epsilon_a_subset("SO", "OCC")
        self.E_vir = self.wfn.epsilon_a_subset("SO", "VIR")
        self.prec = None

        # ground state symmetry
        self.G_gs = 0

        # ground state spin multiplicity
        self.mult_gs = wfn.molecule().multiplicity()

        # excited state symmetry
        self.G_es = None

        # symmetry of transition
        self.occpi = self.wfn.nalphapi()
        self.virpi = self.wfn.nmopi() - self.occpi
        self.nsopi = self.wfn.nsopi()
        self.reset_for_state_symm(0)

    ## API required by "engine" (see p4util.solvers.davidson/hamiltonian_solver)

    def new_vector(self, name=""):
        """Obtain a blank matrix object with the correct symmetry"""
        return core.Matrix(name, self.occpi, self.virpi, self.G_trans)

    def reset_for_state_symm(self, symmetry):
        """Reset internal quantities so the object is prepared to deal with transition to state with symmetry given
        """
        self.G_es = symmetry
        self._build_prec()
        self.product_cache.reset()

    def compute_products(self, vectors):
        """Given a set of vectors X Compute products
        if ptype == rpa:
           Returns pair (A+B)X, (A-B)X
        if ptype == tda:
           Returns AX
        """

        n_old = self.product_cache.count()
        n_new = len(vectors)

        if n_new <= n_old:
            self.product_cache.reset()
            compute_vectors = vectors
        else:
            compute_vectors = vectors[n_old:]

        n_prod = len(compute_vectors)

        # Build base one and two electron quantities
        Fx = self.wfn.onel_Hx(compute_vectors)
        twoel = self.wfn.twoel_Hx(compute_vectors, False, "SO")
        Jx, Kx = self._split_twoel(twoel)

        # Switch between rpa and tda
        if self.ptype == 'rpa':
            H1X_new, H2X_new = self._combine_H1_H2(Fx, Jx, Kx)
            for H1x in H1X_new:
                self.vector_scale(-1.0, H1x)
            for H2x in H2X_new:
                self.vector_scale(-1.0, H2x)

            H1X_all = self.product_cache.add("H1", H1X_new)
            H2X_all = self.product_cache.add("H2", H2X_new)
            return H1X_all, H2X_all, n_prod

        else:
            AX_new = self._combine_A(Fx, Jx, Kx)
            for Ax in AX_new:
                self.vector_scale(-1.0, Ax)
            AX_all = self.product_cache.add("A", AX_new)
            return AX_all, n_prod

    def precondition(self, Rvec, shift):
        """Applies the preconditioner with a shift to a residual vector

        value = R / (shift - preconditioner)
        """
        for h in range(self.wfn.nirrep()):
            den = shift - self.prec.nph[h]
            den[np.abs(den) < 0.0001] = 1.0
            Rvec.nph[h][:] /= den
        return Rvec

    def generate_guess(self, nguess):
        """Generate a set of guess vectors based on orbital energy differences
        """
        deltas = []
        guess_vectors = []
        for ho in range(self.wfn.nirrep()):
            hv = ho ^ self.G_trans
            for i, ei in enumerate(self.E_occ.nph[ho]):
                for a, ea in enumerate(self.E_vir.nph[hv]):
                    deltas.append((ea - ei, i, a, ho))
        deltas_sorted = sorted(deltas, key=lambda x: x[0])
        nguess = min(nguess, len(deltas_sorted))
        for i in range(nguess):
            v = self.new_vector()
            oidx = deltas_sorted[i][1]
            vidx = deltas_sorted[i][2]
            h = deltas_sorted[i][3]
            v.set(h, oidx, vidx, 1.0)
            guess_vectors.append(v)
        return guess_vectors

    def residue(self, X, so_prop_ints):
        # return zeros if spin multiplicity of GS and ES differ
        if not self.singlet and (self.mult_gs == 1):
            return np.zeros(len(so_prop_ints))

        prop = [core.triplet(self.Co, x, self.Cv, True, False, False) for x in so_prop_ints]

        return np.sqrt(2.0) * np.array([X.vector_dot(u) for u in prop])

    ## Helper functions

    def _combine_H1_H2(self, Fx, Jx, Kx=None):
        """Build the combinations:
        Singlet:
            H1 X =  [(Ea - Ei) + 4J - K - K^T]X
            H2 X =  [(Ea - Ei) - K  + K^T]X

        Triplet:
            H1 X =  [(Ea - Ei) - K - K^T]X
            H2 X =  [(Ea - Ei) - K + K^T]X
        """

        H1X = []
        H2X = []
        if Kx is not None:
            for Fxi, Jxi, Kxi in zip(Fx, Jx, Kx):
                Kxit = self.vector_transpose(Kxi)
                # H1x = -K singlet/triplet
                H1X_so = self.vector_copy(Kxi)
                H1X_so = self.vector_scale(-1.0, H1X_so)
                # H1X -= K^T singlet/triplet
                H1X_so = self.vector_axpy(-1.0, Kxit, H1X_so)
                # H2x = K^T  - K singlet/triplet
                H2X_so = self.vector_axpy(-1.0, Kxi, Kxit)
                if self.singlet:
                    # H1x += 4*J (singlet only)
                    H1X_so = self.vector_axpy(4.0, Jxi, H1X_so)

                # transform + add Ea-Ei
                H1X.append(self.vector_axpy(1.0, Fxi, self._so_to_mo(H1X_so)))
                H2X.append(self.vector_axpy(1.0, Fxi, self._so_to_mo(H2X_so)))

        else:
            for Fxi, Jxi in zip(Fx, Jx):
                if self.singlet:
                    H1X_so = self.vector_scale(4.0, Jxi)
                    H1X.append(self.vector_axpy(1.0, Fxi, self._so_to_mo(H1X_so)))
                else:
                    H1X.append(self.vector_copy(Fxi))
                H2X.append(self.vector_copy(Fxi))
        return H1X, H2X

    def _combine_A(self, Fx, Jx, Kx=None):
        """Build the combinations
        Singlet:
           A X = [(Ea - Ei) + 2 J - K] X

        Triplet:
           A X = [(Ea - Ei) - K] X
        """
        Ax = []
        if Kx is not None:
            for Fxi, Jxi, Kxi in zip(Fx, Jx, Kx):
                Ax_so = self.vector_scale(-1.0, self.vector_copy(Kxi))
                if self.singlet:
                    Ax_so = self.vector_axpy(2.0, Jxi, Ax_so)
                Ax.append(self.vector_axpy(1.0, Fxi, self._so_to_mo(Ax_so)))

        else:
            for Fxi, Jxi in zip(Fx, Jx):
                if self.singlet:
                    Ax.append(self.vector_axpy(1.0, Fxi, self._so_to_mo(self.vector_scale(2.0, Jxi))))
                else:
                    Ax.append(self.vector_copy(Fxi))
        return Ax

    def _so_to_mo(self, X):
        """Transform (C_occ)^T X C_vir"""
        return core.triplet(self.Co, X, self.Cv, True, False, False)

    def _split_twoel(self, twoel):
        """Unpack J and K matrices
        """
        if self.needs_K_like:
            Jx = twoel[0::2]
            Kx = twoel[1::2]
        else:
            Jx = twoel
            Kx = None
        return Jx, Kx

    @property
    def G_trans(self):
        """The symmetry of the transition vector"""
        return self.G_gs ^ self.G_es

    def _build_prec(self):
        """Builds energy denominator
        """
        self.prec = self.new_vector()
        for h in range(self.wfn.nirrep()):
            self.prec.nph[h][:] = self.E_vir.nph[h ^ self.G_trans] - self.E_occ.nph[h].reshape(-1, 1)


class TDUSCFEngine(PairedMatPerVector):
    """Engine for U(HF/KS) products

    Fulfills the API required by :class:`~psi4.driver.p4util.solvers.SolverEngine`

    Parameters
    ----------
    wfn : :py:class:`psi4.core.Wavefunction`
        The converged SCF wfn
    ptype : str {'rpa', 'tda'}
        The product type to be evaluated. When ``ptype == 'rpa'``. The return of `compute_products` will be as
        expected by :func:`~psi4.driver.p4util.solvers.hamiltonian_solver`, when ``ptype == 'tda'`` the return of
        compute_products will be as expected by :func:`~psi4.driver.p4util.solvers.davidson_solver`.

    """

    def __init__(self, wfn, *, ptype):

        # Primary data
        self.wfn = wfn
        self.ptype = ptype

        # Find product type
        if ptype == 'rpa':
            self.product_cache = ProductCache("H1", "H2")
        elif ptype == 'hess':
            self.product_cache = ProductCache("H1")
        else:
            self.product_cache = ProductCache("A")

        # Save orbitals and eigenvalues
        self.Co = [wfn.Ca_subset("SO", "OCC"), wfn.Cb_subset("SO", "OCC")]
        self.Cv = [wfn.Ca_subset("SO", "VIR"), wfn.Cb_subset("SO", "VIR")]
        self.E_occ = [wfn.epsilon_a_subset("SO", "OCC"), wfn.epsilon_b_subset("SO", "OCC")]
        self.E_vir = [wfn.epsilon_a_subset("SO", "VIR"), wfn.epsilon_b_subset("SO", "VIR")]
        self.needs_K_like = self.wfn.functional().is_x_hybrid() or self.wfn.functional().is_x_lrc()

        # dimensions
        self.occpi = [self.wfn.nalphapi(), self.wfn.nbetapi()]
        self.virpi = [self.wfn.nmopi() - self.occpi[0], self.wfn.nmopi() - self.occpi[1]]
        self.nsopi = self.wfn.nsopi()

        # Orbital energy differences
        self.prec = [None, None]

        # Ground state symmetry
        self.G_gs = 0
        for h in range(self.wfn.nirrep()):
            for i in range(self.occpi[0][h] - self.occpi[1][h]):
                self.G_gs = self.G_gs ^ h
        # Excited state symmetry
        self.G_es = None
        self.reset_for_state_symm(0)

    ## API Required by "engine" (see p4util.solvers.davidson/hamiltonian_solver)

    def precondition(self, Rvec, shift):
        """Applies the preconditioner with a shift to a residual vector

        value = R / (shift - preconditioner)
        """
        for h in range(self.wfn.nirrep()):

            den = shift - self.prec[0].nph[h]
            den[abs(den) < 0.0001] = 1.0
            Rvec[0].nph[h][:] /= den

            den = shift - self.prec[1].nph[h]
            den[abs(den) < 0.0001] = 1.0
            Rvec[1].nph[h][:] /= den
        return Rvec

    def compute_products(self, vectors):
        """Compute Products for a list of guess vectors (X).
        if ptype == 'rpa':
                           H1 ,   H2
           returns pair (A+B)X, (A-B)X products
        if ptype == 'hess':
           returns (A+B)X products
        if ptype == 'tda':
           returns Ax products.
        """

        n_old = self.product_cache.count()
        n_new = len(vectors)
        if n_new <= n_old:
            self.product_cache.reset()
            compute_vectors = vectors
        else:
            compute_vectors = vectors[n_old:]

        n_prod = len(compute_vectors)

        # flatten list of [(A,B)_i, ...] to [A_i, B_i, ...]
        vec_flat = sum(compute_vectors, [])

        Fx = self._pair_onel(self.wfn.onel_Hx(vec_flat))
        twoel = self.wfn.twoel_Hx(vec_flat, False, "SO")
        Jx, Kx = self._split_twoel(twoel)

        if self.ptype == "rpa":
            H1X_new, H2X_new = self._combine_H1_H2(Fx, Jx, Kx)
            for H1x in H1X_new:
                self.vector_scale(-1.0, H1x)
            for H2x in H2X_new:
                self.vector_scale(-1.0, H2x)
            H1X_all = self.product_cache.add("H1", H1X_new)
            H2X_all = self.product_cache.add("H2", H2X_new)
            return H1X_all, H2X_all, n_prod
        elif self.ptype == "hess":
            H1X_new = self._combine_H1(Fx, Jx, Kx)
            for H1x in H1X_new:
                self.vector_scale(-1.0, H1x)
            H1X_all = self.product_cache.add("H1", H1X_new)
            return H1X_all, n_prod
        else:
            AX_new = self._combine_A(Fx, Jx, Kx)
            for Ax in AX_new:
                self.vector_scale(-1.0, Ax)
            AX_all = self.product_cache.add("A", AX_new)
            return AX_all, n_prod

    def generate_guess(self, nguess):
        """Generate a set of guess vectors based on orbital energy differences
        """
        guess_vectors = []
        deltas = []
        for ho in range(self.wfn.nirrep()):
            hv = self.G_trans ^ ho
            for i, ei in enumerate(self.E_occ[0].nph[ho]):
                for a, ea in enumerate(self.E_vir[0].nph[hv]):
                    deltas.append((ea - ei, 0, i, a, ho))
            for i, ei in enumerate(self.E_occ[1].nph[ho]):
                for a, ea in enumerate(self.E_vir[1].nph[hv]):
                    deltas.append((ea - ei, 1, i, a, ho))

        deltas_sorted = sorted(deltas, key=lambda x: x[0])
        nguess = min(nguess, len(deltas_sorted))
        for i in range(nguess):
            v = self.new_vector()
            spin = deltas_sorted[i][1]
            oidx = deltas_sorted[i][2]
            vidx = deltas_sorted[i][3]
            h = deltas_sorted[i][4]
            v[spin].set(h, oidx, vidx, 1.0)
            guess_vectors.append(v)
        return guess_vectors

    def new_vector(self, name=""):
        """Build a new object with shape symmetry like a trial vector """
        return [
            core.Matrix(name + 'a', self.occpi[0], self.virpi[0], self.G_trans),
            core.Matrix(name + 'b', self.occpi[1], self.virpi[1], self.G_trans)
        ]

    def reset_for_state_symm(self, symmetry):
        """Reset internal quantities so the object is prepared to deal with transition to state with symmetry given
        """
        self.G_es = symmetry
        self._build_prec()
        self.product_cache.reset()

    def residue(self, X, so_prop_ints):
        prop_a = [core.triplet(self.Co[0], x, self.Cv[0], True, False, False) for x in so_prop_ints]
        prop_b = [core.triplet(self.Co[1], x, self.Cv[1], True, False, False) for x in so_prop_ints]

        return np.array([X[0].vector_dot(u[0]) + X[1].vector_dot(u[1]) for u in zip(prop_a, prop_b)])

    ## Helper Functions

    def _combine_H1(self, Fx, Jx, Kx=None):
        """Build the combination:
            H1 X =  [(Ea - Ei) + 2J - K - K^T]X
        """

        H1X = []
        if Kx is not None:
            for Fxi, Jxi, Kxi in zip(Fx, Jx, Kx):
                H1X_so = self.vector_scale(2.0, Jxi)
                Kxit = self.vector_transpose(Kxi)
                H1X_so = self.vector_axpy(-1.0, Kxi, H1X_so)
                H1X_so = self.vector_axpy(-1.0, Kxit, H1X_so)
                H1X.append(self.vector_axpy(1.0, Fxi, self._so_to_mo(H1X_so)))
        else:
            for Fxi, Jxi in zip(Fx, Jx):
                H1X_so = self.vector_scale(2.0, Jxi)
                H1X.append(self.vector_axpy(1.0, Fxi, self._so_to_mo(H1X_so)))

        return H1X
    def _combine_H1_H2(self, Fx, Jx, Kx=None):
        """Build the combinations:
            H1 X =  [(Ea - Ei) + 2J - K - K^T]X
            H2 X =  [(Ea - Ei) - K  + K^T]X
        """

        H1X = []
        H2X = []
        if Kx is not None:
            for Fxi, Jxi, Kxi in zip(Fx, Jx, Kx):
                H1X_so = self.vector_scale(2.0, Jxi)
                Kxit = self.vector_transpose(Kxi)
                H1X_so = self.vector_axpy(-1.0, Kxi, H1X_so)
                H1X_so = self.vector_axpy(-1.0, Kxit, H1X_so)
                H1X.append(self.vector_axpy(1.0, Fxi, self._so_to_mo(H1X_so)))

                H2X_so = self.vector_axpy(-1.0, Kxi, Kxit)
                H2X.append(self.vector_axpy(1.0, Fxi, self._so_to_mo(H2X_so)))
        else:
            for Fxi, Jxi in zip(Fx, Jx):
                H1X_so = self.vector_scale(2.0, Jxi)
                H1X.append(self.vector_axpy(1.0, Fxi, self._so_to_mo(H1X_so)))
                H2X.append(self.vector_copy(Fxi))

        return H1X, H2X

    def _combine_A(self, Fx, Jx, Kx):
        """Build the combination

        A X = [(Ea-Ei) + J  - K] X
        """
        Ax = []
        if Kx is not None:
            for Fxi, Jxi, Kxi in zip(Fx, Jx, Kx):
                Ax_so = self.vector_axpy(-1.0, Kxi, Jxi)
                Ax.append(self.vector_axpy(1.0, Fxi, self._so_to_mo(Ax_so)))
        else:
            for Fxi, Jxi in zip(Fx, Jx):
                Ax.append(self.vector_axpy(1.0, Fxi, self._so_to_mo(Jxi)))
        return Ax

    def _so_to_mo(self, X):
        """Transform (C_occ)^T X C_vir"""
        return [core.triplet(self.Co[i], X[i], self.Cv[i], True, False, False) for i in (0, 1)]

    def _pair_onel(self, onel):
        """Pair up A/B from onel_Hx return"""
        return list(zip(onel[0::2], onel[1::2]))

    def _split_twoel(self, twoel):
        """Unpack J and K matrices and pair alpha/beta
        """
        if self.needs_K_like:
            Jx = list(zip(twoel[0::4], twoel[1::4]))
            Kx = list(zip(twoel[2::4], twoel[3::4]))
        else:
            Jx = list(zip(twoel[0::2], twoel[1::2]))
            Kx = None
        return Jx, Kx

    @property
    def G_trans(self):
        """Symmetry of transition vector"""
        return self.G_gs ^ self.G_es

    def _build_prec(self):
        """Builds energy denominator
        """
        self.prec = self.new_vector()
        for h in range(self.wfn.nirrep()):
            self.prec[0].nph[h][:] = self.E_vir[0].nph[h ^ self.G_trans] - self.E_occ[0].nph[h].reshape(-1, 1)
            self.prec[1].nph[h][:] = self.E_vir[1].nph[h ^ self.G_trans] - self.E_occ[1].nph[h].reshape(-1, 1)
