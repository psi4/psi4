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

import numpy as np

from psi4 import core
from psi4.driver.p4util.exceptions import *


"""
- H1/ H2 docstrings to keep them straight
- RHF equations

Aia,jb  = (E_a - E_i) +     J    -    K
        = (E_a - E_i) + 4(ia|jb) - (ij|ab)

Bia,jb  =    J    -   K^T
        = (ai|bj) - (aj|bi)

H2ia,jb =                A                   -           B
        = [(E_a - E_i) +    J    -    K   ]  - [   J     -    K^T ]
        = [(E_a - E_i) + (ia|jb) - (ij|ab)]  - [(ai|bj)  - (aj|bi)]
        = (E_a - E_i) -    K    +    K^T
        = (E_a - E_i) - (ij|ab) + (aj|bi)

H1ia,jb =              A                    +       B
        = [(E_a - E_i) +     J    -    K   ] + [J    -   K^T]
        = [(E_a - E_i) +  (ia|jb) - (ij|ab)] + [(ai|bj)  - (aj|bi)]
        = [(E_a - E_i) + 4   J    -    K     -   K^T]
        = [(E_a - E_i) + 4(ia|jb) - (ij|ab)  - (aj|bi)]

"""


class SingleMatPerVector:
    """Matrix operations for RHF-like systems
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
    """Matrix operations for UHF-like systems
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
        """

        """
        lens = [len(self._products[pkey]) for pkey in self._products.keys()]
        all_same = all(lens[0] == x for x in lens)

        if all_same:
            return lens[0]
        else:
            raise ValueError("Cache lengths are not the same, invalid cache error. Call a developer.")


class TDRSCFEngine(SingleMatPerVector):
    def __init__(self, wfn, ptype, triplet=False):

        # primary data
        self.wfn = wfn
        self.ptype = ptype.lower()
        self.needs_K_like = self.wfn.functional().is_x_hybrid() or self.wfn.functional().is_x_lrc()

        if self.ptype not in ["rpa", "tda"]:
            raise KeyError("Product type {} not understood".format(self.ptype))

        # product type
        self.triplet = triplet
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

        # excited state symmetry
        self.G_es = None

        # symmetry of transition
        self.occpi = self.wfn.nalphapi()
        self.virpi = self.wfn.nmopi() - self.occpi
        self.nsopi = self.wfn.nsopi()
        self.reset_for_state_symm(0)


## API required by engine

## Helper functions

    def new_vector(self, name=""):
        """Obtain a blank matrix object with the correct symmetry"""
        return core.Matrix(name, self.occpi, self.virpi, self.G_trans)

    def reset_for_state_symm(self, symmetry):
        """Reset Wavefuncitn to a new symmetry
        """
        self.G_es = symmetry
        self._build_prec()
        self.product_cache.reset()

    def combine_H1_H2(self, Fx, Jx, Kx=None):
        """
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
                if self.triplet:
                    H1X_so = self.vector_scale(0.0, Jxi)
                else:
                    H1X_so = self.vector_scale(4.0, Jxi)
                H1X_so = self.vector_axpy(-1.0, Kxi, H1X_so)
                H1X_so = self.vector_axpy(-1.0, Kxit, H1X_so)
                H1X.append(self.vector_axpy(1.0, Fxi, self.so_to_mo(H1X_so)))

                H2X_so = self.vector_axpy(-1.0, Kxi, Kxit)
                H2X.append(self.vector_axpy(1.0, Fxi, self.so_to_mo(H2X_so)))

        else:
            for Fxi, Jxi in zip(Fx, Jx):
                if self.triplet:
                    H1X.append(self.vector_copy(Fxi))
                else:
                    H1X_so = self.vector_scale(4.0, Jxi)
                    H1X.append(self.vector_axpy(1.0, Fxi, self.so_to_mo(H1X_so)))
                H2X.append(self.vector_copy(Fxi))
        return H1X, H2X

    def combine_A(self, Fx, Jx, Kx=None):
        """
        Singlet:
           A X = [(Ea - Ei) + 2 J - K] X

        Triplet:
           A X = [(Ea - Ei) - K] X
        """
        Ax = []
        if Kx is not None:
            for Fxi, Jxi, Kxi in zip(Fx, Jx, Kx):
                if self.triplet:
                    Ax.append(self.vector_axpy(1.0, Fxi, self.so_to_mo(self.vector_scale(-1.0, Kxi))))
                else:
                    Ax_so = self.vector_scale(2.0, Jxi)
                    Ax_so = self.vector_axpy(-1.0, Kxi, Ax_so)
                    Ax.append(self.vector_axpy(1.0, Fxi, self.so_to_mo(Ax_so)))

        else:
            for Fxi, Jxi in zip(Fx, Jx):
                if self.triplet:
                    Ax.append(self.vector_copy(Fxi))
                else:
                    Ax.append(self.vector_axpy(1.0, Fxi, self.so_to_mo(self.vector_scale(2.0, Jxi))))
        return Ax

    def compute_products(self, vectors):
        """Given a set of vectors X Compute products
        if ptype == rpa:
           Returns pair (A+B)X, (A-B)X
        if ptype == tda:
           Returns AX
        """

        self.product_cache.reset()
        n_old = self.product_cache.count()
        n_new = len(vectors)

        if n_new <= n_old:
            self.product_cache.reset()
            compute_vectors = vectors
        else:
            compute_vectors = vectors[n_old:]

        # Build base one and two electron quantities
        Fx = self.wfn.onel_Hx(compute_vectors)
        twoel = self.wfn.twoel_Hx(compute_vectors, False, "SO")
        Jx, Kx = self.split_twoel(twoel)

        # Switch between rpa and tda
        if self.ptype == 'rpa':
            H1X_new, H2X_new = self.combine_H1_H2(Fx, Jx, Kx)
            for H1x in H1X_new:
                self.vector_scale(-1.0, H1x)
            for H2x in H2X_new:
                self.vector_scale(-1.0, H2x)

            H1X_all = self.product_cache.add("H1", H1X_new)
            H2X_all = self.product_cache.add("H2", H2X_new)
            return H1X_all, H2X_all

        else:
            AX_new = self.combine_A(Fx, Jx, Kx)
            for Ax in AX_new:
                self.vector_scale(-1.0, Ax)
            AX_all = self.product_cache.add("A", AX_new)
            return AX_all

    def so_to_mo(self, X):
        return core.Matrix.triplet(self.Co, X, self.Cv, True, False, False)

    def split_twoel(self, twoel):
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
        return self.G_gs ^ self.G_es

    def _build_prec(self):
        """
        Builds energy denominator
        """
        self.prec = self.new_vector()
        for h in range(self.wfn.nirrep()):
            self.prec.nph[h][:] = self.E_vir.nph[h ^ self.G_trans] - self.E_occ.nph[h].reshape(-1, 1)

    def precondition(self, Rvec, shift):
        """
        Applies the preconditioner with a shift

        value = R / (shift - preconditioner)
        """
        for h in range(self.wfn.nirrep()):
            den = shift - self.prec.nph[h]
            den[np.abs(den) < 0.0001] = 1.0
            Rvec.nph[h][:] /= den
        return Rvec

    def generate_guess(self, nguess):
        deltas = []
        guess_vectors = []
        for h in range(self.wfn.nirrep()):
            for i, ei in enumerate(self.E_occ.nph[h]):
                for a, ea in enumerate(self.E_vir.nph[h ^ self.G_trans]):
                    deltas.append((ea - ei, i, a, h))
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


class TDUSCFEngine(PairedMatPerVector):
    def __init__(self, wfn, ptype):

        # Primary data
        self.wfn = wfn
        self.ptype = ptype

        # Find product type
        if ptype == 'rpa':
            self.product_cache = ProductCache("H1", "H2")
        else:
            self.product_cache = ProductCache("A")

        # Save orbitals and eigenvalues
        self.Co = [wfn.Ca_subset("SO", "OCC"), wfn.Cb_subset("SO", "OCC")]
        self.Cv = [wfn.Ca_subset("SO", "VIR"), wfn.Cb_subset("SO", "VIR")]
        self.E_occ = [wfn.epsilon_a_subset("SO", "OCC"), wfn.epsilon_b_subset("SO", "OCC")]
        self.E_vir = [wfn.epsilon_a_subset("SO", "VIR"), wfn.epsilon_b_subset("SO", "VIR")]
        self.needs_K_like = self.wfn.functional().is_x_hybrid() or self.wfn.functional().is_x_lrc()

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

    def new_vector(self, name=""):
        return [
            core.Matrix(name + 'a', self.occpi[0], self.virpi[0], self.G_trans),
            core.Matrix(name + 'b', self.occpi[1], self.virpi[1], self.G_trans)
        ]

    def reset_for_state_symm(self, symmetry):
        self.G_es = symmetry
        self._build_prec()
        self.product_cache.reset()

    def combine_H1_H2(self, Fx, Jx, Kx=None):
        """Combine Fx, Jx, Kx (if hybrid dft) to make (A+B)x products"""
        H1X = []
        H2X = []
        if Kx is not None:
            # H1X = Fx + 2Jx - Kx - Kx^T
            # H2X = Fx - Kx + Kx^T
            for Fxi, Jxi, Kxi in zip(Fx, Jx, Kx):
                H1X_so = self.vector_scale(2.0, Jxi)
                Kxit = self.vector_transpose(Kxi)
                H1X_so = self.vector_axpy(-1.0, Kxi, H1X_so)
                H1X_so = self.vector_axpy(-1.0, Kxit, H1X_so)
                H1X.append(self.vector_axpy(1.0, Fxi, self.so_to_mo(H1X_so)))

                H2X_so = self.vector_axpy(-1.0, Kxi, Kxit)
                H2X.append(self.vector_axpy(1.0, Fxi, self.so_to_mo(H2X_so)))
        else:
            # H1X = Fx + 2Jx - Kx - Kx^T
            # H2X = Fx - Kx + Kx^T
            for Fxi, Jxi in zip(Fx, Jx):
                H1X_so = self.vector_scale(2.0, Jxi)
                H1X.append(self.vector_axpy(1.0, Fxi, self.so_to_mo(H1X_so)))
                H2X.append(self.vector_copy(Fxi))

        return H1X, H2X

    def combine_A(self, Fx, Jx, Kx):
        Ax = []
        if Kx is not None:
            # Ax = Fx + J - K
            for Fxi, Jxi, Kxi in zip(Fx, Jx, Kx):
                Ax_so = self.vector_axpy(-1.0, Kxi, Jxi)
                Ax.append(self.vector_axpy(1.0, Fxi, self.so_to_mo(Ax_so)))
        else:
            for Fxi, Jxi in zip(Fx, Jx):
                Ax.append(self.vector_axpy(1.0, Fxi, self.so_to_mo(Jxi)))
        return Ax

    def compute_products(self, vectors):
        """Compute Products for a list of guess vectors (X).
        if ptype == 'rpa':
                           H1 ,   H2
           returns pair (A+B)X, (A-B)X products
        if ptype == 'tda':
           returns Ax products.
        """

        def _elemwise_prod(x, y):
            p = self.vector_copy(x)
            for h in range(self.wfn.nirrep()):
                p[0].nph[h][:] = -1.0 * x[0].nph[h] * y[0].nph[h]
                p[1].nph[h][:] = -1.0 * x[1].nph[h] * y[1].nph[h]
            return p

        #TODO: product cache
        self.product_cache.reset()
        vec_flat = []
        compute_vectors = vectors
        for vec_a, vec_b in compute_vectors:
            vec_flat.append(vec_a)
            vec_flat.append(vec_b)
        #Fx = self.pair_onel(self.wfn.onel_Hx(vec_flat))
        Fx = [_elemwise_prod(v, self.prec) for v in vectors]
        twoel = self.wfn.twoel_Hx(vec_flat, False, "SO")
        Jx, Kx = self.split_twoel(twoel)
        if self.ptype == "rpa":
            H1X_new, H2X_new = self.combine_H1_H2(Fx, Jx, Kx)
            for H1x in H1X_new:
                self.vector_scale(-1.0, H1x)
            for H2x in H2X_new:
                self.vector_scale(-1.0, H2x)
            H1X_all = self.product_cache.add("H1", H1X_new)
            H2X_all = self.product_cache.add("H2", H2X_new)
            return H1X_all, H2X_all
        else:
            AX_new = self.combine_A(Fx, Jx, Kx)
            for Ax in AX_new:
                self.vector_scale(-1.0, Ax)
            AX_all = self.product_cache.add("A", AX_new)
            return AX_all

    def so_to_mo(self, X):
        return [core.Matrix.triplet(self.Co[i], X[i], self.Cv[i], True, False, False) for i in (0, 1)]

    def pair_onel(self, onel):
        return list(zip(onel[0::2], onel[1::2]))

    def split_twoel(self, twoel):
        if self.needs_K_like:
            Jx = list(zip(twoel[0::4], twoel[1::4]))
            Kx = list(zip(twoel[2::4], twoel[3::4]))
        else:
            Jx = list(zip(twoel[0::2], twoel[1::2]))
            Kx = None
        return Jx, Kx

    @property
    def G_trans(self):
        return self.G_gs ^ self.G_es

    def _build_prec(self):
        self.prec = self.new_vector()
        for h in range(self.wfn.nirrep()):
            self.prec[0].nph[h][:] = self.E_vir[0].nph[h ^ self.G_trans] - self.E_occ[0].nph[h].reshape(-1, 1)
            self.prec[1].nph[h][:] = self.E_vir[1].nph[h ^ self.G_trans] - self.E_occ[1].nph[h].reshape(-1, 1)

    def precondition(self, Rvec, shift):
        for h in range(self.wfn.nirrep()):

            den = self.prec[0].nph[h] - shift
            den[np.abs(den) < 0.0001] = 1.0
            Rvec[0].nph[h][:] /= den

            den = self.prec[1].nph[h] - shift
            den[np.abs(den) < 0.0001] = 1.0
            Rvec[1].nph[h][:] /= den
        return Rvec

    def generate_guess(self, nguess):
        guess_vectors = []
        deltas = []
        for h in range(self.wfn.nirrep()):
            for i, ei in enumerate(self.E_occ[0].nph[h]):
                for a, ea in enumerate(self.E_vir[0].nph[h ^ self.G_trans]):
                    deltas.append((ea - ei, 0, i, a, h))
            for i, ei in enumerate(self.E_occ[1].nph[h]):
                for a, ea in enumerate(self.E_vir[1].nph[h ^ self.G_trans]):
                    deltas.append((ea - ei, 1, i, a, h))

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
