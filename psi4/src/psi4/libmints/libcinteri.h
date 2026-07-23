/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi4_libmints_libcinteri_h
#define _psi4_libmints_libcinteri_h

#include <memory>
#include <vector>

#include "psi4/libmints/twobody.h"

namespace psi {

class BasisSet;
class GaussianShell;
class IntegralFactory;

/*! \ingroup MINTS
 *  \class LibcintTwoElectronInt
 *  \brief Two-electron repulsion integrals via libcint (Sun, J. Comput. Chem. 2015).
 *
 *  Experimental optional backend, selected with INTEGRAL_PACKAGE LIBCINT.
 *
 *  Design (validated against libint2 to ~1e-14 in-tree; ERIs, erf-ERIs and SCF
 *  energies for s..f, cc-pVDZ/TZ and def2-SVP):
 *   - Requests spherical integrals from libcint (int2e_sph).
 *   - Basis functions are reproduced exactly by feeding libcint the *original*
 *     coefficients scaled by CINTgto_norm, then rescaling each contraction to
 *     unit self-overlap via libcint's own int1e_ovlp (normalize_shells) --
 *     matching libint2's embed_normalization, including split general
 *     contractions.
 *   - Component ordering: for spherical shells libcint uses m = -l..+l for
 *     l != 1 and the cartesian order (px,py,pz) for l == 1, while psi4 (libint2,
 *     Gaussian solid-harmonic ordering) uses m = 0,+1,-1,+2,-2,...; for
 *     cartesian shells (int2e_cart) libcint's order already equals psi4's
 *     CartesianIter order, so the map is the identity. All relative signs +1.
 *   - Cartesian normalization: libcint's cartesian shell differs from psi4's by
 *     a single per-shell scale; normalizing the axial (l,0,0) self-overlap to 1
 *     via int1e_ovlp_cart (see normalize_shells) cancels it.
 *   - libcint writes col-major (first index fastest); psi4 wants row-major
 *     (first index slowest). Both are handled in the same repack.
 *
 *  Current scope: 4-center (ab|cd), spherical and cartesian basis sets,
 *  deriv=0, plus the range-separated erf/erfc variants (env[PTR_RANGE_OMEGA]).
 *  Density-fitting
 *  (2-/3-center) works too: psi4 passes the absent center as a dummy s-shell
 *  (l=0, exp=0) which is fed to libcint as a bare constant matching libint2's
 *  unit shell, so int2e_sph over the dummy yields (ij|k)/(i|k) integrals
 *  bit-identical to the Libint2 path. (A first-class native 2c/3c interface is
 *  the intended long-term replacement for that dummy-shell route.)
 */
class LibcintTwoElectronInt : public TwoBodyAOInt {
   protected:
    /// libcint environment describing all shells of the (deduplicated) bases.
    std::vector<int> atm_;
    std::vector<int> bas_;
    std::vector<double> env_;
    /// libcint CINTOpt* (held opaquely so the header stays free of <cint.h>).
    void *opt_;

    /// Global libcint bas-index of shell 0 of each of the four psi4 bases.
    int bas_start_[4];
    /// Deduplicated bases and their libcint bas starts (for sieve lookup).
    std::vector<std::pair<const BasisSet *, int>> basis_starts_;

    /// Owned buffer holding one shell quartet in psi4 order.
    std::vector<double> target_store_;
    /// Scratch for libcint's (col-major) output before repack into target_.
    std::vector<double> cint_buf_;

    /// Range-separation parameter written to env[PTR_RANGE_OMEGA]
    /// (0 = full Coulomb, >0 = erf/long-range, <0 = erfc/short-range).
    double omega_;

    /// Whether the (non-dummy) bases are cartesian (int2e_cart) or spherical.
    bool cart_;

    void common_init();

    /// Build the libcint atm/bas/env arrays from the four psi4 basis sets.
    void build_environment();

    /// Append one psi4 basis set's shells to bas_/env_, return its bas start.
    int append_basis(const BasisSet &bs);

    /// Rescale each shell's contraction to unit self-overlap (matches libint2).
    void normalize_shells();

    /// Compute one shell quartet (given libcint bas indices and angular momenta)
    /// into target_full_ in psi4 order. Returns the number of integrals computed.
    size_t compute_quartet(int g1, int g2, int g3, int g4, int l1, int l2, int l3, int l4);

    size_t compute_shell_for_sieve(const std::shared_ptr<BasisSet> bs, int s1, int s2, int s3, int s4,
                                   bool is_bra) override;

   public:
    LibcintTwoElectronInt(const IntegralFactory *integral, int deriv = 0, double omega = 0.0,
                          bool use_shell_pairs = false, bool needs_exchange = false);
    LibcintTwoElectronInt(const LibcintTwoElectronInt &rhs);
    ~LibcintTwoElectronInt() override;

    size_t compute_shell(const AOShellCombinationsIterator &) override;
    size_t compute_shell(int s1, int s2, int s3, int s4) override;

    size_t compute_shell_deriv1(int, int, int, int) override;
    size_t compute_shell_deriv2(int, int, int, int) override;

    /// Single-shell blocks, mirroring the Libint2 backend.
    void compute_shell_blocks(int shellpair12, int shellpair34, int npair12 = -1, int npair34 = -1) override;
};

class LibcintERI : public LibcintTwoElectronInt {
   public:
    LibcintERI(const IntegralFactory *integral, int deriv = 0, bool use_shell_pairs = false,
               bool needs_exchange = false);
    LibcintERI(const LibcintERI &rhs);
    LibcintERI *clone() const override { return new LibcintERI(*this); }
};

class LibcintErfERI : public LibcintTwoElectronInt {
   public:
    LibcintErfERI(double omega, const IntegralFactory *integral, int deriv = 0, bool use_shell_pairs = false,
                  bool needs_exchange = false);
    LibcintErfERI(const LibcintErfERI &rhs);
    LibcintErfERI *clone() const override { return new LibcintErfERI(*this); }
};

class LibcintErfComplementERI : public LibcintTwoElectronInt {
   public:
    LibcintErfComplementERI(double omega, const IntegralFactory *integral, int deriv = 0,
                            bool use_shell_pairs = false, bool needs_exchange = false);
    LibcintErfComplementERI(const LibcintErfComplementERI &rhs);
    LibcintErfComplementERI *clone() const override { return new LibcintErfComplementERI(*this); }
};

}  // namespace psi

#endif
