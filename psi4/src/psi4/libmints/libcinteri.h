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
#include "psi4/libmints/integral.h"

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
 *  Design (validated against libint2 to ~1e-13 in-tree; see tests/pytests/test_libcint.py):
 *   - Requests spherical (int2e_sph) or cartesian (int2e_cart) integrals to match
 *     the bases. When they mix (e.g., cartesian primary, spherical auxiliary), all
 *     are requested cartesian and the spherical shells transformed with psi4's
 *     solid harmonics (transform_mixed).
 *   - Basis functions are exactly those of the Libint2 path: each shell's Libint2
 *     contraction coefficients (BasisSet::l2_shell) times a per-l factor between
 *     the two conventions (libint2_to_libcint, measured with libcint's own
 *     overlap). Bases with deliberately non-normalized coefficients (diffuse
 *     external charges, SAP potentials) thus carry over unchanged.
 *   - Component ordering: for spherical shells libcint uses m = -l..+l for
 *     l != 1 and the cartesian order (px,py,pz) for l == 1, while psi4 (libint2,
 *     Gaussian solid-harmonic ordering) uses m = 0,+1,-1,+2,-2,...; for
 *     cartesian shells (int2e_cart) libcint's order already equals psi4's
 *     CartesianIter order, so the map is the identity. All relative signs +1.
 *     (This assumes libcint built without -DPYPZPX.)
 *   - libcint writes col-major (first index fastest); psi4 wants row-major
 *     (first index slowest). Both are handled in the same repack.
 *
 *  Current scope: 4-center (ab|cd), spherical, cartesian, and mixed basis sets,
 *  deriv=0, plus the range-separated erf/erfc variants (env[PTR_RANGE_OMEGA]).
 *  Density-fitting (2-/3-center) works too: psi4 passes the absent center as a dummy s-shell
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
    /// libcint's working cache, sized for the most demanding quartet.
    std::vector<double> cache_;

    /// Range-separation parameter written to env[PTR_RANGE_OMEGA]
    /// (0 = full Coulomb, >0 = erf/long-range, <0 = erfc/short-range).
    double omega_;

    /// Whether libcint is asked for cartesian (int2e_cart) or spherical (int2e_sph) integrals.
    bool cart_;
    /// Whether the bases mix cartesian and spherical shells (e.g., a cartesian primary
    /// with a spherical auxiliary). Then cart_ is true and the spherical shells are
    /// transformed here, shell by shell, with psi4's own solid harmonics.
    bool mixed_;
    /// Per libcint bas index, whether psi4 wants the shell spherical (used when mixed_).
    std::vector<char> bas_pure_;
    /// psi4 cartesian -> spherical transforms by angular momentum (used when mixed_). Our
    /// own, rather than TwoBodyAOInt::pure_transform's from the IntegralFactory, since the
    /// sieve's quartets come from one basis, not bs1_..bs4_, and the factory's stop at l=8.
    std::vector<SphericalTransform> sph_trans_;

    /// Set up the libcint optimizer and buffers (fresh objects and clones alike).
    void common_init();

    /// Build the libcint atm/bas/env arrays from the four psi4 basis sets.
    void build_environment();

    /// Append one psi4 basis set's atoms and shells to atm_/bas_/env_, return its bas start.
    /// coef_ratio[l] converts Libint2 contraction coefficients to libcint's convention.
    int append_basis(const BasisSet &bs, const std::vector<double> &coef_ratio);

    /// Transform a row-major cartesian quartet in place of target_full_ for the spherical
    /// shells among g1..g4 (mixed_ only). Returns the final number of integrals.
    size_t transform_mixed(const int *g, const int *l);

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

    void initialize_sieve() override;

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
