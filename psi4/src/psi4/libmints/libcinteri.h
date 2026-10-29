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

#include <cstdint>
#include <memory>
#include <unordered_map>
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
 *  Design (validated against libint2 to ~1e-14 in-tree; ERIs, erf-ERIs and SCF
 *  energies for s..f, cc-pVDZ/TZ and def2-SVP):
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
 *   - Cartesian normalization: libcint's cartesian shell differs from psi4's by
 *     a single per-l scale, which libint2_to_libcint measures on the axial
 *     (l,0,0) component, as Libint2 normalizes.
 *   - libcint writes col-major (first index fastest); psi4 wants row-major
 *     (first index slowest). Both are handled in the same repack.
 *
 *  Current scope: 4-center (ab|cd), spherical, cartesian, and mixed basis sets, with
 *  first and second nuclear derivatives, plus the range-separated erf/erfc
 *  variants (env[PTR_RANGE_OMEGA]).
 *   - Derivatives: libcint differentiates the function in its first slot
 *     (int2e_ip1: (nabla a b|c d); int2e_ipip1: (nabla nabla a b|c d);
 *     int2e_ipvip1: (nabla a nabla b|c d); int2e_ip1ip2: (nabla a b|nabla c d)),
 *     so each center's derivative comes from one call with the quartet permuted
 *     to put that center first. A nuclear derivative is minus the electronic
 *     one, so first derivatives change sign and second derivatives do not. The
 *     buffers follow Libint2's layout exactly: derivatives with respect to the
 *     non-dummy centers only, in center order (12/9/6 first derivatives for
 *     (ab|cd), (a|cd) or (ab|c), (a|c)), second derivatives as the row-major
 *     upper triangle over those coordinates (78/45/21), and the remaining
 *     buffers of the 12/78 psi4 exposes pointing at zeros.
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
    /// libcint CINTOpt* for each kernel (held opaquely so the header stays free
    /// of <cint.h>): (ab|cd), and when deriv_ allows, ip1, ipip1, ipvip1, ip1ip2.
    /// libcint only reads them while computing, so clones share them.
    std::shared_ptr<void> opt_;
    std::shared_ptr<void> opt_ip1_;
    std::shared_ptr<void> opt_ipip1_;
    std::shared_ptr<void> opt_ipvip1_;
    std::shared_ptr<void> opt_ip1ip2_;

    /// Global libcint bas-index of shell 0 of each of the four psi4 bases.
    int bas_start_[4];
    /// Deduplicated bases and their libcint bas starts (for sieve lookup).
    std::vector<std::pair<const BasisSet *, int>> basis_starts_;

    /// Owned buffer holding one shell quartet in psi4 order.
    std::vector<double> target_store_;
    /// Scratch for libcint's (col-major) output before repack into target_.
    std::vector<double> cint_buf_;
    /// libcint's working memory, which it would otherwise allocate and free on every call. libcint never checks the
    /// size of what it is handed, so cache_sizes_ keeps libcint's own count of what each kind of quartet needs, keyed
    /// by kernel and by each slot's angular momentum and number of primitives.
    std::vector<double> cache_;
    std::unordered_map<uint64_t, size_t> cache_sizes_;
    /// Derivative buffers (12 or 78 quartet-sized blocks), a block of zeros for
    /// the buffers Libint2 leaves empty, and a sink for redundant components.
    std::vector<double> deriv_store_;
    std::vector<double> zeros_;
    std::vector<double> discard_;
    /// Size of one quartet block of the factory's four bases.
    size_t quartet_size_;

    /// The centers (0-3) whose basis is not the density-fitting dummy, in order.
    std::vector<int> active_;

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
    /// Scratch for transform_mixed, as large as a cartesian quartet block.
    std::vector<double> mixed_scratch_;

    void common_init();
    /// Size the quartet buffers and, for a mix, the solid-harmonic transforms.
    void allocate_buffers();

    /// Build the libcint atm/bas/env arrays from the four psi4 basis sets.
    void build_environment();

    /// Append one psi4 basis set's atoms and shells to atm_/bas_/env_, return its bas start.
    /// coef_ratio[l] converts Libint2 contraction coefficients to libcint's convention.
    int append_basis(const BasisSet &bs, const std::vector<double> &coef_ratio);
    /// Raise env_[PTR_EXPCUTOFF] if libcint's pair screening would drop the
    /// integrals of this basis's diffuse high-l primitives.
    void set_expcutoff();

    /// The libcint two-electron kernels used: (ab|cd) and its derivatives.
    enum class Kernel { ERI, IP1, IPIP1, IPVIP1, IP1IP2 };

    /// Evaluate one libcint kernel (sph or cart per cart_) for the shells at
    /// libcint bas indices g with angular momenta l, placing psi4 center perm[t]
    /// in libcint slot t, and scatter component c, scaled by sign, into the
    /// psi4-ordered block out[c]. Returns whether libcint found any nonzero
    /// integral (the blocks are zeroed otherwise).
    bool run_kernel(Kernel kernel, const int g[4], const int l[4], const int perm[4], double sign,
                    double *const *out);

    /// Transform a row-major cartesian quartet in psi4 center order, in place, for the
    /// spherical shells among the libcint bas indices g (mixed_ only).
    void transform_mixed(double *data, const int g[4], const int l[4]);
    /// Number of psi4 functions of the shell at libcint bas index g, of angular momentum l.
    int nfunction(int g, int l) const;

    /// Compute one shell quartet (given libcint bas indices and angular momenta)
    /// into target_full_ in psi4 order. Returns the number of integrals computed.
    size_t compute_quartet(int g1, int g2, int g3, int g4, int l1, int l2, int l3, int l4);

    /// libcint bas indices and angular momenta of a psi4 shell quartet.
    void quartet(int s1, int s2, int s3, int s4, int g[4], int l[4]) const;

    size_t compute_shell_for_sieve(const std::shared_ptr<BasisSet> bs, int s1, int s2, int s3, int s4,
                                   bool is_bra) override;

   public:
    LibcintTwoElectronInt(const IntegralFactory *integral, int deriv = 0, double omega = 0.0,
                          bool use_shell_pairs = false, bool needs_exchange = false);
    LibcintTwoElectronInt(const LibcintTwoElectronInt &rhs);
    ~LibcintTwoElectronInt() override = default;

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
