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

#include "psi4/libmints/libcinteri.h"

#include <algorithm>

#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/gshell.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libpsi4util/exception.h"

#include <cmath>

extern "C" {
#include <cint.h>
int int2e_sph(double *out, int *dims, int *shls, int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt,
              double *cache);
void int2e_optimizer(CINTOpt **opt, int *atm, int natm, int *bas, int nbas, double *env);
int int1e_ovlp_sph(double *out, int *dims, int *shls, int *atm, int natm, int *bas, int nbas, double *env,
                   CINTOpt *opt, double *cache);
double CINTgto_norm(int n, double a);
}

namespace psi {

namespace {

/// Is this psi4 basis the density-fitting dummy (BasisSet::zero_ao_basis_set)?
bool is_dummy_basis(const BasisSet &bs) {
    return bs.nshell() == 1 && bs.nprimitive() == 1 && bs.shell(0).am() == 0 && bs.shell(0).exp(0) == 0.0;
}

/// psi4 / libint2-Gaussian spherical slot for magnetic quantum number m
/// (m == 0 -> 0, m > 0 -> 2m-1, m < 0 -> -2m). All relative signs are +1
/// (machine-verified against libint2 in Gaussian ordering, l = 0..4).
inline int psi4_slot(int m) { return m == 0 ? 0 : (m > 0 ? 2 * m - 1 : -2 * m); }

/// Magnetic quantum number of libcint spherical slot i for an l-shell. libcint
/// orders l >= 2 as m = -l..+l (so m = i - l), but p functions (l == 1) are the
/// special cartesian-like order (px, py, pz) = m(+1, -1, 0). Verified against
/// libint2 (Gaussian ordering) to ~1e-14 for l = 0..4.
inline int libcint_m(int i, int l) {
    if (l == 1) return i == 0 ? +1 : (i == 1 ? -1 : 0);
    return i - l;
}

}  // namespace

LibcintTwoElectronInt::LibcintTwoElectronInt(const IntegralFactory *integral, int deriv, double omega,
                                             bool use_shell_pairs, bool needs_exchange)
    : TwoBodyAOInt(integral, deriv), opt_(nullptr), omega_(omega) {
    if (deriv_ != 0)
        throw PSIEXCEPTION("LIBCINT backend: integral derivatives are not implemented (use INTEGRAL_PACKAGE LIBINT2).");

    // Spherical harmonics only for now; cartesian normalization differs from
    // psi4's and would need a separate per-component rescale. Density fitting
    // is supported by a simple jerry-rig: psi4 passes the "absent" center as a
    // dummy s-shell (l=0, exp=0, coef=1) -- exactly libint2's unit shell, a
    // constant function -- so the 4-center int2e_sph path with that dummy shell
    // already yields the (ij|k) 3-center and (i|k) 2-center integrals.
    for (const auto *bs : {original_bs1_.get(), original_bs2_.get(), original_bs3_.get(), original_bs4_.get()}) {
        if (is_dummy_basis(*bs)) continue;  // handled shell-by-shell below
        for (int s = 0; s < bs->nshell(); ++s)
            if (bs->shell(s).is_cartesian())
                throw PSIEXCEPTION("LIBCINT backend: cartesian basis sets are not yet supported (spherical only).");
    }

    build_environment();
    common_init();
}

LibcintTwoElectronInt::LibcintTwoElectronInt(const LibcintTwoElectronInt &rhs)
    : TwoBodyAOInt(rhs),
      atm_(rhs.atm_),
      bas_(rhs.bas_),
      env_(rhs.env_),
      opt_(nullptr),
      basis_starts_(rhs.basis_starts_),
      omega_(rhs.omega_) {
    bas_start_[0] = rhs.bas_start_[0];
    bas_start_[1] = rhs.bas_start_[1];
    bas_start_[2] = rhs.bas_start_[2];
    bas_start_[3] = rhs.bas_start_[3];
    common_init();
}

LibcintTwoElectronInt::~LibcintTwoElectronInt() {
    if (opt_) {
        CINTOpt *o = static_cast<CINTOpt *>(opt_);
        CINTdel_optimizer(&o);
    }
}

void LibcintTwoElectronInt::common_init() {
    // Range separation: env[PTR_RANGE_OMEGA] > 0 -> erf (long-range), < 0 -> erfc.
    env_[PTR_RANGE_OMEGA] = omega_;

    const int nbas = static_cast<int>(bas_.size() / BAS_SLOTS);
    const int natm = static_cast<int>(atm_.size() / ATM_SLOTS);
    CINTOpt *o = nullptr;
    int2e_optimizer(&o, atm_.data(), natm, bas_.data(), nbas, env_.data());
    opt_ = o;

    // Size scratch for the largest quartet block we may compute. This must
    // cover not only the factory's (bs1,bs2,bs3,bs4) block but also the Schwarz
    // sieve, which computes (MN|MN) over a single basis -- so a block up to
    // (2*maxam+1)^4 for the largest-am basis (e.g. the auxiliary in DF, whose
    // (MN|MN) blocks dwarf the aux*dummy*prim*prim factory block).
    int maxam = 0;
    for (const auto *bs : {original_bs1_.get(), original_bs2_.get(), original_bs3_.get(), original_bs4_.get()})
        maxam = std::max(maxam, bs->max_am());
    const size_t nmax = static_cast<size_t>(2 * maxam + 1);
    const size_t blk = nmax * nmax * nmax * nmax;
    target_store_.assign(blk, 0.0);
    cint_buf_.assign(blk, 0.0);
    target_full_ = target_store_.data();
    source_full_ = nullptr;
    buffers_.resize(1, target_full_);

    setup_sieve();
    create_blocks();
}

int LibcintTwoElectronInt::append_basis(const BasisSet &bs) {
    const int start = static_cast<int>(bas_.size() / BAS_SLOTS);
    for (int s = 0; s < bs.nshell(); ++s) {
        const GaussianShell &sh = bs.shell(s);
        const int nprim = sh.nprimitive();

        const int exp_off = static_cast<int>(env_.size());
        for (int p = 0; p < nprim; ++p) env_.push_back(sh.exp(p));

        // libcint takes env coefficients as multipliers of the *unnormalized*
        // primitive, expecting the primitive normalization baked in via
        // CINTgto_norm (the standard pyscf convention). psi4's *original*
        // coefficients are the same normalized-contraction coefficients psi4
        // hands libint2 (basisset.cc update_l2_shells), so
        //   env_coef_p = original_coef_p * CINTgto_norm(l, a_p)
        // reproduces libint2's embed-normalized basis function exactly
        // (verified in-tree: unit shell self-overlap, ERIs to ~1e-13).
        const int l = sh.am();
        const int coef_off = static_cast<int>(env_.size());
        const bool dummy = (nprim == 1 && l == 0 && sh.exp(0) == 0.0);
        if (dummy) {
            // Unit shell = the bare constant function "1", matching
            // libint2::Shell::unit() (its renorm() skips alpha==0, leaving
            // coeff=1). CINTgto_norm is singular at exp=0 and cannot be used, so
            // it is not applied here -- but libcint's spherical transform still
            // folds the Y_00 = 1/sqrt(4*pi) normalization into every s shell.
            // Feed sqrt(4*pi) to cancel it, so the dummy is a true constant 1 and
            // the DF metric (P|Q) and 3-center (Q|mn) are *identical* to the
            // Libint2 path, not merely proportional. A global scale on the metric
            // would otherwise change which vectors a fixed linear-dependence
            // threshold discards when it is inverted.
            env_.push_back(std::sqrt(4.0 * M_PI));
        } else {
            for (int p = 0; p < nprim; ++p) env_.push_back(sh.original_coef(p) * CINTgto_norm(l, sh.exp(p)));
        }

        bas_.push_back(bs.shell_to_center(s));  // ATOM_OF
        bas_.push_back(sh.am());                // ANG_OF
        bas_.push_back(nprim);                  // NPRIM_OF
        bas_.push_back(1);                      // NCTR_OF
        bas_.push_back(0);                      // KAPPA_OF
        bas_.push_back(exp_off);                // PTR_EXP
        bas_.push_back(coef_off);               // PTR_COEFF
        bas_.push_back(0);                      // RESERVE
    }
    return start;
}

void LibcintTwoElectronInt::build_environment() {
    env_.assign(PTR_ENV_START, 0.0);

    // Atoms: all four bases share one molecule (primary/auxiliary on the same
    // nuclei), so one atm/coord table indexed by shell_to_center() suffices.
    auto mol = original_bs1_->molecule();
    const int natom = mol->natom();
    for (int a = 0; a < natom; ++a) {
        const int coord_off = static_cast<int>(env_.size());
        env_.push_back(mol->x(a));
        env_.push_back(mol->y(a));
        env_.push_back(mol->z(a));
        atm_.push_back(static_cast<int>(mol->true_atomic_number(a)));  // CHARGE_OF
        atm_.push_back(coord_off);                                     // PTR_COORD
        atm_.push_back(0);                                             // NUC_MOD_OF
        atm_.push_back(0);                                             // PTR_ZETA
        atm_.push_back(0);
        atm_.push_back(0);
    }

    const BasisSet *bs[4] = {original_bs1_.get(), original_bs2_.get(), original_bs3_.get(), original_bs4_.get()};
    for (int i = 0; i < 4; ++i) {
        int start = -1;
        for (const auto &kv : basis_starts_)
            if (kv.first == bs[i]) start = kv.second;
        if (start < 0) {
            start = append_basis(*bs[i]);
            basis_starts_.emplace_back(bs[i], start);
        }
        bas_start_[i] = start;
    }

    normalize_shells();
}

void LibcintTwoElectronInt::normalize_shells() {
    // Rescale each shell's contraction to unit self-overlap, matching libint2's
    // embed_normalization. Doing this through libcint's own int1e_ovlp keeps us
    // convention-agnostic: generally-contracted bases that psi4 splits into
    // segmented shells arrive here un-normalized, and this restores unit norm
    // exactly as the Libint2 path does. Single-primitive shells are already
    // unit-normalized via CINTgto_norm, so they are left effectively unchanged.
    const int nbas = static_cast<int>(bas_.size() / BAS_SLOTS);
    const int natm = static_cast<int>(atm_.size() / ATM_SLOTS);
    std::vector<double> ov;
    for (int g = 0; g < nbas; ++g) {
        const int l = bas_[g * BAS_SLOTS + ANG_OF];
        const int nprim = bas_[g * BAS_SLOTS + NPRIM_OF];
        const int coef_off = bas_[g * BAS_SLOTS + PTR_COEFF];
        const int exp_off = bas_[g * BAS_SLOTS + PTR_EXP];
        // Skip the dummy/unit shell (exp=0): its self-overlap diverges and it
        // must stay the bare constant to mirror libint2's unit shell.
        if (nprim == 1 && l == 0 && env_[exp_off] == 0.0) continue;
        const int nf = 2 * l + 1;
        ov.assign(static_cast<size_t>(nf) * nf, 0.0);
        int shls[2] = {g, g};
        int1e_ovlp_sph(ov.data(), nullptr, shls, atm_.data(), natm, bas_.data(), nbas, env_.data(), nullptr, nullptr);
        const double s = ov[0];  // diagonal self-overlap (equal for all components)
        if (s > 0.0) {
            const double inv = 1.0 / std::sqrt(s);
            for (int p = 0; p < nprim; ++p) env_[coef_off + p] *= inv;
        }
    }
}

size_t LibcintTwoElectronInt::compute_quartet(int g1, int g2, int g3, int g4, int n1, int n2, int n3, int n4) {
    const int l1 = (n1 - 1) / 2, l2 = (n2 - 1) / 2, l3 = (n3 - 1) / 2, l4 = (n4 - 1) / 2;
    const size_t n = static_cast<size_t>(n1) * n2 * n3 * n4;
    if (cint_buf_.size() < n) cint_buf_.assign(n, 0.0);

    int shls[4] = {g1, g2, g3, g4};
    const int nbas = static_cast<int>(bas_.size() / BAS_SLOTS);
    const int natm = static_cast<int>(atm_.size() / ATM_SLOTS);
    int ok = int2e_sph(cint_buf_.data(), nullptr, shls, atm_.data(), natm, bas_.data(), nbas, env_.data(),
                       static_cast<CINTOpt *>(opt_), nullptr);
    double *tgt = target_full_;
    if (!ok) {
        std::fill(tgt, tgt + n, 0.0);
        return 0;
    }

    // Repack: libcint is col-major (index 1 fastest) in m = -l..+l order;
    // psi4 is row-major (index 1 slowest) in Gaussian-m order. Signs all +1.
    for (int a = 0; a < n1; ++a) {
        const int pa = psi4_slot(libcint_m(a, l1));
        for (int b = 0; b < n2; ++b) {
            const int pb = psi4_slot(libcint_m(b, l2));
            for (int c = 0; c < n3; ++c) {
                const int pc = psi4_slot(libcint_m(c, l3));
                for (int d = 0; d < n4; ++d) {
                    const int pd = psi4_slot(libcint_m(d, l4));
                    const size_t iC =
                        static_cast<size_t>(a) + n1 * (static_cast<size_t>(b) + n2 * (static_cast<size_t>(c) + static_cast<size_t>(n3) * d));
                    const size_t iP = (((static_cast<size_t>(pa) * n2 + pb) * n3 + pc) * n4 + pd);
                    tgt[iP] = cint_buf_[iC];
                }
            }
        }
    }
    return n;
}

size_t LibcintTwoElectronInt::compute_shell(const AOShellCombinationsIterator &shellIter) {
    return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
}

size_t LibcintTwoElectronInt::compute_shell(int s1, int s2, int s3, int s4) {
    const int n1 = original_bs1_->shell(s1).nfunction();
    const int n2 = original_bs2_->shell(s2).nfunction();
    const int n3 = original_bs3_->shell(s3).nfunction();
    const int n4 = original_bs4_->shell(s4).nfunction();
    curr_buff_size_ = n1 * n2 * n3 * n4;

    const size_t ret = compute_quartet(bas_start_[0] + s1, bas_start_[1] + s2, bas_start_[2] + s3, bas_start_[3] + s4,
                                       n1, n2, n3, n4);
    buffers_[0] = target_full_;
    return ret;
}

size_t LibcintTwoElectronInt::compute_shell_for_sieve(const std::shared_ptr<BasisSet> bs, int s1, int s2, int s3,
                                                      int s4, bool /*is_bra*/) {
    int start = -1;
    for (const auto &kv : basis_starts_)
        if (kv.first == bs.get()) start = kv.second;
    if (start < 0) throw PSIEXCEPTION("LIBCINT backend: sieve basis not found in libcint environment.");

    const int n1 = bs->shell(s1).nfunction();
    const int n2 = bs->shell(s2).nfunction();
    const int n3 = bs->shell(s3).nfunction();
    const int n4 = bs->shell(s4).nfunction();
    curr_buff_size_ = n1 * n2 * n3 * n4;

    const size_t ret = compute_quartet(start + s1, start + s2, start + s3, start + s4, n1, n2, n3, n4);
    buffers_[0] = target_full_;
    return ret;
}

void LibcintTwoElectronInt::compute_shell_blocks(int shellpair12, int shellpair34, int /*npair12*/, int /*npair34*/) {
    // This engine does not batch shells, so each block is a single quartet.
    const int s1 = blocks12_[shellpair12][0].first;
    const int s2 = blocks12_[shellpair12][0].second;
    const int s3 = blocks34_[shellpair34][0].first;
    const int s4 = blocks34_[shellpair34][0].second;
    compute_shell(s1, s2, s3, s4);
}

size_t LibcintTwoElectronInt::compute_shell_deriv1(int, int, int, int) {
    throw PSIEXCEPTION("LIBCINT backend: gradients are not implemented (use INTEGRAL_PACKAGE LIBINT2).");
}

size_t LibcintTwoElectronInt::compute_shell_deriv2(int, int, int, int) {
    throw PSIEXCEPTION("LIBCINT backend: Hessians are not implemented (use INTEGRAL_PACKAGE LIBINT2).");
}

//////////////////////////////////////////////////////////////////////////////

LibcintERI::LibcintERI(const IntegralFactory *integral, int deriv, bool use_shell_pairs, bool needs_exchange)
    : LibcintTwoElectronInt(integral, deriv, 0.0, use_shell_pairs, needs_exchange) {}
LibcintERI::LibcintERI(const LibcintERI &rhs) : LibcintTwoElectronInt(rhs) {}

LibcintErfERI::LibcintErfERI(double omega, const IntegralFactory *integral, int deriv, bool use_shell_pairs,
                             bool needs_exchange)
    : LibcintTwoElectronInt(integral, deriv, omega, use_shell_pairs, needs_exchange) {}
LibcintErfERI::LibcintErfERI(const LibcintErfERI &rhs) : LibcintTwoElectronInt(rhs) {}

LibcintErfComplementERI::LibcintErfComplementERI(double omega, const IntegralFactory *integral, int deriv,
                                                 bool use_shell_pairs, bool needs_exchange)
    : LibcintTwoElectronInt(integral, deriv, -omega, use_shell_pairs, needs_exchange) {}
LibcintErfComplementERI::LibcintErfComplementERI(const LibcintErfComplementERI &rhs) : LibcintTwoElectronInt(rhs) {}

}  // namespace psi
