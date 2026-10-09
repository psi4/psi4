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
#include "psi4/libqt/qt.h"

#include <cmath>

#include <libint2/shell.h>

extern "C" {
#include <cint_funcs.h>
}

// The atm/bas/shls arrays here are plain int, as is libcint's FINT unless it was built with -DI8.
static_assert(sizeof(FINT) == sizeof(int), "LIBCINT backend requires libcint built without I8 (32-bit FINT).");

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

inline int ncart(int l) { return (l + 1) * (l + 2) / 2; }

/// libcint's default primitive screening cutoff (EXPCUTOFF in its config.h, not installed).
constexpr double EXPCUTOFF_DEFAULT = 60.0;

/// Factor converting a Libint2 contraction coefficient into the libcint coefficient
/// describing the same function, for angular momentum l. Both conventions scale a
/// primitive's coefficient as a^((2l+3)/4), so the ratio of the unit-normalized
/// single-primitive coefficients (here at exponent 1) holds for every exponent.
/// libcint's is measured with its own overlap integral (cartesian: the axial
/// component, as Libint2 normalizes), keeping us agnostic of its l=0,1 factors.
double libint2_to_libcint(int l, bool cart) {
    std::vector<int> atm = {0, PTR_ENV_START, 0, 0, 0, 0};
    std::vector<int> bas = {0, l, 1, 1, 0, PTR_ENV_START + 3, PTR_ENV_START + 4, 0};
    std::vector<double> env(PTR_ENV_START + 5, 0.0);
    env[PTR_ENV_START + 3] = 1.0;
    env[PTR_ENV_START + 4] = CINTgto_norm(l, 1.0);
    std::vector<double> ov(static_cast<size_t>(ncart(l)) * ncart(l));
    int shls[2] = {0, 0};
    if (cart)
        int1e_ovlp_cart(ov.data(), nullptr, shls, atm.data(), 1, bas.data(), 1, env.data(), nullptr, nullptr);
    else
        int1e_ovlp_sph(ov.data(), nullptr, shls, atm.data(), 1, bas.data(), 1, env.data(), nullptr, nullptr);
    const double libcint_unit = env[PTR_ENV_START + 4] / std::sqrt(ov[0]);
    const libint2::Shell unit_prim(libint2::svector<double>{1.0}, {libint2::Shell::Contraction{l, !cart, {1.0}}},
                                   {{0.0, 0.0, 0.0}});
    const double libint2_unit = unit_prim.contr[0].coeff[0];
    return libcint_unit / libint2_unit;
}

}  // namespace

LibcintTwoElectronInt::LibcintTwoElectronInt(const IntegralFactory *integral, int deriv, double omega,
                                             bool use_shell_pairs, bool needs_exchange)
    : TwoBodyAOInt(integral, deriv), opt_(nullptr), omega_(omega) {
    if (deriv_ != 0)
        throw PSIEXCEPTION("LIBCINT backend: integral derivatives are not implemented (use INTEGRAL_PACKAGE LIBINT2).");

    // Determine whether the (non-dummy) bases are cartesian, spherical, or a mix.
    // Bases can differ (e.g., a cartesian 6-31G* primary with a spherical RI
    // auxiliary in DF-MP2), but libcint picks int2e_cart or int2e_sph per call, not
    // per shell. A mix is computed cartesian and the spherical shells transformed
    // here (transform_mixed). The l=0 dummy shell is identical either way. Density
    // fitting is handled by the dummy s-shell (l=0, exp=0) that psi4 places in the
    // absent center's slot -- see append_basis.
    bool any_cart = false, any_pure = false;
    for (const auto *bs : {original_bs1_.get(), original_bs2_.get(), original_bs3_.get(), original_bs4_.get()}) {
        if (is_dummy_basis(*bs)) continue;
        for (int s = 0; s < bs->nshell(); ++s) {
            if (bs->shell(s).am() == 0) continue;
            if (bs->shell(s).is_cartesian())
                any_cart = true;
            else
                any_pure = true;
        }
    }
    cart_ = any_cart;
    mixed_ = any_cart && any_pure;

    timer_on("LibcintTwoElectronInt::LibcintTwoElectronInt");
    build_environment();
    common_init();
    // Clones copy the sieve and blocks from rhs (TwoBodyAOInt copy constructor),
    // so only a fresh object computes them.
    setup_sieve();
    create_blocks();
    timer_off("LibcintTwoElectronInt::LibcintTwoElectronInt");
}

LibcintTwoElectronInt::LibcintTwoElectronInt(const LibcintTwoElectronInt &rhs)
    : TwoBodyAOInt(rhs),
      atm_(rhs.atm_),
      bas_(rhs.bas_),
      env_(rhs.env_),
      opt_(nullptr),
      basis_starts_(rhs.basis_starts_),
      omega_(rhs.omega_),
      cart_(rhs.cart_),
      mixed_(rhs.mixed_),
      bas_pure_(rhs.bas_pure_) {
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
    const size_t nmax = static_cast<size_t>(cart_ ? ncart(maxam) : 2 * maxam + 1);
    const size_t blk = nmax * nmax * nmax * nmax;
    target_store_.assign(blk, 0.0);
    cint_buf_.assign(blk, 0.0);
    if (mixed_) {
        mixed_buf_.assign(blk, 0.0);
        sph_trans_.clear();
        for (int l = 0; l <= maxam; ++l) sph_trans_.emplace_back(l);
    }
    target_full_ = target_store_.data();
    source_full_ = nullptr;
    buffers_.resize(1, target_full_);

    // libcint mallocs its working cache on every call unless handed one. Size one
    // up front as the largest any quartet needs; as in pyscf, the (gg|gg) diagonal
    // quartets bound all others. (A null `out` makes libcint return the size.)
    size_t cache_size = 0;
    for (int g = 0; g < nbas; ++g) {
        int shls[4] = {g, g, g, g};
        const size_t sz = cart_ ? int2e_cart(nullptr, nullptr, shls, atm_.data(), natm, bas_.data(), nbas,
                                             env_.data(), nullptr, nullptr)
                                : int2e_sph(nullptr, nullptr, shls, atm_.data(), natm, bas_.data(), nbas,
                                            env_.data(), nullptr, nullptr);
        cache_size = std::max(cache_size, sz);
    }
    cache_.assign(cache_size, 0.0);
}

void LibcintTwoElectronInt::initialize_sieve() {
    // Manual initialization, for objects built with SCREENING NONE. Unlike
    // Libint2, there are no engine-side shell-pair data to rebuild.
    create_sieve_pair_info_manager();
    create_blocks();
}

int LibcintTwoElectronInt::append_basis(const BasisSet &bs, const std::vector<double> &coef_ratio) {
    // Atoms: each basis contributes its own molecule's centers, so bases on
    // different molecules (or a dummy basis) never index the wrong coordinates.
    const int atom_start = static_cast<int>(atm_.size() / ATM_SLOTS);
    auto mol = bs.molecule();
    for (int a = 0; a < mol->natom(); ++a) {
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

    const int start = static_cast<int>(bas_.size() / BAS_SLOTS);
    for (int s = 0; s < bs.nshell(); ++s) {
        const GaussianShell &sh = bs.shell(s);
        const int nprim = sh.nprimitive();

        const int exp_off = static_cast<int>(env_.size());
        for (int p = 0; p < nprim; ++p) env_.push_back(sh.exp(p));

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
            // Exactly the coefficients the Libint2 path uses -- normalization embedded,
            // or deliberately not (diffuse external charges and SAP potentials; see
            // BasisSet::negative_gaussian_normalization_to_coefficients) -- in
            // libcint's convention. Re-deriving the normalization from the original
            // coefficients instead would silently override such bases.
            const auto &l2coef = bs.l2_shell(s).contr[0].coeff;
            for (int p = 0; p < nprim; ++p) env_.push_back(l2coef[p] * coef_ratio[l]);
        }

        bas_pure_.push_back(sh.is_pure() ? 1 : 0);
        bas_.push_back(atom_start + bs.shell_to_center(s));  // ATOM_OF
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

    const BasisSet *bs[4] = {original_bs1_.get(), original_bs2_.get(), original_bs3_.get(), original_bs4_.get()};
    int maxam = 0;
    for (const auto *b : bs) maxam = std::max(maxam, b->max_am());
    std::vector<double> coef_ratio;
    for (int l = 0; l <= maxam; ++l) coef_ratio.push_back(libint2_to_libcint(l, cart_));

    for (int i = 0; i < 4; ++i) {
        int start = -1;
        for (const auto &kv : basis_starts_)
            if (kv.first == bs[i]) start = kv.second;
        if (start < 0) {
            start = append_basis(*bs[i], coef_ratio);
            basis_starts_.emplace_back(bs[i], start);
        }
        bas_start_[i] = start;
    }

    // libcint drops a primitive quartet when its log-magnitude estimate falls below
    // -expcutoff (default 60). The estimate approximates the polynomial factor as
    // (d+1)^(li+lj) where its own derivation (CINTset_pairdata) has
    // (d+1/sqrt(aij))^(li+lj), so for diffuse high-l functions it underestimates by
    // up to (li+lj)/2 ln(1/aij) per pair and wrongly drops significant integrals
    // (e.g., (gg|gg) = 0.01 for a g exponent of 2.6e-4; scf-auto-cholesky). Raise the
    // cutoff by that bound over the bases' smallest exponent and largest l so libcint
    // never screens what libint2 keeps; psi4's own sieve does the real screening.
    double amin = 1.0;
    int lmax = 0;
    for (const auto *b : bs) {
        if (is_dummy_basis(*b)) continue;
        lmax = std::max(lmax, b->max_am());
        for (int s = 0; s < b->nshell(); ++s)
            for (int p = 0; p < b->shell(s).nprimitive(); ++p) amin = std::min(amin, b->shell(s).exp(p));
    }
    env_[PTR_EXPCUTOFF] = EXPCUTOFF_DEFAULT + 4 * lmax * 0.5 * std::log(1.0 / amin);
}

size_t LibcintTwoElectronInt::compute_quartet(int g1, int g2, int g3, int g4, int l1, int l2, int l3, int l4) {
    const int n1 = cart_ ? ncart(l1) : 2 * l1 + 1;
    const int n2 = cart_ ? ncart(l2) : 2 * l2 + 1;
    const int n3 = cart_ ? ncart(l3) : 2 * l3 + 1;
    const int n4 = cart_ ? ncart(l4) : 2 * l4 + 1;
    const size_t n = static_cast<size_t>(n1) * n2 * n3 * n4;
    if (cint_buf_.size() < n) cint_buf_.assign(n, 0.0);

    int shls[4] = {g1, g2, g3, g4};
    const int nbas = static_cast<int>(bas_.size() / BAS_SLOTS);
    const int natm = static_cast<int>(atm_.size() / ATM_SLOTS);
    int ok = cart_ ? int2e_cart(cint_buf_.data(), nullptr, shls, atm_.data(), natm, bas_.data(), nbas, env_.data(),
                                static_cast<CINTOpt *>(opt_), cache_.data())
                   : int2e_sph(cint_buf_.data(), nullptr, shls, atm_.data(), natm, bas_.data(), nbas, env_.data(),
                               static_cast<CINTOpt *>(opt_), cache_.data());
    double *tgt = target_full_;
    if (!ok) {
        std::fill(tgt, tgt + n, 0.0);
        return 0;
    }

    // Repack libcint (col-major, index 1 fastest) into psi4 (row-major, index 1
    // slowest), mapping each shell's components. Spherical: libcint m-order ->
    // psi4 Gaussian-m order (all +1 signs). Cartesian: identity (libcint's
    // cartesian order already equals psi4's CartesianIter order).
    auto slot = [&](int i, int l) { return cart_ ? i : psi4_slot(libcint_m(i, l)); };
    for (int a = 0; a < n1; ++a) {
        const int pa = slot(a, l1);
        for (int b = 0; b < n2; ++b) {
            const int pb = slot(b, l2);
            for (int c = 0; c < n3; ++c) {
                const int pc = slot(c, l3);
                for (int d = 0; d < n4; ++d) {
                    const int pd = slot(d, l4);
                    const size_t iC =
                        static_cast<size_t>(a) + n1 * (static_cast<size_t>(b) + n2 * (static_cast<size_t>(c) + static_cast<size_t>(n3) * d));
                    const size_t iP = (((static_cast<size_t>(pa) * n2 + pb) * n3 + pc) * n4 + pd);
                    tgt[iP] = cint_buf_[iC];
                }
            }
        }
    }
    if (mixed_) {
        const int g[4] = {g1, g2, g3, g4};
        const int l[4] = {l1, l2, l3, l4};
        return transform_mixed(g, l);
    }
    return n;
}

size_t LibcintTwoElectronInt::transform_mixed(const int *g, const int *l) {
    // target_full_ holds a row-major cartesian quartet. Transform, one index at a
    // time, each shell psi4 wants spherical, using psi4's solid harmonics (the
    // same cartesian -> Gaussian-ordered spherical transform the Libint2 path
    // effectively applies, given libcint cartesians normalized like libint2's).
    size_t d[4];
    for (int k = 0; k < 4; ++k) d[k] = ncart(l[k]);
    double *data = target_full_;
    for (int k = 0; k < 4; ++k) {
        if (!bas_pure_[g[k]] || l[k] == 0) continue;
        const size_t np = 2 * l[k] + 1;
        size_t pre = 1, post = 1;
        for (int j = 0; j < k; ++j) pre *= d[j];
        for (int j = k + 1; j < 4; ++j) post *= d[j];
        double *out = mixed_buf_.data();
        std::fill(out, out + pre * np * post, 0.0);
        const SphericalTransform &st = sph_trans_[l[k]];
        for (int c = 0; c < st.n(); ++c) {
            const size_t ci = st.cartindex(c), pi = st.pureindex(c);
            const double coef = st.coef(c);
            for (size_t i = 0; i < pre; ++i) {
                const double *src = data + (i * d[k] + ci) * post;
                double *dst = out + (i * np + pi) * post;
                for (size_t j = 0; j < post; ++j) dst[j] += coef * src[j];
            }
        }
        d[k] = np;
        std::copy(out, out + pre * np * post, data);
    }
    return d[0] * d[1] * d[2] * d[3];
}

size_t LibcintTwoElectronInt::compute_shell(const AOShellCombinationsIterator &shellIter) {
    return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
}

size_t LibcintTwoElectronInt::compute_shell(int s1, int s2, int s3, int s4) {
    const int l1 = original_bs1_->shell(s1).am(), l2 = original_bs2_->shell(s2).am();
    const int l3 = original_bs3_->shell(s3).am(), l4 = original_bs4_->shell(s4).am();
    curr_buff_size_ = original_bs1_->shell(s1).nfunction() * original_bs2_->shell(s2).nfunction() *
                      original_bs3_->shell(s3).nfunction() * original_bs4_->shell(s4).nfunction();

    const size_t ret = compute_quartet(bas_start_[0] + s1, bas_start_[1] + s2, bas_start_[2] + s3, bas_start_[3] + s4,
                                       l1, l2, l3, l4);
    buffers_[0] = target_full_;
    return ret;
}

size_t LibcintTwoElectronInt::compute_shell_for_sieve(const std::shared_ptr<BasisSet> bs, int s1, int s2, int s3,
                                                      int s4, bool /*is_bra*/) {
    int start = -1;
    for (const auto &kv : basis_starts_)
        if (kv.first == bs.get()) start = kv.second;
    if (start < 0) throw PSIEXCEPTION("LIBCINT backend: sieve basis not found in libcint environment.");

    const int l1 = bs->shell(s1).am(), l2 = bs->shell(s2).am(), l3 = bs->shell(s3).am(), l4 = bs->shell(s4).am();
    curr_buff_size_ = bs->shell(s1).nfunction() * bs->shell(s2).nfunction() * bs->shell(s3).nfunction() *
                      bs->shell(s4).nfunction();

    const size_t ret = compute_quartet(start + s1, start + s2, start + s3, start + s4, l1, l2, l3, l4);
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
