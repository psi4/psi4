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
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/psi4-dec.h"

#include <cmath>

#include <libint2/shell.h>

extern "C" {
#include <cint_funcs.h>
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

inline int ncart(int l) { return (l + 1) * (l + 2) / 2; }

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

/// Does the libcint in use drop the integrals of diffuse primitives of high
/// angular momentum? Its estimate of the size of a pair of primitives leaves out
/// the factor a_ij^(-(l_i+l_j)/2) of its own bound on the pair,
/// (d + 1/sqrt(a_ij))^(l_i+l_j) (pi/a_ij)^(3/2), so that, for instance, the
/// one-center (gg|gg) integrals of exponents below 3e-4 come out zero (libcint
/// and qcint 6.1.3). The derivative kernels are affected the same way.
/// Integrals of normalized primitives of a common exponent a on one center scale
/// exactly as sqrt(a): compare (gg|gg) at a = 1e-4 with a = 1 under the default
/// screening. Evaluated once per process.
bool libcint_drops_diffuse_pairs() {
    static const bool drops = []() {
        const int l = 4;
        auto gggg = [](double a) {
            int atm[ATM_SLOTS] = {0};
            atm[PTR_COORD] = PTR_ENV_START;
            int bas[BAS_SLOTS] = {0};
            bas[ANG_OF] = l;
            bas[NPRIM_OF] = 1;
            bas[NCTR_OF] = 1;
            bas[PTR_EXP] = PTR_ENV_START + 3;
            bas[PTR_COEFF] = PTR_ENV_START + 4;
            std::vector<double> env(PTR_ENV_START + 5, 0.0);
            env[PTR_ENV_START + 3] = a;
            env[PTR_ENV_START + 4] = CINTgto_norm(l, a);
            int shls[4] = {0, 0, 0, 0};
            const int n = 2 * l + 1;
            std::vector<double> buf(static_cast<size_t>(n) * n * n * n);
            int2e_sph(buf.data(), nullptr, shls, atm, 1, bas, 1, env.data(), nullptr, nullptr);
            return buf[0];  // (00|00) of the m = -l function
        };
        const double a = 1e-4;
        const double expected = std::sqrt(a) * gggg(1.0);
        const bool wrong = std::fabs(gggg(a) - expected) > 1e-6 * expected;
        if (wrong)
            outfile->Printf(
                "  Note: this libcint drops integrals of diffuse primitives of high angular momentum;\n"
                "        raising its screening cutoff for basis sets that need it.\n");
        return wrong;
    }();
    return drops;
}

}  // namespace

LibcintTwoElectronInt::LibcintTwoElectronInt(const IntegralFactory *integral, int deriv, double omega,
                                             bool use_shell_pairs, bool needs_exchange)
    : TwoBodyAOInt(integral, deriv),
      opt_(nullptr),
      opt_ip1_(nullptr),
      opt_ipip1_(nullptr),
      opt_ipvip1_(nullptr),
      opt_ip1ip2_(nullptr),
      omega_(omega) {
    if (deriv_ < 0 || deriv_ > 2)
        throw PSIEXCEPTION("LIBCINT backend: only first and second integral derivatives are implemented.");

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

    // Derivatives are taken with respect to the real centers only, as in
    // Libint2's xs_xx / xx_xs / xs_xs brakets.
    const BasisSet *bs[4] = {original_bs1_.get(), original_bs2_.get(), original_bs3_.get(), original_bs4_.get()};
    for (int c = 0; c < 4; ++c)
        if (!is_dummy_basis(*bs[c])) active_.push_back(c);

    build_environment();
    common_init();
}

LibcintTwoElectronInt::LibcintTwoElectronInt(const LibcintTwoElectronInt &rhs)
    : TwoBodyAOInt(rhs),
      atm_(rhs.atm_),
      bas_(rhs.bas_),
      env_(rhs.env_),
      opt_(nullptr),
      opt_ip1_(nullptr),
      opt_ipip1_(nullptr),
      opt_ipvip1_(nullptr),
      opt_ip1ip2_(nullptr),
      basis_starts_(rhs.basis_starts_),
      active_(rhs.active_),
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
    for (void *opt : {opt_, opt_ip1_, opt_ipip1_, opt_ipvip1_, opt_ip1ip2_}) {
        if (!opt) continue;
        CINTOpt *o = static_cast<CINTOpt *>(opt);
        CINTdel_optimizer(&o);
    }
}

void LibcintTwoElectronInt::common_init() {
    // Range separation: env[PTR_RANGE_OMEGA] > 0 -> erf (long-range), < 0 -> erfc.
    env_[PTR_RANGE_OMEGA] = omega_;

    const int nbas = static_cast<int>(bas_.size() / BAS_SLOTS);
    const int natm = static_cast<int>(atm_.size() / ATM_SLOTS);
    auto make_opt = [&](CINTOptimizerFunction *f) {
        CINTOpt *o = nullptr;
        f(&o, atm_.data(), natm, bas_.data(), nbas, env_.data());
        return static_cast<void *>(o);
    };
    opt_ = make_opt(int2e_optimizer);
    if (deriv_ >= 1) opt_ip1_ = make_opt(int2e_ip1_optimizer);
    if (deriv_ >= 2) {
        opt_ipip1_ = make_opt(int2e_ipip1_optimizer);
        opt_ipvip1_ = make_opt(int2e_ipvip1_optimizer);
        opt_ip1ip2_ = make_opt(int2e_ip1ip2_optimizer);
    }

    // Size scratch for the largest quartet block we may compute. This must
    // cover not only the factory's (bs1,bs2,bs3,bs4) block but also the Schwarz
    // sieve, which computes (MN|MN) over a single basis -- so a block up to
    // (2*maxam+1)^4 for the largest-am basis (e.g. the auxiliary in DF, whose
    // (MN|MN) blocks dwarf the aux*dummy*prim*prim factory block).
    int maxam = 0;
    for (const auto *bs : {original_bs1_.get(), original_bs2_.get(), original_bs3_.get(), original_bs4_.get()})
        maxam = std::max(maxam, bs->max_am());
    const size_t nmax = static_cast<size_t>(cart_ ? ncart(maxam) : 2 * maxam + 1);
    const size_t sieve_blk = nmax * nmax * nmax * nmax;
    quartet_size_ = 1;
    for (const auto *bs : {original_bs1_.get(), original_bs2_.get(), original_bs3_.get(), original_bs4_.get()})
        quartet_size_ *= static_cast<size_t>(cart_ ? ncart(bs->max_am()) : 2 * bs->max_am() + 1);
    // Derivative kernels return 3 (first) or 9 (second) components per call.
    const size_t ncomp = deriv_ >= 2 ? 9 : (deriv_ == 1 ? 3 : 1);
    const size_t nbuf = deriv_ >= 2 ? 78 : (deriv_ == 1 ? 12 : 1);
    target_store_.assign(std::max(sieve_blk, quartet_size_), 0.0);
    if (mixed_) {
        mixed_scratch_.assign(std::max(sieve_blk, quartet_size_), 0.0);
        sph_trans_.clear();
        for (int l = 0; l <= maxam; ++l) sph_trans_.emplace_back(l);
    }
    cint_buf_.assign(std::max(sieve_blk, ncomp * quartet_size_), 0.0);
    if (deriv_ > 0) {
        deriv_store_.assign(nbuf * quartet_size_, 0.0);
        zeros_.assign(quartet_size_, 0.0);
        discard_.assign(quartet_size_, 0.0);
    }
    target_full_ = target_store_.data();
    source_full_ = nullptr;
    buffers_.assign(nbuf, target_full_);

    setup_sieve();
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

    set_expcutoff();
}

void LibcintTwoElectronInt::set_expcutoff() {
    // libcint skips a pair of primitives when its estimate of the size of the
    // pair falls below exp(-expcutoff), 60 by default; in a four-center integral
    // the estimates of the two pairs are added. If the library leaves out the
    // factor for diffuse primitives (see libcint_drops_diffuse_pairs), a
    // calculation in a basis with such functions collapses far below the
    // variational limit. With the factor of the most diffuse pair missing from
    // the estimates of both pairs, the effective cutoff is expcutoff - 2*deficit.
    // Keep it at least 40, the smallest cutoff libcint accepts. This leaves the
    // default untouched for common basis sets, such as aug-cc-pV5Z, and screens
    // less only where the estimates fall short.
    if (!libcint_drops_diffuse_pairs()) return;
    // Smallest exponent for each angular momentum. The density-fitting dummy
    // (l = 0, exponent 0) pairs with every shell, so its zero is kept.
    std::vector<double> amin;
    const int nbas = static_cast<int>(bas_.size() / BAS_SLOTS);
    for (int g = 0; g < nbas; ++g) {
        const int l = bas_[g * BAS_SLOTS + ANG_OF];
        const int nprim = bas_[g * BAS_SLOTS + NPRIM_OF];
        const int exp_off = bas_[g * BAS_SLOTS + PTR_EXP];
        if (static_cast<int>(amin.size()) <= l) amin.resize(l + 1, -1.0);
        for (int p = 0; p < nprim; ++p)
            if (amin[l] < 0.0 || env_[exp_off + p] < amin[l]) amin[l] = env_[exp_off + p];
    }
    double deficit = 0.0;
    for (size_t li = 0; li < amin.size(); ++li)
        for (size_t lj = 0; lj < amin.size(); ++lj) {
            if (amin[li] < 0.0 || amin[lj] < 0.0) continue;
            const double aij = amin[li] + amin[lj];
            if (aij > 0.0) deficit = std::max(deficit, -0.5 * static_cast<double>(li + lj) * std::log(aij));
        }
    if (40.0 + 2.0 * deficit > 60.0) env_[PTR_EXPCUTOFF] = 40.0 + 2.0 * deficit;
}

bool LibcintTwoElectronInt::run_kernel(Kernel kernel, const int g[4], const int l[4], const int perm[4], double sign,
                                       double *const *out) {
    CINTIntegralFunction *f = nullptr;
    void *opt = nullptr;
    int ncomp = 1;
    switch (kernel) {
        case Kernel::ERI:
            f = cart_ ? int2e_cart : int2e_sph;
            opt = opt_;
            break;
        case Kernel::IP1:
            f = cart_ ? int2e_ip1_cart : int2e_ip1_sph;
            opt = opt_ip1_;
            ncomp = 3;
            break;
        case Kernel::IPIP1:
            f = cart_ ? int2e_ipip1_cart : int2e_ipip1_sph;
            opt = opt_ipip1_;
            ncomp = 9;
            break;
        case Kernel::IPVIP1:
            f = cart_ ? int2e_ipvip1_cart : int2e_ipvip1_sph;
            opt = opt_ipvip1_;
            ncomp = 9;
            break;
        case Kernel::IP1IP2:
            f = cart_ ? int2e_ip1ip2_cart : int2e_ip1ip2_sph;
            opt = opt_ip1ip2_;
            ncomp = 9;
            break;
    }

    // Function counts in psi4 center order (N) and in libcint slot order (n).
    int N[4], n[4], shls[4];
    for (int c = 0; c < 4; ++c) N[c] = cart_ ? ncart(l[c]) : 2 * l[c] + 1;
    for (int t = 0; t < 4; ++t) {
        n[t] = N[perm[t]];
        shls[t] = g[perm[t]];
    }
    const size_t nint = static_cast<size_t>(N[0]) * N[1] * N[2] * N[3];
    if (cint_buf_.size() < ncomp * nint) cint_buf_.assign(ncomp * nint, 0.0);

    const int nbas = static_cast<int>(bas_.size() / BAS_SLOTS);
    const int natm = static_cast<int>(atm_.size() / ATM_SLOTS);
    const bool nonzero = f(cint_buf_.data(), nullptr, shls, atm_.data(), natm, bas_.data(), nbas, env_.data(),
                           static_cast<CINTOpt *>(opt), nullptr);
    if (!nonzero) {
        for (int c = 0; c < ncomp; ++c) std::fill(out[c], out[c] + nint, 0.0);
        return false;
    }

    // Repack libcint (col-major, slot 0 fastest; component slowest) into psi4
    // (row-major, center 0 slowest), undoing the permutation and mapping each
    // shell's components. Spherical: libcint m-order -> psi4 Gaussian-m order
    // (all +1 signs). Cartesian: identity (libcint's cartesian order already
    // equals psi4's CartesianIter order).
    auto slot = [&](int i, int lc) { return cart_ ? i : psi4_slot(libcint_m(i, lc)); };
    const size_t stride[4] = {static_cast<size_t>(N[1]) * N[2] * N[3], static_cast<size_t>(N[2]) * N[3],
                              static_cast<size_t>(N[3]), 1};
    for (int c = 0; c < ncomp; ++c) {
        const double *src = cint_buf_.data() + c * nint;
        double *dst = out[c];
        size_t iC = 0;
        for (int i3 = 0; i3 < n[3]; ++i3) {
            const size_t o3 = slot(i3, l[perm[3]]) * stride[perm[3]];
            for (int i2 = 0; i2 < n[2]; ++i2) {
                const size_t o2 = o3 + slot(i2, l[perm[2]]) * stride[perm[2]];
                for (int i1 = 0; i1 < n[1]; ++i1) {
                    const size_t o1 = o2 + slot(i1, l[perm[1]]) * stride[perm[1]];
                    for (int i0 = 0; i0 < n[0]; ++i0, ++iC)
                        dst[o1 + slot(i0, l[perm[0]]) * stride[perm[0]]] = sign * src[iC];
                }
            }
        }
        if (mixed_) transform_mixed(dst, g, l);
    }
    return true;
}

void LibcintTwoElectronInt::transform_mixed(double *data, const int g[4], const int l[4]) {
    // data holds a row-major cartesian quartet in psi4 center order. Transform, one
    // index at a time, each shell psi4 wants spherical, using psi4's solid harmonics
    // (the same cartesian -> Gaussian-ordered spherical transform the Libint2 path
    // effectively applies, given libcint cartesians normalized like libint2's).
    // Alternate between data and mixed_scratch_, which is as large.
    size_t d[4];
    for (int k = 0; k < 4; ++k) d[k] = ncart(l[k]);
    double *in = data;
    double *out = mixed_scratch_.data();
    for (int k = 0; k < 4; ++k) {
        if (!bas_pure_[g[k]] || l[k] == 0) continue;
        const size_t np = 2 * l[k] + 1;
        size_t pre = 1, post = 1;
        for (int j = 0; j < k; ++j) pre *= d[j];
        for (int j = k + 1; j < 4; ++j) post *= d[j];
        std::fill(out, out + pre * np * post, 0.0);
        const SphericalTransform &st = sph_trans_[l[k]];
        for (int c = 0; c < st.n(); ++c) {
            const size_t ci = st.cartindex(c), pi = st.pureindex(c);
            const double coef = st.coef(c);
            for (size_t i = 0; i < pre; ++i) {
                const double *src = in + (i * d[k] + ci) * post;
                double *dst = out + (i * np + pi) * post;
                for (size_t j = 0; j < post; ++j) dst[j] += coef * src[j];
            }
        }
        d[k] = np;
        std::swap(in, out);
    }
    if (in != data) std::copy(in, in + d[0] * d[1] * d[2] * d[3], data);
}

int LibcintTwoElectronInt::nfunction(int g, int l) const {
    return cart_ && !(mixed_ && bas_pure_[g]) ? ncart(l) : 2 * l + 1;
}

size_t LibcintTwoElectronInt::compute_quartet(int g1, int g2, int g3, int g4, int l1, int l2, int l3, int l4) {
    const int g[4] = {g1, g2, g3, g4};
    const int l[4] = {l1, l2, l3, l4};
    const int identity[4] = {0, 1, 2, 3};
    double *out[1] = {target_full_};
    if (!run_kernel(Kernel::ERI, g, l, identity, 1.0, out)) return 0;
    size_t n = 1;
    for (int c = 0; c < 4; ++c) n *= nfunction(g[c], l[c]);
    return n;
}

void LibcintTwoElectronInt::quartet(int s1, int s2, int s3, int s4, int g[4], int l[4]) const {
    const int s[4] = {s1, s2, s3, s4};
    const BasisSet *bs[4] = {original_bs1_.get(), original_bs2_.get(), original_bs3_.get(), original_bs4_.get()};
    for (int c = 0; c < 4; ++c) {
        g[c] = bas_start_[c] + s[c];
        l[c] = bs[c]->shell(s[c]).am();
    }
}

size_t LibcintTwoElectronInt::compute_shell(const AOShellCombinationsIterator &shellIter) {
    return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
}

size_t LibcintTwoElectronInt::compute_shell(int s1, int s2, int s3, int s4) {
    int g[4], l[4];
    quartet(s1, s2, s3, s4, g, l);
    curr_buff_size_ = original_bs1_->shell(s1).nfunction() * original_bs2_->shell(s2).nfunction() *
                      original_bs3_->shell(s3).nfunction() * original_bs4_->shell(s4).nfunction();

    const size_t ret = compute_quartet(g[0], g[1], g[2], g[3], l[0], l[1], l[2], l[3]);
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

namespace {

/// Slot orders placing center c (0-3) in libcint's first slot while keeping it
/// a valid ERI permutation: (ab|cd) = (ba|cd) = (cd|ab) = (dc|ab).
constexpr int first_slot_perm[4][4] = {{0, 1, 2, 3}, {1, 0, 2, 3}, {2, 3, 0, 1}, {3, 2, 0, 1}};

/// Position of the coordinate pair (u, v), u <= v, in the row-major upper
/// triangle of an n x n matrix (Libint2's second-derivative order).
inline int upper_triangle(int u, int v, int n) { return u * n - u * (u - 1) / 2 + (v - u); }

}  // namespace

size_t LibcintTwoElectronInt::compute_shell_deriv1(int s1, int s2, int s3, int s4) {
    if (deriv_ < 1) throw PSIEXCEPTION("LIBCINT backend: integral object was not created for first derivatives.");
    int g[4], l[4];
    quartet(s1, s2, s3, s4, g, l);
    size_t nint = 1;
    for (int c = 0; c < 4; ++c) nint *= nfunction(g[c], l[c]);
    curr_buff_size_ = nint;

    // d/dX (ab|cd) = -(nabla x ...|...) with center X permuted into slot 0.
    bool nonzero = false;
    for (size_t k = 0; k < active_.size(); ++k) {
        double *out[3];
        for (int a = 0; a < 3; ++a) out[a] = deriv_store_.data() + (3 * k + a) * quartet_size_;
        nonzero |= run_kernel(Kernel::IP1, g, l, first_slot_perm[active_[k]], -1.0, out);
    }
    for (size_t i = 0; i < 12; ++i)
        buffers_[i] = i < 3 * active_.size() ? deriv_store_.data() + i * quartet_size_ : zeros_.data();
    return nonzero ? 12 * nint : 0;
}

size_t LibcintTwoElectronInt::compute_shell_deriv2(int s1, int s2, int s3, int s4) {
    if (deriv_ < 2) throw PSIEXCEPTION("LIBCINT backend: integral object was not created for second derivatives.");
    int g[4], l[4];
    quartet(s1, s2, s3, s4, g, l);
    size_t nint = 1;
    for (int c = 0; c < 4; ++c) nint *= nfunction(g[c], l[c]);
    curr_buff_size_ = nint;

    const int ncoord = 3 * static_cast<int>(active_.size());
    auto block = [&](int u, int v) { return deriv_store_.data() + upper_triangle(u, v, ncoord) * quartet_size_; };

    // Second nuclear derivatives equal the electronic ones (two sign flips).
    // Component 3a+b of each kernel is d/dX_a d/dY_b, with X the center in
    // libcint slot 0 and Y the other differentiated one.
    bool nonzero = false;
    for (size_t kx = 0; kx < active_.size(); ++kx) {
        const int X = active_[kx];
        for (size_t ky = kx; ky < active_.size(); ++ky) {
            const int Y = active_[ky];
            Kernel kernel;
            int perm[4];
            if (X == Y) {
                // (nabla nabla x ...|...): symmetric in a, b; keep a <= b.
                kernel = Kernel::IPIP1;
                std::copy(first_slot_perm[X], first_slot_perm[X] + 4, perm);
            } else if ((X < 2) == (Y < 2)) {
                // X, Y on the same electron: (nabla x nabla y|...), i.e. (ab| or |cd).
                kernel = Kernel::IPVIP1;
                std::copy(first_slot_perm[X], first_slot_perm[X] + 4, perm);
            } else {
                // X in the bra, Y in the ket: (nabla x ..|nabla y ..).
                kernel = Kernel::IP1IP2;
                const int x_partner = 1 - X, y_partner = 5 - Y;  // the undifferentiated center beside each
                const int order[4] = {X, x_partner, Y, y_partner};
                std::copy(order, order + 4, perm);
            }
            double *out[9];
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) {
                    const int u = 3 * static_cast<int>(kx) + a, v = 3 * static_cast<int>(ky) + b;
                    out[3 * a + b] = (X == Y && a > b) ? discard_.data() : block(u, v);
                }
            nonzero |= run_kernel(kernel, g, l, perm, 1.0, out);
        }
    }
    const size_t nfilled = static_cast<size_t>(ncoord) * (ncoord + 1) / 2;
    for (size_t i = 0; i < 78; ++i)
        buffers_[i] = i < nfilled ? deriv_store_.data() + i * quartet_size_ : zeros_.data();
    return nonzero ? 78 * nint : 0;
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
