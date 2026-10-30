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

#include "psi4/libmints/libcintonebody.h"

#include <algorithm>

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/gshell.h"
#include "psi4/libmints/libcint_util.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libpsi4util/exception.h"

extern "C" {
#include <cint_funcs.h>
}

namespace psi {

using libcint_util::ncart;

namespace {

/// Index of (u, v), u <= v, in the row-major upper triangle of an n x n matrix: the
/// order of Libint2's second-derivative results.
inline size_t upper_triangle(int u, int v, int n) { return static_cast<size_t>(u) * (2 * n - u - 1) / 2 + v; }

/// Number of psi4 functions in a shell pair, the stride between components in store_.
inline size_t pair_size(const LibcintOneBodyEngine &engine, int s1, int s2) {
    return static_cast<size_t>(engine.nfunction1(s1)) * engine.nfunction2(s2);
}

/// Size store_ and point buffers_ at its first nchunk blocks of the largest shell pair.
void point_buffers(const BasisSet &bs1, const BasisSet &bs2, size_t nchunk, std::vector<double> &store,
                   std::vector<const double *> &buffers) {
    const size_t block = static_cast<size_t>(bs1.max_function_per_shell()) * bs2.max_function_per_shell();
    store.assign(nchunk * block, 0.0);
    buffers.resize(nchunk);
}

/// Overlap and kinetic integrals depend on the bra and ket centers A and B only
/// through A - B, so d/dB = -d/dA. d/dA moves the bra function, -nabla of the
/// electron's coordinate in libcint's first slot. Fill Libint2's six first
/// derivatives (A x, y, z, then B x, y, z).
void two_center_deriv1(LibcintOneBodyEngine &engine, LibcintOneBodyEngine::Kernel ip, int s1, int s2,
                       std::vector<double> &store, std::vector<const double *> &buffers) {
    const size_t n = pair_size(engine, s1, s2);
    double *dA[3] = {store.data(), store.data() + n, store.data() + 2 * n};
    engine.compute(ip, s1, s2, false, -1.0, false, dA);
    for (size_t k = 0; k < 3 * n; ++k) store[3 * n + k] = -store[k];
    buffers.resize(6);
    for (int c = 0; c < 6; ++c) buffers[c] = store.data() + c * n;
}

/// Second derivatives of overlap or kinetic integrals: d2/dAi dAj is nabla_i nabla_j
/// on the bra function (component 3i+j of the ipip kernel), and with d/dB = -d/dA,
/// d2/dAi dBj = -d2/dAi dAj and d2/dBi dBj = d2/dAi dAj. Fill Libint2's 21: the
/// upper triangle of the 6 x 6 matrix of the coordinates (A x, y, z, B x, y, z).
void two_center_deriv2(LibcintOneBodyEngine &engine, LibcintOneBodyEngine::Kernel ipip, int s1, int s2,
                       std::vector<double> &store, std::vector<const double *> &buffers) {
    const size_t n = pair_size(engine, s1, s2);
    double *aa = store.data() + 21 * n;  // the nine d2/dAi dAj, after the 21 results
    double *out[9];
    for (int c = 0; c < 9; ++c) out[c] = aa + c * n;
    engine.compute(ipip, s1, s2, false, 1.0, false, out);
    for (int u = 0; u < 6; ++u)
        for (int v = u; v < 6; ++v) {
            const double *src = aa + (3 * (u % 3) + v % 3) * n;
            const double sign = (u < 3) == (v < 3) ? 1.0 : -1.0;
            double *dst = store.data() + upper_triangle(u, v, 6) * n;
            for (size_t k = 0; k < n; ++k) dst[k] = sign * src[k];
        }
    buffers.resize(21);
    for (int c = 0; c < 21; ++c) buffers[c] = store.data() + c * n;
}

}  // namespace

// ==> LibcintOneBodyEngine <== //

LibcintOneBodyEngine::LibcintOneBodyEngine(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2)
    : bs1_(bs1), bs2_(bs2) {
    bool any_cart = false, any_pure = false;
    for (const auto *bs : {bs1_.get(), bs2_.get()})
        for (int s = 0; s < bs->nshell(); ++s) {
            if (bs->shell(s).am() == 0) continue;
            if (bs->shell(s).is_cartesian())
                any_cart = true;
            else
                any_pure = true;
        }
    cart_ = any_cart;
    mixed_ = any_cart && any_pure;
    if (mixed_)
        for (int l = 0; l <= std::max(bs1_->max_am(), bs2_->max_am()); ++l) sph_trans_.emplace_back(l);
    build_environment();
}

void LibcintOneBodyEngine::build_environment() {
    atm_.clear();
    bas_.clear();
    bas_pure_.clear();
    env_.assign(PTR_ENV_START, 0.0);

    for (const auto &q : charges_) {
        const int coord_off = static_cast<int>(env_.size());
        env_.insert(env_.end(), q.second.begin(), q.second.end());
        env_.push_back(q.first);
        atm_.push_back(0);                // CHARGE_OF, unused under FRAC_CHARGE_NUC
        atm_.push_back(coord_off);        // PTR_COORD
        atm_.push_back(FRAC_CHARGE_NUC);  // NUC_MOD_OF
        atm_.push_back(0);                // PTR_ZETA
        atm_.push_back(coord_off + 3);    // PTR_FRAC_CHARGE
        atm_.push_back(0);
    }

    std::vector<double> coef_ratio;
    for (int l = 0; l <= std::max(bs1_->max_am(), bs2_->max_am()); ++l)
        coef_ratio.push_back(libcint_util::libint2_to_libcint(l, cart_));
    bas_start1_ = libcint_util::append_basis(*bs1_, coef_ratio, atm_, bas_, env_, bas_pure_);
    bas_start2_ =
        bs2_ == bs1_ ? bas_start1_ : libcint_util::append_basis(*bs2_, coef_ratio, atm_, bas_, env_, bas_pure_);
}

void LibcintOneBodyEngine::set_charges(const ChargeField &Zxyz) {
    charges_ = Zxyz;
    build_environment();
}

void LibcintOneBodyEngine::set_rinv_origin(const std::array<double, 3> &R) {
    std::copy(R.begin(), R.end(), env_.begin() + PTR_RINV_ORIG);
}

int LibcintOneBodyEngine::nfunction(int g) const {
    const int l = bas_[g * BAS_SLOTS + ANG_OF];
    return cart_ && !(mixed_ && bas_pure_[g]) ? ncart(l) : 2 * l + 1;
}

int LibcintOneBodyEngine::nfunction1(int s1) const { return nfunction(bas_start1_ + s1); }
int LibcintOneBodyEngine::nfunction2(int s2) const { return nfunction(bas_start2_ + s2); }

int LibcintOneBodyEngine::compute(Kernel k, int s1, int s2, bool swap, double scale, bool accumulate,
                                  double *const *out) {
    CINTIntegralFunction *f = nullptr;
    int ncomp = 1;
    int natm = static_cast<int>(atm_.size() / ATM_SLOTS);
    switch (k) {
        case Kernel::OVLP:
            f = cart_ ? int1e_ovlp_cart : int1e_ovlp_sph;
            break;
        case Kernel::KIN:
            f = cart_ ? int1e_kin_cart : int1e_kin_sph;
            break;
        case Kernel::NUC:
            f = cart_ ? int1e_nuc_cart : int1e_nuc_sph;
            natm = static_cast<int>(charges_.size());
            break;
        case Kernel::IPOVLP:
            f = cart_ ? int1e_ipovlp_cart : int1e_ipovlp_sph;
            ncomp = 3;
            break;
        case Kernel::IPKIN:
            f = cart_ ? int1e_ipkin_cart : int1e_ipkin_sph;
            ncomp = 3;
            break;
        case Kernel::IPIPOVLP:
            f = cart_ ? int1e_ipipovlp_cart : int1e_ipipovlp_sph;
            ncomp = 9;
            break;
        case Kernel::IPIPKIN:
            f = cart_ ? int1e_ipipkin_cart : int1e_ipipkin_sph;
            ncomp = 9;
            break;
        case Kernel::IPRINV:
            f = cart_ ? int1e_iprinv_cart : int1e_iprinv_sph;
            ncomp = 3;
            break;
        case Kernel::IPIPRINV:
            f = cart_ ? int1e_ipiprinv_cart : int1e_ipiprinv_sph;
            ncomp = 9;
            break;
        case Kernel::IPRINVIP:
            f = cart_ ? int1e_iprinvip_cart : int1e_iprinvip_sph;
            ncomp = 9;
            break;
    }

    const int g1 = bas_start1_ + s1, g2 = bas_start2_ + s2;
    const int l1 = bas_[g1 * BAS_SLOTS + ANG_OF], l2 = bas_[g2 * BAS_SLOTS + ANG_OF];
    // Function counts libcint returns, and psi4's after any transform.
    const int n1 = cart_ ? ncart(l1) : 2 * l1 + 1, n2 = cart_ ? ncart(l2) : 2 * l2 + 1;
    const int m1 = nfunction(g1), m2 = nfunction(g2);
    const size_t nint = static_cast<size_t>(n1) * n2;
    if (cint_buf_.size() < ncomp * nint) cint_buf_.resize(ncomp * nint);
    if (pair_.size() < nint) pair_.resize(nint);
    if (scratch_.size() < nint) scratch_.resize(nint);

    int shls[2] = {swap ? g2 : g1, swap ? g1 : g2};
    const int nbas = static_cast<int>(bas_.size() / BAS_SLOTS);
    // libcint zeroes the output of a pair it screens out.
    f(cint_buf_.data(), nullptr, shls, atm_.data(), natm, bas_.data(), nbas, env_.data(), nullptr, nullptr);

    for (int c = 0; c < ncomp; ++c) {
        // Repack libcint (col-major, slot 0 fastest) into psi4 (row-major, bra slowest),
        // mapping each shell's components.
        const double *src = cint_buf_.data() + c * nint;
        for (int i2 = 0; i2 < n2; ++i2) {
            const size_t p2 = libcint_util::psi4_index(i2, l2, cart_);
            for (int i1 = 0; i1 < n1; ++i1) {
                const size_t p1 = libcint_util::psi4_index(i1, l1, cart_);
                pair_[p1 * n2 + p2] = swap ? src[i2 + static_cast<size_t>(n2) * i1] : src[i1 + static_cast<size_t>(n1) * i2];
            }
        }

        // Transform the shells psi4 wants spherical, one index at a time.
        double *in = pair_.data();
        if (mixed_) {
            double *tmp = scratch_.data();
            int d1 = n1, d2 = n2;
            if (bas_pure_[g1] && l1 > 0) {
                std::fill(tmp, tmp + static_cast<size_t>(2 * l1 + 1) * d2, 0.0);
                const SphericalTransform &st = sph_trans_[l1];
                for (int t = 0; t < st.n(); ++t)
                    for (int j = 0; j < d2; ++j)
                        tmp[st.pureindex(t) * d2 + j] += st.coef(t) * in[st.cartindex(t) * d2 + j];
                d1 = 2 * l1 + 1;
                std::swap(in, tmp);
            }
            if (bas_pure_[g2] && l2 > 0) {
                std::fill(tmp, tmp + static_cast<size_t>(d1) * (2 * l2 + 1), 0.0);
                const SphericalTransform &st = sph_trans_[l2];
                for (int i = 0; i < d1; ++i)
                    for (int t = 0; t < st.n(); ++t)
                        tmp[i * (2 * l2 + 1) + st.pureindex(t)] += st.coef(t) * in[i * d2 + st.cartindex(t)];
                std::swap(in, tmp);
            }
        }

        const size_t m = static_cast<size_t>(m1) * m2;
        double *dst = out[c];
        if (accumulate)
            for (size_t i = 0; i < m; ++i) dst[i] += scale * in[i];
        else
            for (size_t i = 0; i < m; ++i) dst[i] = scale * in[i];
    }
    return ncomp;
}

// ==> LibcintOverlapInt <== //

LibcintOverlapInt::LibcintOverlapInt(std::vector<SphericalTransform> &st, std::shared_ptr<BasisSet> bs1,
                                     std::shared_ptr<BasisSet> bs2, int deriv)
    : OverlapInt(st, bs1, bs2, deriv), engine_(bs1, bs2) {
    // 21 second derivatives, and the nine that make them up
    point_buffers(*bs1, *bs2, deriv == 2 ? 30 : (deriv == 1 ? 6 : 1), store_, buffers_);
}

void LibcintOverlapInt::compute_shell(int s1, int s2) {
    double *out = store_.data();
    engine_.compute(LibcintOneBodyEngine::Kernel::OVLP, s1, s2, false, 1.0, false, &out);
    buffers_.resize(1);
    buffers_[0] = store_.data();
}

void LibcintOverlapInt::compute_shell_deriv1(int s1, int s2) {
    two_center_deriv1(engine_, LibcintOneBodyEngine::Kernel::IPOVLP, s1, s2, store_, buffers_);
    set_chunks(6);
}

void LibcintOverlapInt::compute_shell_deriv2(int s1, int s2) {
    two_center_deriv2(engine_, LibcintOneBodyEngine::Kernel::IPIPOVLP, s1, s2, store_, buffers_);
    set_chunks(21);
}

// ==> LibcintKineticInt <== //

LibcintKineticInt::LibcintKineticInt(std::vector<SphericalTransform> &st, std::shared_ptr<BasisSet> bs1,
                                     std::shared_ptr<BasisSet> bs2, int deriv)
    : KineticInt(st, bs1, bs2, deriv), engine_(bs1, bs2) {
    point_buffers(*bs1, *bs2, deriv == 2 ? 30 : (deriv == 1 ? 6 : 1), store_, buffers_);
}

void LibcintKineticInt::compute_shell(int s1, int s2) {
    double *out = store_.data();
    engine_.compute(LibcintOneBodyEngine::Kernel::KIN, s1, s2, false, 1.0, false, &out);
    buffers_.resize(1);
    buffers_[0] = store_.data();
}

void LibcintKineticInt::compute_shell_deriv1(int s1, int s2) {
    two_center_deriv1(engine_, LibcintOneBodyEngine::Kernel::IPKIN, s1, s2, store_, buffers_);
    set_chunks(6);
}

void LibcintKineticInt::compute_shell_deriv2(int s1, int s2) {
    two_center_deriv2(engine_, LibcintOneBodyEngine::Kernel::IPIPKIN, s1, s2, store_, buffers_);
    set_chunks(21);
}

// ==> LibcintPotentialInt <== //

LibcintPotentialInt::LibcintPotentialInt(std::vector<SphericalTransform> &st, std::shared_ptr<BasisSet> bs1,
                                         std::shared_ptr<BasisSet> bs2, int deriv)
    : PotentialInt(st, bs1, bs2, deriv), engine_(bs1, bs2) {
    // As PotentialInt: the nuclei of the bra basis's molecule.
    LibcintOneBodyEngine::ChargeField Zxyz;
    auto mol = bs1->molecule();
    for (int A = 0; A < mol->natom(); ++A) Zxyz.push_back({mol->Z(A), {mol->x(A), mol->y(A), mol->z(A)}});
    set_charge_field(Zxyz);
}

void LibcintPotentialInt::set_charge_field(const std::vector<std::pair<double, std::array<double, 3>>> &Zxyz) {
    PotentialInt::set_charge_field(Zxyz);
    engine_.set_charges(Zxyz);
    const size_t ncoord = 3 * (2 + Zxyz.size());
    // Second derivatives, and per charge the 27 that make them up
    const size_t nchunk = deriv_ == 2 ? ncoord * (ncoord + 1) / 2 + 27 : (deriv_ == 1 ? ncoord + 3 : 1);
    point_buffers(*bs1_, *bs2_, nchunk, store_, buffers_);
}

void LibcintPotentialInt::compute_shell(int s1, int s2) {
    double *out = store_.data();
    engine_.compute(LibcintOneBodyEngine::Kernel::NUC, s1, s2, false, 1.0, false, &out);
    buffers_.resize(1);
    buffers_[0] = store_.data();
}

void LibcintPotentialInt::compute_shell_deriv1(int s1, int s2) {
    // The potential of charge k, V_k = -q_k/|r - C_k|, depends on the three centers
    // only through their differences, so d/dC_k = -(d/dA + d/dB) for its term. With
    // ip = <nabla a|1/|r - C_k||b>, d/dA V_k = q_k ip and d/dB V_k = q_k ip of the
    // swapped pair.
    using K = LibcintOneBodyEngine::Kernel;
    const auto &charges = engine_.charges();
    const size_t ncoord = 3 * (2 + charges.size());
    const size_t n = pair_size(engine_, s1, s2);
    std::fill(store_.begin(), store_.begin() + (ncoord + 3) * n, 0.0);
    auto block = [&](size_t coord) { return store_.data() + coord * n; };

    double *dA[3] = {block(0), block(1), block(2)};
    double *dB[3] = {block(3), block(4), block(5)};
    double *ip[3] = {block(ncoord), block(ncoord + 1), block(ncoord + 2)};  // scratch
    for (size_t k = 0; k < charges.size(); ++k) {
        const double q = charges[k].first;
        if (q == 0.0) continue;
        engine_.set_rinv_origin(charges[k].second);
        double *dC[3] = {block(6 + 3 * k), block(7 + 3 * k), block(8 + 3 * k)};
        engine_.compute(K::IPRINV, s1, s2, false, q, false, ip);
        for (int x = 0; x < 3; ++x)
            for (size_t i = 0; i < n; ++i) {
                dA[x][i] += ip[x][i];
                dC[x][i] = -ip[x][i];
            }
        engine_.compute(K::IPRINV, s1, s2, true, q, false, ip);
        for (int x = 0; x < 3; ++x)
            for (size_t i = 0; i < n; ++i) {
                dB[x][i] += ip[x][i];
                dC[x][i] -= ip[x][i];
            }
    }
    buffers_.resize(ncoord);
    for (size_t c = 0; c < ncoord; ++c) buffers_[c] = block(c);
    set_chunks(ncoord);
}

void LibcintPotentialInt::compute_shell_deriv2(int s1, int s2) {
    // Per charge k, the second derivatives of its term with respect to the bra
    // and ket centers A and B are, with V_k = -q_k/|r - C_k|,
    //   AA = <nabla nabla a|V_k|b>, AB = <nabla a|V_k|nabla b>, BB = <a|V_k|nabla nabla b>,
    // and with d/dC_k = -(d/dA + d/dB) on its term,
    //   AC = -(AA + AB), BC = -(BA + BB), CC = AA + AB + BA + BB,
    // where BA_ij = AB_ji. Different charges don't mix.
    using K = LibcintOneBodyEngine::Kernel;
    const auto &charges = engine_.charges();
    const int ncoord = static_cast<int>(3 * (2 + charges.size()));
    const size_t nresult = static_cast<size_t>(ncoord) * (ncoord + 1) / 2;
    const size_t n = pair_size(engine_, s1, s2);
    std::fill(store_.begin(), store_.begin() + (nresult + 27) * n, 0.0);
    auto block = [&](int u, int v) {
        return store_.data() + (u <= v ? upper_triangle(u, v, ncoord) : upper_triangle(v, u, ncoord)) * n;
    };
    double *aa = store_.data() + nresult * n, *ab = aa + 9 * n, *bb = ab + 9 * n;  // scratch
    double *out_aa[9], *out_ab[9], *out_bb[9];
    for (int c = 0; c < 9; ++c) {
        out_aa[c] = aa + c * n;
        out_ab[c] = ab + c * n;
        out_bb[c] = bb + c * n;
    }

    for (size_t k = 0; k < charges.size(); ++k) {
        const double q = charges[k].first;
        if (q == 0.0) continue;
        engine_.set_rinv_origin(charges[k].second);
        engine_.compute(K::IPIPRINV, s1, s2, false, -q, false, out_aa);
        engine_.compute(K::IPRINVIP, s1, s2, false, -q, false, out_ab);
        engine_.compute(K::IPIPRINV, s1, s2, true, -q, false, out_bb);
        const int C = 6 + 3 * static_cast<int>(k);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                const double *AA = out_aa[3 * i + j], *AB = out_ab[3 * i + j], *BA = out_ab[3 * j + i],
                             *BB = out_bb[3 * i + j];
                // Each unordered pair of coordinates once: within a center, i <= j.
                const bool diag = i <= j;
                double *dAA = diag ? block(i, j) : nullptr, *dBB = diag ? block(3 + i, 3 + j) : nullptr,
                       *dCC = diag ? block(C + i, C + j) : nullptr;
                double *dAB = block(i, 3 + j), *dAC = block(i, C + j), *dBC = block(3 + i, C + j);
                for (size_t p = 0; p < n; ++p) {
                    if (diag) {
                        dAA[p] += AA[p];
                        dBB[p] += BB[p];
                        dCC[p] = AA[p] + AB[p] + BA[p] + BB[p];
                    }
                    dAB[p] += AB[p];
                    dAC[p] = -(AA[p] + AB[p]);
                    dBC[p] = -(BA[p] + BB[p]);
                }
            }
    }
    buffers_.resize(nresult);
    for (size_t c = 0; c < nresult; ++c) buffers_[c] = store_.data() + c * n;
    set_chunks(static_cast<int>(nresult));
}

}  // namespace psi
