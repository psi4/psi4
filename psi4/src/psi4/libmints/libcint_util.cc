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

#include "psi4/libmints/libcint_util.h"

#include <cmath>

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/gshell.h"
#include "psi4/libmints/molecule.h"

#include <libint2/shell.h>

extern "C" {
#include <cint_funcs.h>
}

// psi4_slot hard-codes the Gaussian solid-harmonic ordering (psi4's default), while the transforms of mixed
// cartesian and spherical bases use psi4's SphericalTransform, which follows psi4_SHGSHELL_ORDERING. In a build
// with the Standard ordering, the spherical repack would put components in the wrong slots, so stop such a build
// until psi4_slot follows the ordering too.
#if psi4_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_GAUSSIAN
#error "LIBCINT backend: Gaussian solid-harmonic ordering only. Teach psi4_slot Standard, or ENABLE_cint=OFF."
#endif

namespace psi {
namespace libcint_util {

/// Both conventions scale a primitive's coefficient as a^((2l+3)/4), so the ratio
/// of the unit-normalized single-primitive coefficients (here at exponent 1) holds
/// for every exponent. libcint's is measured with its own overlap integral
/// (cartesian: the axial component, as Libint2 normalizes), keeping us agnostic of
/// its l=0,1 factors.
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

int append_basis(const BasisSet &bs, const std::vector<double> &coef_ratio, std::vector<int> &atm,
                 std::vector<int> &bas, std::vector<double> &env, std::vector<char> &bas_pure) {
    // Atoms: each basis contributes its own molecule's centers, so bases on
    // different molecules (or a dummy basis) never index the wrong coordinates.
    const int atom_start = static_cast<int>(atm.size() / ATM_SLOTS);
    auto mol = bs.molecule();
    for (int a = 0; a < mol->natom(); ++a) {
        const int coord_off = static_cast<int>(env.size());
        env.push_back(mol->x(a));
        env.push_back(mol->y(a));
        env.push_back(mol->z(a));
        atm.push_back(static_cast<int>(mol->true_atomic_number(a)));  // CHARGE_OF
        atm.push_back(coord_off);                                     // PTR_COORD
        atm.push_back(0);                                             // NUC_MOD_OF
        atm.push_back(0);                                             // PTR_ZETA
        atm.push_back(0);
        atm.push_back(0);
    }

    const int start = static_cast<int>(bas.size() / BAS_SLOTS);
    for (int s = 0; s < bs.nshell(); ++s) {
        const GaussianShell &sh = bs.shell(s);
        const int nprim = sh.nprimitive();

        const int exp_off = static_cast<int>(env.size());
        for (int p = 0; p < nprim; ++p) env.push_back(sh.exp(p));

        const int l = sh.am();
        const int coef_off = static_cast<int>(env.size());
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
            env.push_back(std::sqrt(4.0 * M_PI));
        } else {
            // Exactly the coefficients the Libint2 path uses -- normalization embedded,
            // or deliberately not (diffuse external charges and SAP potentials; see
            // BasisSet::negative_gaussian_normalization_to_coefficients) -- in
            // libcint's convention. Re-deriving the normalization from the original
            // coefficients instead would silently override such bases.
            const auto &l2coef = bs.l2_shell(s).contr[0].coeff;
            for (int p = 0; p < nprim; ++p) env.push_back(l2coef[p] * coef_ratio[l]);
        }

        bas_pure.push_back(sh.is_pure() ? 1 : 0);
        bas.push_back(atom_start + bs.shell_to_center(s));  // ATOM_OF
        bas.push_back(sh.am());                             // ANG_OF
        bas.push_back(nprim);                               // NPRIM_OF
        bas.push_back(1);                                   // NCTR_OF
        bas.push_back(0);                                   // KAPPA_OF
        bas.push_back(exp_off);                             // PTR_EXP
        bas.push_back(coef_off);                            // PTR_COEFF
        bas.push_back(0);                                   // RESERVE
    }
    return start;
}

}  // namespace libcint_util
}  // namespace psi
