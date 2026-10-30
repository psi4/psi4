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

#ifndef _psi4_libmints_libcintonebody_h
#define _psi4_libmints_libcintonebody_h

#include <array>
#include <memory>
#include <utility>
#include <vector>

#include "psi4/libmints/integral.h"
#include "psi4/libmints/kinetic.h"
#include "psi4/libmints/overlap.h"
#include "psi4/libmints/potential.h"

namespace psi {

class BasisSet;

/*! \ingroup MINTS
 *  \class LibcintOneBodyEngine
 *  \brief libcint's description of a bra and a ket basis set, and of a field of
 *         point charges, from which the libcint one-electron integral classes
 *         compute shell pairs.
 *
 *  The basis functions are those of the Libint2 path (see libcint_util), and each
 *  shell pair is returned in psi4's layout: row-major (bra slowest), psi4's
 *  component order within each shell, cartesian or spherical as each basis set
 *  asks. When the bases mix cartesian and spherical shells, libcint computes all
 *  cartesian and the spherical shells are transformed here with psi4's solid
 *  harmonics, as in LibcintTwoElectronInt.
 *
 *  The point charges come first in libcint's atom table, ahead of the centers of
 *  the basis functions, so that int1e_nuc, which sums over the first natm atoms,
 *  sums over the charges alone. They use libcint's fractional-charge model, which
 *  takes any charge, negative ones included.
 */
class LibcintOneBodyEngine {
   public:
    /// The libcint kernels used, by the integral each computes. ip differentiates the
    /// function in libcint's first slot with respect to the electron's coordinates;
    /// rinv is 1/|r - R| about the origin set_rinv_origin sets.
    enum class Kernel { OVLP, KIN, NUC, IPOVLP, IPKIN, IPIPOVLP, IPIPKIN, IPRINV, IPIPRINV, IPRINVIP };

    using ChargeField = std::vector<std::pair<double, std::array<double, 3>>>;

    LibcintOneBodyEngine(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2);

    /// Set the point charges int1e_nuc sums over.
    void set_charges(const ChargeField &Zxyz);
    const ChargeField &charges() const { return charges_; }
    /// Set the center R of the rinv kernels.
    void set_rinv_origin(const std::array<double, 3> &R);

    /// Number of psi4 functions in bra shell s1 (of bs1) and ket shell s2 (of bs2).
    int nfunction1(int s1) const;
    int nfunction2(int s2) const;

    /// Compute each component c of kernel k for bra shell s1 and ket shell s2, into
    /// out[c] in psi4's layout, scaled by scale; add to out[c] if accumulate. If
    /// swap, the ket function takes libcint's first slot, the one ip
    /// differentiates, but out[c] stays bra-major. Returns the number of components.
    int compute(Kernel k, int s1, int s2, bool swap, double scale, bool accumulate, double *const *out);

   private:
    std::shared_ptr<BasisSet> bs1_;
    std::shared_ptr<BasisSet> bs2_;
    ChargeField charges_;

    /// libcint's atom, shell and environment arrays: the charges' atoms, then those
    /// of bs1 and of bs2 (once, if they are the same basis set).
    std::vector<int> atm_;
    std::vector<int> bas_;
    std::vector<double> env_;
    /// libcint bas index of shell 0 of bs1 and of bs2
    int bas_start1_;
    int bas_start2_;

    /// Whether libcint is asked for cartesian (int1e_*_cart) or spherical (int1e_*_sph) integrals.
    bool cart_;
    /// Whether the bases mix cartesian and spherical shells; then cart_ is true and the
    /// spherical shells are transformed here.
    bool mixed_;
    /// Per libcint bas index, whether psi4 wants the shell spherical (used when mixed_).
    std::vector<char> bas_pure_;
    /// psi4 cartesian -> spherical transforms by angular momentum (used when mixed_).
    std::vector<SphericalTransform> sph_trans_;

    /// libcint's output, and a shell pair being repacked and transformed.
    std::vector<double> cint_buf_;
    std::vector<double> pair_;
    std::vector<double> scratch_;

    /// (Re)build atm_/bas_/env_ from the charges and the two basis sets.
    void build_environment();
    int nfunction(int g) const;
};

/*! \ingroup MINTS
 *  \class LibcintOverlapInt
 *  \brief Overlap integrals and their first and second derivatives, computed
 *         with libcint, in OverlapInt's layout.
 */
class LibcintOverlapInt : public OverlapInt {
    LibcintOneBodyEngine engine_;
    std::vector<double> store_;

   public:
    LibcintOverlapInt(std::vector<SphericalTransform> &, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
                      int deriv = 0);

    void compute_shell(int, int) override;
    void compute_shell_deriv1(int, int) override;
    void compute_shell_deriv2(int, int) override;
};

/*! \ingroup MINTS
 *  \class LibcintKineticInt
 *  \brief Kinetic energy integrals and their first and second derivatives,
 *         computed with libcint, in KineticInt's layout.
 */
class LibcintKineticInt : public KineticInt {
    LibcintOneBodyEngine engine_;
    std::vector<double> store_;

   public:
    LibcintKineticInt(std::vector<SphericalTransform> &, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
                      int deriv = 0);

    void compute_shell(int, int) override;
    void compute_shell_deriv1(int, int) override;
    void compute_shell_deriv2(int, int) override;
};

/*! \ingroup MINTS
 *  \class LibcintPotentialInt
 *  \brief Potential integrals of a field of point charges (by default, the
 *         nuclei of the bra basis's molecule) and their first and second
 *         derivatives, computed with libcint, in PotentialInt's layout: the
 *         derivatives are with respect to the bra center, the ket center, and
 *         then each charge in turn.
 */
class LibcintPotentialInt : public PotentialInt {
    LibcintOneBodyEngine engine_;
    std::vector<double> store_;

   public:
    LibcintPotentialInt(std::vector<SphericalTransform> &, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
                        int deriv = 0);

    void set_charge_field(const std::vector<std::pair<double, std::array<double, 3>>> &Zxyz) override;

    void compute_shell(int, int) override;
    void compute_shell_deriv1(int, int) override;
    void compute_shell_deriv2(int, int) override;
};

}  // namespace psi

#endif
