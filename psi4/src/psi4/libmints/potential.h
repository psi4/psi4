/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi_src_lib_libmints_potential_h_
#define _psi_src_lib_libmints_potential_h_

#include <vector>
#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/osrecur.h"

namespace psi {
    class BasisSet;
    class GaussianShell;
    class IntegralFactory;
    class SphericalTransform;
    class CdSalcList;

/*! \ingroup MINTS
 *  \class PotentialInt
 *  \brief Computes potential integrals.
 * Use an IntegralFactory to create this object.
 */
class PotentialInt : public OneBodyAOInt
{
    /// Computes integrals between two shell objects.
    void compute_pair(const GaussianShell&, const GaussianShell&);
    /// Computes integrals between two shell objects.
    void compute_pair_deriv1_no_charge_term(const GaussianShell&, const GaussianShell& );
    void compute_pair_deriv1(const GaussianShell&, const GaussianShell& );
    void compute_pair_deriv2(const GaussianShell&, const GaussianShell& );

protected:
    /// Recursion object that does the heavy lifting.
    ObaraSaikaTwoCenterVIRecursion* potential_recur_;

    /// Matrix of coordinates/charges of partial charges
    SharedMatrix Zxyz_;

public:
    /// Constructor. Assumes nuclear centers/charges as the potential
    PotentialInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, int deriv=0);
    virtual ~PotentialInt();

    /// Computes the first derivatives and stores them in result
    virtual void compute_deriv1(std::vector<SharedMatrix > &result);

    /// Computes the first derivatives and stores them in result
    virtual void compute_deriv1_no_charge_term(std::vector<SharedMatrix > &result);
    /// Computes the first derivatives, but neglects the derivatives on the third center.
    /// This code is used for gradients in the presence of an external potential.
    void compute_shell_deriv1_no_charge_term(int, int);

    /// Computes the second derivatives and store them in result
    virtual void compute_deriv2(std::vector<SharedMatrix>& result);

    /// Set the field of charges
    void set_charge_field(SharedMatrix Zxyz) { Zxyz_ = Zxyz; }

    /// Get the field of charges
    SharedMatrix charge_field() const { return Zxyz_; }

    /// Does the method provide first derivatives?
    bool has_deriv1() { return true; }
};

class PotentialSOInt : public OneBodySOInt
{
    int natom_;
public:
    PotentialSOInt(const std::shared_ptr<OneBodyAOInt>& , const std::shared_ptr<IntegralFactory> &);
    PotentialSOInt(const std::shared_ptr<OneBodyAOInt>& , const IntegralFactory*);

    /**
     * Computes one-electron integral derivative matrices.
     * Specifically handles CdSalc SO potential integral derivatives.
     *
     * \param result Where the integral derivatives are going.
     * \param cdsalcs The Cartesian displacement SALCs that you are interested in.
     */
    void compute_deriv1(std::vector<SharedMatrix > result,
                        const CdSalcList& cdsalcs);
};

}

#endif
