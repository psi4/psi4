/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#ifndef _psi_src_lib_libmints_sointegral_onebody_h
#define _psi_src_lib_libmints_sointegral_onebody_h

#include "psi4/libmints/typedefs.h"

namespace psi {
class SOBasisSet;
class CdSalcList;
class OneBodyAOInt;
class IntegralFactory;

class PSI_API OneBodySOInt {
   protected:
    std::shared_ptr<OneBodyAOInt> ob_;
    const IntegralFactory* integral_;
    int deriv_;

    std::shared_ptr<SOBasisSet> b1_;
    std::shared_ptr<SOBasisSet> b2_;

    void common_init();

   public:
    OneBodySOInt(const std::shared_ptr<OneBodyAOInt>&, const std::shared_ptr<IntegralFactory>&);
    OneBodySOInt(const std::shared_ptr<OneBodyAOInt>&, const IntegralFactory*);
    virtual ~OneBodySOInt();

    std::shared_ptr<SOBasisSet> basis() const;
    std::shared_ptr<SOBasisSet> basis1() const;
    std::shared_ptr<SOBasisSet> basis2() const;

    /**
     * Returns the underlying AO integral engine being used.
     */
    std::shared_ptr<OneBodyAOInt> ob() const;

    /**
     * Computes a one-electron integral matrix. Only works for symmetric
     * operators (multipole operators will not work).
     *
     * \param result Where the integrals are going.
     */
    virtual void compute(SharedMatrix result);

    /**
     * Computes one-electron integral matrices. Should be able to handle
     * multipole operators
     *
     * \param results Where the integrals are going.
     */
    virtual void compute(std::vector<SharedMatrix> results);

    /**
     * Computes one-electron integral derivative matrices.
     *
     * \param result Where the integral derivatives are going.
     * \param cdsalcs The Cartesian displacement SALCs that you are interested
     *                in.
     */
    virtual void compute_deriv1(std::vector<SharedMatrix> result, const CdSalcList& cdsalcs);
};

}  // namespace psi

#endif
