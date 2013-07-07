/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#ifndef _psi_src_lib_libmints_sointegral_onebody_h
#define _psi_src_lib_libmints_sointegral_onebody_h

#include "typedefs.h"

namespace psi {

class OneBodySOInt
{
protected:
    boost::shared_ptr<OneBodyAOInt> ob_;
    const IntegralFactory* integral_;
    int deriv_;

    boost::shared_ptr<SOBasisSet> b1_;
    boost::shared_ptr<SOBasisSet> b2_;

    void common_init();

public:
    OneBodySOInt(const boost::shared_ptr<OneBodyAOInt>&, 
                 const boost::shared_ptr<IntegralFactory> &);
    OneBodySOInt(const boost::shared_ptr<OneBodyAOInt>&, 
                 const IntegralFactory*);
    virtual ~OneBodySOInt();

    boost::shared_ptr<SOBasisSet> basis() const;
    boost::shared_ptr<SOBasisSet> basis1() const;
    boost::shared_ptr<SOBasisSet> basis2() const;

    /**
      * Returns the underlying AO integral engine being used.
      */
    boost::shared_ptr<OneBodyAOInt> ob() const;

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
    virtual void compute(std::vector<SharedMatrix > results);

    /**
     * Computes one-electron integral derivative matrices.
     *
     * \param result Where the integral derivatives are going.
     * \param cdsalcs The Cartesian displacement SALCs that you are interested
     *                in.
     */
    virtual void compute_deriv1(std::vector<SharedMatrix > result,
                                const CdSalcList& cdsalcs);
};


} // end namespace

#endif

