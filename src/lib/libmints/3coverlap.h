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

#ifndef _psi_src_lib_libmints_3coverlap_h
#define _psi_src_lib_libmints_3coverlap_h

#include "pybuffer.h"

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class IntegralFactory;
class BasisSet;
class GaussianShell;
class ObaraSaikaThreeCenterRecursion;

/** \ingroup MINTS
    \class ThreeCenterOverlapInt
    \brief Three center overlap integral.
 */
class ThreeCenterOverlapInt
{
protected:
    ObaraSaikaThreeCenterRecursion overlap_recur_;

    boost::shared_ptr<BasisSet> bs1_;
    boost::shared_ptr<BasisSet> bs2_;
    boost::shared_ptr<BasisSet> bs3_;

    /// Buffer to hold the source integrals.
    double *buffer_;

    /// Buffer for spherical harmonics
    double *temp_;

    /// Vector of Sphericaltransforms
    std::vector<SphericalTransform> st_;

    /// Whether or not to activate the PyBuffer 
    bool enable_pybuffer_;

    void compute_pair(const GaussianShell& s1,
                      const GaussianShell& s2,
                      const GaussianShell& s3);
    
    /// The PyBuffer object used for sharing the target_ buffer without copying data
    PyBuffer<double> pybuffer_;

public:
    ThreeCenterOverlapInt(std::vector<SphericalTransform>&,
                          boost::shared_ptr<BasisSet> bs1,
                          boost::shared_ptr<BasisSet> bs2,
                          boost::shared_ptr<BasisSet> bs3);

    virtual ~ThreeCenterOverlapInt();

    /// Basis set on center one.
    boost::shared_ptr<BasisSet> basis();
    /// Basis set on center one.
    boost::shared_ptr<BasisSet> basis1();
    /// Basis set on center two.
    boost::shared_ptr<BasisSet> basis2();
    /// Basis set on center three.
    boost::shared_ptr<BasisSet> basis3();

    /// Buffer where the integrals are placed.
    const double *buffer() const { return buffer_; }

    const PyBuffer<double>* py_buffer_object() const {
        if(!enable_pybuffer_) {
            throw PSIEXCEPTION("py_buffer object not enabled.  Used set_enable_pybuffer() first.");
        }
    	return &pybuffer_;
    }

    void set_enable_pybuffer(bool enable = true) {
        enable_pybuffer_ = enable;
    }

    /// Compute the integrals of the form (a|c|b).
    virtual void compute_shell(int, int, int);

    /// Normalize Cartesian functions based on angular momentum
    void normalize_am(const GaussianShell&,
                      const GaussianShell&,
                      const GaussianShell&);

    /// Perform pure (spherical) transform.
    void pure_transform(const GaussianShell&,
                        const GaussianShell&,
                        const GaussianShell&);
};

}

#endif
