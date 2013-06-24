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

#ifndef _psi_src_lib_libmints_onebody_h_
#define _psi_src_lib_libmints_onebody_h_

#include <vector>
#include <exception.h>
#include "typedefs.h"
#include "vector3.h"
#include <boost/python/list.hpp>
#include "pybuffer.h"
#include "exception.h"

namespace psi {

class Matrix;
class IntegralFactory;
class BasisSet;
class GaussianShell;
class SphericalTransform;

/*! \ingroup MINTS
 *  \class OneBodyInt
 *  \brief Basis class for all one-electron integrals.
 */
class OneBodyAOInt
{
protected:
    boost::shared_ptr<BasisSet> bs1_;
    boost::shared_ptr<BasisSet> bs2_;
    std::vector<SphericalTransform>& spherical_transforms_;

    Vector3 origin_;

    double *buffer_;
    double *target_;
    double *tformbuf_;

    /// Whether we want to always generate Cartesian integrals;
    bool force_cartesian_;

    unsigned int count_;
    int deriv_;
    int natom_;

    int nchunk_;

    int buffer_size_;

    /// Whether or not to use the PyBuffer
    bool enable_pybuffer_;

    /// The PyBuffer object used for sharing the target_ buffer without copying data
    PyBuffer<double> pybuffer_;

    OneBodyAOInt(std::vector<SphericalTransform>&, boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2, int deriv=0);

    virtual void compute_pair(const GaussianShell& s1, const GaussianShell& s2) = 0;
    virtual void compute_pair_deriv1(const GaussianShell& s1, const GaussianShell& s2);
    virtual void compute_pair_deriv2(const GaussianShell& s1, const GaussianShell& s2);

    void set_chunks(int nchunk) { nchunk_ = nchunk; }
    void pure_transform(const GaussianShell&, const GaussianShell&, int=1);

    /// Normalize Cartesian functions based on angular momentum
    void normalize_am(const GaussianShell&, const GaussianShell&, int nchunk=1);

public:
    virtual ~OneBodyAOInt();

    /// Basis set on center one.
    boost::shared_ptr<BasisSet> basis();
    /// Basis set on center one.
    boost::shared_ptr<BasisSet> basis1();
    /// Basis set on center two.
    boost::shared_ptr<BasisSet> basis2();

    /// Number of chunks. Normally 1, but dipoles (3) quadrupoles (6).
    int nchunk() const { return nchunk_; }

    /// Sets whether we're forcing this object to always generate Cartesian integrals
    void set_force_cartesian(bool t_f) { force_cartesian_ = t_f; }

    /// Buffer where the integrals are placed.
    const double *buffer() const;

    /// Get a python list version of the current buffer
    const boost::python::list py_buffer() const {
        boost::python::list ret_val;
        //TODO: Just change over the pointer rather than copying, or at least do a memcopy
        for(int i = 0; i < buffer_size_; ++i)
            ret_val.append(buffer_[i]);
        return ret_val;
    }

    const PyBuffer<double>* py_buffer_object() const {
        if(!enable_pybuffer_) {
            throw PSIEXCEPTION("py_buffer object not enabled.  Used set_enable_pybuffer() first.");
        }
    	return &pybuffer_;
    }

    void set_enable_pybuffer(bool enable = true) {
        enable_pybuffer_ = enable;
    }

    /// Compute the integrals between basis function in the given shell pair.
    void compute_shell(int, int);

    /*! @{
     * Computes all integrals and stores them in result
     * @param result Shared matrix object that will hold the results.
     */
    void compute(SharedMatrix& result);
    /*! @} */

    /// Computes all integrals and stores them in result by default this method throws
    virtual void compute(std::vector<SharedMatrix > &result);

    /// Does the method provide first derivatives?
    virtual bool has_deriv1() { return false; }

    /// Does the method provide second derivatives?
    virtual bool has_deriv2() { return false; }

    /// What order of derivative was requested?
    int deriv() const { return deriv_; }

    /// Computes the first derivatives and stores them in result
    virtual void compute_deriv1(std::vector<SharedMatrix > &result);

    /// Computes the second derivatives and stores them in result
    virtual void compute_deriv2(std::vector<SharedMatrix > &result);

    /// Computes the integrals between basis function in the given shell pair
    virtual void compute_shell_deriv1(int, int);
    /// Computes the integrals between basis function in the given shell pair
    virtual void compute_shell_deriv2(int, int);

    /// Return true if the clone member can be called. By default returns false.
    virtual bool cloneable();

    /// Returns a clone of this object. By default throws an exception.
    virtual OneBodyAOInt* clone();

    /// Returns the origin (useful for properties)
    Vector3 origin() const { return origin_; }

    /// Set the origin (useful for properties)
    void set_origin(const Vector3& _origin) { origin_ = _origin; }
};

typedef boost::shared_ptr<OneBodyAOInt> SharedOneBodyAOInt;
}

#endif
