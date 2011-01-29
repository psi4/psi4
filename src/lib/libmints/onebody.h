#ifndef _psi_src_lib_libmints_onebody_h_
#define _psi_src_lib_libmints_onebody_h_

#include <vector>
#include <boost/shared_ptr.hpp>

namespace psi {

class Matrix;
class SimpleMatrix;
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

    double *buffer_;
    double *target_;
    double *tformbuf_;

    unsigned int count_;
    int deriv_;
    int natom_;

    int nchunk_;

    OneBodyAOInt(std::vector<SphericalTransform>&, boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2, int deriv=0);
    virtual void compute_pair(const boost::shared_ptr<GaussianShell>& s1, const boost::shared_ptr<GaussianShell>& s2) = 0;

    void set_chunks(int nchunk) { nchunk_ = nchunk; }
    void pure_transform(const boost::shared_ptr<GaussianShell>&, const boost::shared_ptr<GaussianShell>&, int=1);

    /// Normalize Cartesian functions based on angular momentum
    void normalize_am(const boost::shared_ptr<GaussianShell>&, const boost::shared_ptr<GaussianShell>&, int nchunk=1);

public:
    virtual ~OneBodyAOInt();

    /// Basis set on center one.
    boost::shared_ptr<BasisSet> basis();
    /// Basis set on center one.
    boost::shared_ptr<BasisSet> basis1();
    /// Basis set on center two.
    boost::shared_ptr<BasisSet> basis2();

    /// Buffer where the integrals are placed.
    const double *buffer() const;

    /// Compute the integrals between basis function in the given shell pair.
    void compute_shell(int, int);

    /*! @{
     * Computes all integrals and stores them in result
     * @param result Shared matrix object that will hold the results.
     */
    void compute(boost::shared_ptr<SimpleMatrix> result);
    /*! @} */

    /*! @{
     * Computes all integrals and stores them in result
     * @param result Shared matrix object that will hold the results.
     */
    void compute(boost::shared_ptr<Matrix> result);
    /*! @} */

    /// Computes all integrals and stores them in result by default this method throws
    virtual void compute(std::vector<boost::shared_ptr<Matrix> > &result);
    /// Computes all integrals and stores them in result by default this method throws
    virtual void compute(std::vector<boost::shared_ptr<SimpleMatrix> > &result);

    /// Does the method provide first derivatives?
    virtual bool has_deriv1() { return false; }

    /// Does the method provide second derivatives?
    virtual bool has_deriv2() { return false; }

    /// Computes the first derivatives and stores them in result
    virtual void compute_deriv1(std::vector<boost::shared_ptr<Matrix> > &result);
    /// Computes the first derivatives and stores them in result
    virtual void compute_deriv1(std::vector<boost::shared_ptr<SimpleMatrix> > &result);
    /// Computes the second derivatives and stores them in result
    virtual void compute_deriv2(std::vector<boost::shared_ptr<SimpleMatrix> > &result);

    /// Computes the integrals between basis function in the given shell pair
    virtual void compute_shell_deriv1(int, int);
    /// Computes the integrals between basis function in the given shell pair
    virtual void compute_shell_deriv2(int, int);

    /// Return true if the clone member can be called. By default returns false.
    virtual bool cloneable();

    /// Returns a clone of this object. By default throws an exception.
    virtual OneBodyAOInt* clone();
};

}

#endif
