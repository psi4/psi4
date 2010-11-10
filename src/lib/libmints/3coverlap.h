#ifndef _psi_src_lib_libmints_3coverlap_h
#define _psi_src_lib_libmints_3coverlap_h

#include <boost/shared_ptr.hpp>

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

    /// Buffer to hold the final integrals.
    double *buffer_;

    void compute_pair(shared_ptr<GaussianShell> s1,
                      shared_ptr<GaussianShell> s2,
                      shared_ptr<GaussianShell> s3);
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

    /// Compute the integrals of the form (a|c|b).
    virtual void compute_shell(int, int, int);

    /// Normalize Cartesian functions based on angular momentum
    void normalize_am(boost::shared_ptr<GaussianShell>&,
                      boost::shared_ptr<GaussianShell>&,
                      boost::shared_ptr<GaussianShell>&);

    /// Perform pure (spherical) transform.
    void pure_transform(boost::shared_ptr<GaussianShell>,
                        boost::shared_ptr<GaussianShell>,
                        boost::shared_ptr<GaussianShell>);
};

}

#endif
