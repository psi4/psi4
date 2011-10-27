#ifndef _psi_src_lib_libmints_twobody_h
#define _psi_src_lib_libmints_twobody_h

#include <boost/shared_ptr.hpp>

namespace psi {

class IntegralFactory;
class AOShellCombinationsIterator;
class BasisSet;
class GaussianShell;

/*! \ingroup MINTS
 *  \class TwoBodyInt
 *  \brief Two body integral base class.
 */
class TwoBodyAOInt
{
protected:
    const IntegralFactory* integral_;

    boost::shared_ptr<BasisSet> bs1_;
    boost::shared_ptr<BasisSet> bs2_;
    boost::shared_ptr<BasisSet> bs3_;
    boost::shared_ptr<BasisSet> bs4_;

    const boost::shared_ptr<BasisSet> original_bs1_;
    const boost::shared_ptr<BasisSet> original_bs2_;
    const boost::shared_ptr<BasisSet> original_bs3_;
    const boost::shared_ptr<BasisSet> original_bs4_;

    /// Buffer to hold the final integrals.
    double *target_;
    /// Buffer to hold the transformation intermediates.
    double *tformbuf_;
    /// Buffer to hold the initially computed integrals.
    double *source_;
    /// Maximum number of unique quartets needed to compute a set of SO's
    int max_unique_quartets_;
    /// Number of atoms.
    int natom_;
    /// Derivative level.
    int deriv_;
    /// Whether to force integrals to be generated in the Cartesian (AO) basis;
    bool force_cartesian_;
    //! The order of the derivative integral buffers, after permuting shells
    int *buffer_offsets_;
    //! The number of cartesian functions in this shell quartet
    size_t ncart_;

    void permute_target(double *s, double *t, int sh1, int sh2, int sh3, int sh4, bool p12, bool p34, bool p13p24);
    void permute_1234_to_1243(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_2134(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_2143(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_3412(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_4312(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_3421(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_4321(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);

//    TwoBodyInt(boost::shared_ptr<BasisSet> bs1,
//               boost::shared_ptr<BasisSet> bs2,
//               boost::shared_ptr<BasisSet> bs3,
//               boost::shared_ptr<BasisSet> bs4,
//               int deriv = 0);

    TwoBodyAOInt(const IntegralFactory* intsfactory, int deriv=0);

public:
    virtual ~TwoBodyAOInt();

    /// Basis set on center one
    boost::shared_ptr<BasisSet> basis();
    /// Basis set on center one
    boost::shared_ptr<BasisSet> basis1();
    /// Basis set on center two
    boost::shared_ptr<BasisSet> basis2();
    /// Basis set on center three
    boost::shared_ptr<BasisSet> basis3();
    /// Basis set on center four
    boost::shared_ptr<BasisSet> basis4();

    /// Sets whether we're forcing this object to always generate Cartesian integrals
    void set_force_cartesian(bool t_f) { force_cartesian_ = t_f; }

    /// Returns the derivative level this object is setup for.
    int deriv() const { return deriv_; }

    /// Buffer where the integrals are placed
    const double *buffer() const { return target_; }

    /// Returns the integral factory used to create this object
    const IntegralFactory* integral() const { return integral_; }

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    virtual void compute_shell(const AOShellCombinationsIterator&) = 0;

    /// Compute the integrals
    virtual void compute_shell(int, int, int, int) = 0;

    virtual double get_derivative_integral(int center, int xyz, size_t index);

    ///Is the shell zero?
    virtual int shell_is_zero(int,int,int,int) { return 0; }

    /// Compute the integrals
    virtual void compute_shell_deriv1(int, int, int, int) = 0;

    /// Normalize Cartesian functions based on angular momentum
    void normalize_am(boost::shared_ptr<GaussianShell>, boost::shared_ptr<GaussianShell>, boost::shared_ptr<GaussianShell>, boost::shared_ptr<GaussianShell>, int nchunk=1);

    /// Return true if the clone member can be called. By default returns false.
    virtual bool cloneable();

    /// Returns a clone of this object. By default throws an exception
    virtual TwoBodyAOInt* clone();

    /// Results go back to buffer_
    void pure_transform(int, int, int, int, int nchunk);
};

typedef boost::shared_ptr<TwoBodyAOInt> SharedTwoBodyAOInt;

}

#endif
