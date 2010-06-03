#ifndef _psi_src_lib_libmints_twobody_h
#define _psi_src_lib_libmints_twobody_h

#include <libmints/matrix.h>
#include <libmints/integral.h>

namespace psi {
    
class IntegralFactory;
class BasisSet;
class GaussianShell;

/*! \ingroup MINTS
 *  \class TwoBodyInt
 *  \brief Two body integral base class.
 */
class TwoBodyInt
{
protected:
    shared_ptr<BasisSet> bs1_;
    shared_ptr<BasisSet> bs2_;
    shared_ptr<BasisSet> bs3_;
    shared_ptr<BasisSet> bs4_;

    shared_ptr<BasisSet> original_bs1_;
    shared_ptr<BasisSet> original_bs2_;
    shared_ptr<BasisSet> original_bs3_;
    shared_ptr<BasisSet> original_bs4_;
    
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
    
    void permute_target(double *s, double *t, int sh1, int sh2, int sh3, int sh4, bool p12, bool p34, bool p13p24);
    void permute_1234_to_1243(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_2134(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_2143(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_3412(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_4312(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_3421(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_4321(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    
    TwoBodyInt(
               shared_ptr<BasisSet> bs1,
               shared_ptr<BasisSet> bs2,
               shared_ptr<BasisSet> bs3,
               shared_ptr<BasisSet> bs4,
               int deriv = 0);
               
public:
    virtual ~TwoBodyInt();
    
    /// Basis set on center one
    shared_ptr<BasisSet> basis();
    /// Basis set on center one
    shared_ptr<BasisSet> basis1();
    /// Basis set on center two
    shared_ptr<BasisSet> basis2();
    /// Basis set on center three
    shared_ptr<BasisSet> basis3();
    /// Basis set on center four
    shared_ptr<BasisSet> basis4();

    /// Buffer where the integrals are placed
    const double *buffer() const { return target_; };
    
    /// Compute ERIs between 4 shells. Result is stored in buffer.
    virtual void compute_shell(const ShellCombinationsIterator&) = 0;

    /// Compute the integrals
    virtual void compute_shell(int, int, int, int) = 0;

    /// Compute the integrals
    virtual void compute_shell_deriv1(int, int, int, int) = 0;
    
    /// Normalize Cartesian functions based on angular momentum
    void normalize_am(shared_ptr<GaussianShell>, shared_ptr<GaussianShell>, shared_ptr<GaussianShell>, shared_ptr<GaussianShell>, int nchunk=1);
        
    /// Return true if the clone member can be called. By default returns false.
    virtual bool cloneable();
    
    /// Returns a clone of this object. By default throws an exception
    virtual TwoBodyInt* clone();
    
    /// Results go back to buffer_
    void pure_transform(int, int, int, int, int nchunk);
};

}

#endif
