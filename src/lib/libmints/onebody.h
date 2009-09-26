#ifndef _psi_src_lib_libmints_onebody_h_
#define _psi_src_lib_libmints_onebody_h_

/*!
    \file libmints/onebody.h
    \ingroup MINTS
*/

#include <libutil/ref.h>
#include <libmints/matrix.h>

namespace psi {
    
class IntegralFactory;
class BasisSet;
class GaussianShell;

/// Basis class for all one-electron integrals.
class OneBodyInt
{
protected:
    IntegralFactory *integral_;
    BasisSet* bs1_;
    BasisSet* bs2_;
    
    double *buffer_;
    unsigned int count_;
    int deriv_;
    int natom_;
    
    OneBodyInt(IntegralFactory *integral, BasisSet* bs1, BasisSet* bs2=0, int deriv=0);
    
public:
    virtual ~OneBodyInt();
    
    /// Basis set on center one.
    BasisSet* basis();
    /// Basis set on center one.
    BasisSet* basis1();
    /// Basis set on center two.
    BasisSet* basis2();
    
    /// Buffer where the integrals are placed.
    const double *buffer() const;
    
    /// Compute the integrals between basis function in the given shell pair.
    virtual void compute_shell(int, int) = 0;
    
    /// Computes all integrals and stores them in result
    void compute(Matrix* result);
    void compute(Matrix& result);
    
    /// Computes all integrals and stores them in result by default this method throws
    virtual void compute(Matrix** result);
    /// Computes all integrals and stores them in result by default this method throws
    virtual void compute(SimpleMatrix** result);
    
    /// Does the method provide first derivatives?
    virtual bool has_deriv1() { return false; }
    
    /// Does the method provide second derivatives?
    virtual bool has_deriv2() { return false; }

    /// Computes the first derivatives and stores them in result
    virtual void compute_deriv1(Matrix** result);
    /// Computes the first derivatives and stores them in result
    virtual void compute_deriv1(SimpleMatrix** result);
    /// Computes the second derivatives and stores them in result
    virtual void compute_deriv2(SimpleMatrix** result);
    
    /// Computes the integrals between basis function in the given shell pair
    virtual void compute_shell_deriv1(int, int);
    /// Computes the integrals between basis function in the given shell pair
    virtual void compute_shell_deriv2(int, int);
    
    /// Integral object that created me.
    IntegralFactory *integral() const { return integral_; }
    
    /// Return true if the clone member can be called. By default returns false.
    virtual bool cloneable();
    
    /// Returns a clone of this object. By default throws an exception.
    virtual OneBodyInt* clone();
    
    /// Normalize Cartesian functions based on angular momentum
    void normalize_am(GaussianShell*, GaussianShell*, int nchunk=1);
    
    /// Transform Cartesian integrals to spherical harmonic ones.
    /// Reads from buffer_ and stores results back in buffer_.
    virtual void spherical_transform(GaussianShell*, GaussianShell*);
    
    void do_transform(GaussianShell*, GaussianShell*, int);
    
    /// Accumulates results into a Matrix
    void so_transform(Matrix       *result, int, int, int ichunk=0);
    /// Accumulates results into a SimpleMatrix
    void so_transform(SimpleMatrix *result, int, int, int ichunk=0);
};

}

#endif
