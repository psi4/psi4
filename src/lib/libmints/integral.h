#ifndef _psi_src_lib_libmints_integral_h_
#define _psi_src_lib_libmints_integral_h_

#include <boost/shared_ptr.hpp>
#include <vector>

#include "onebody.h"
#include "twobody.h"

/*! \def INT_NCART(am)
    Gives the number of cartesian functions for an angular momentum.
*/
#define INT_NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)
/*! \def INT_PURE(am)
    Gives the number of spherical functions for an angular momentum.
*/
#define INT_NPURE(am) (2*(am)+1)
/*! \def INT_NFUNC(pu,am)
    Gives the number of functions for an angular momentum based on pu.
*/
#define INT_NFUNC(pu,am) ((pu)?INT_NPURE(am):INT_NCART(am))
/*! \def INT_CARTINDEX(am,i,j)
    Computes offset index for cartesian function.
*/
#define INT_CARTINDEX(am,i,j) (((i) == (am))? 0 : (((((am) - (i) + 1)*((am) - (i)))>>1) + (am) - (i) - (j)))
/*! \def INT_ICART(a, b, c)
    Given a, b, and c compute a cartesian offset.
*/
#define INT_ICART(a, b, c) (((((((a)+(b)+(c)+1)<<1)-(a))*((a)+1))>>1)-(b)-1)
/*! \def INT_IPURE(l, m)
    Given l and m compute a pure function offset.
*/
#define INT_IPURE(l, m) ((l)+(m))

namespace psi {

class BasisSet;
class GaussianShell;
class OneBodyAOInt;
class OneBodySOInt;
class TwoBodyAOInt;
class ThreeCenterOverlapInt;
class Symmetry;
class CartesianIter;
class RedundantCartesianIter;
class RedundantCartesianSubIter;
class ShellRotation;
class SymmetryOperation;
class SOTransform;
class SOBasisSet;
class CorrelationFactor;

/*! \ingroup MINTS */
class SphericalTransformComponent
{
protected:
    int a_, b_, c_;
    int cartindex_, pureindex_;

    double coef_;

public:
    /// Returns the exponent of x.
    int a() const { return a_; }
    /// Returns the exponent of y.
    int b() const { return b_; }
    /// Returns the exponent of z.
    int c() const { return c_; }
    /// Returns the index of the Cartesian basis function
    int cartindex() const { return cartindex_; }
    /// Returns the index of the spherical harmonic basis function
    int pureindex() const { return pureindex_; }
    /// Returns the coefficient of this component of the transformation
    double coef() const { return coef_; }

    void init(int a, int b, int c, double coef, int cartindex, int pureindex);
};

/*! \ingroup MINTS */
/** This is a base class for a container for a sparse Cartesian to solid
    harmonic basis function transformation. */
class SphericalTransform
{
protected:
    std::vector<SphericalTransformComponent> components_;
    int l_; // The angular momentum this transform is for.
    int subl_;

    SphericalTransform();

    virtual void init();
public:
    SphericalTransform(int l, int subl = -1);
    virtual ~SphericalTransform() {}

    /// Returns the Cartesian basis function index of component i
    int cartindex(int i) const { return components_[i].cartindex(); }
    /// Returns the spherical harmonic basis index of component i
    int pureindex(int i) const { return components_[i].pureindex(); }
    /// Returns the transformation coefficient of component i
    double coef(int i) const { return components_[i].coef(); }
    /// Returns the Cartesian basis function's x exponent of component i
    int a(int i) const { return components_[i].a(); }
    /// Returns the Cartesian basis function's y exponent of component i
    int b(int i) const { return components_[i].b(); }
    /// Returns the Cartesian basis function's z exponent of component i
    int c(int i) const { return components_[i].c(); }
    /// Returns the number of components in the transformation
    int n() const { return components_.size(); }
    /// Returns the angular momentum
    int l() const { return l_; }
};

/*! \ingroup MINTS */
/// This describes a solid harmonic to Cartesian transformation.
class ISphericalTransform : public SphericalTransform
{
protected:
    ISphericalTransform();
    virtual void init();
public:
    ISphericalTransform(int l, int subl=-1);
};

class SphericalTransformIter
{
private:
    const SphericalTransform& trans_;
    int i_;

public:
    SphericalTransformIter(const SphericalTransform& trans) : trans_(trans) { i_ = 0; }

    void first() { i_ = 0; }
    void next()  { i_++;   }
    bool is_done() { return i_ < trans_.n() ? false : true; }

    /// Returns how many transforms are in this iterator.
    int n() const { return trans_.n(); }

    /// Returns the Cartesian basis function index of component i
    int cartindex() const { return trans_.cartindex(i_); }
    /// Returns the spherical harmonic basis index of component i
    int pureindex() const { return trans_.pureindex(i_); }
    /// Returns the transformation coefficient of component i
    double coef()   const { return trans_.coef(i_); }
    /// Returns the Cartesian basis function's x exponent of component i
    int a()         const { return trans_.a(i_); }
    /// Returns the Cartesian basis function's y exponent of component i
    int b()         const { return trans_.b(i_); }
    /// Returns the Cartesian basis function's z exponent of component i
    int c()         const { return trans_.c(i_); }
    /// Return a component of the transform.
    int l(int i) { return i?(i==1?b():c()):a(); }
};

/*! \ingroup MINTS */
class AOIntegralsIterator
{
private:
    struct Integral {
        int i;
        int j;
        int k;
        int l;
        unsigned int index;
    };

    Integral current;
    boost::shared_ptr<GaussianShell> usi, usj, usk, usl;

    bool done;

    int ii, iimax, jj, jjmax, kk, kkmax, ll, llmax;
    int ni, nj, nk, nl, fii, fij, fik, fil;

public:
    AOIntegralsIterator(boost::shared_ptr<GaussianShell> s1, boost::shared_ptr<GaussianShell> s2,
                      boost::shared_ptr<GaussianShell> s3, boost::shared_ptr<GaussianShell> s4);

    void first();
    void next();
    bool is_done() { return done; }

    int i() const { return current.i; }
    int j() const { return current.j; }
    int k() const { return current.k; }
    int l() const { return current.l; }
    int index() const { return current.index;}
};


/*! \ingroup MINTS */
class AOShellCombinationsIterator
{
private:
    struct ShellQuartet {
        int P;
        int Q;
        int R;
        int S;
        bool end_of_PK;
    };

    ShellQuartet current;
    int usi_arr[3], usj_arr[3], usk_arr[3], usl_arr[3];
    int usii, usjj, uskk, usll, upk;

    int num_unique_pk;

    bool done;

    boost::shared_ptr<BasisSet> bs1_;
    boost::shared_ptr<BasisSet> bs2_;
    boost::shared_ptr<BasisSet> bs3_;
    boost::shared_ptr<BasisSet> bs4_;

public:
    AOShellCombinationsIterator(boost::shared_ptr<BasisSet>bs1, boost::shared_ptr<BasisSet>bs2,
                              boost::shared_ptr<BasisSet>bs3, boost::shared_ptr<BasisSet>bs4);
    AOShellCombinationsIterator();
    void init(boost::shared_ptr<BasisSet>bs1, boost::shared_ptr<BasisSet>bs2,
            boost::shared_ptr<BasisSet>bs3, boost::shared_ptr<BasisSet>bs4);

    void first();
    void next();
    bool is_done() { return done; }

    int p() const { return current.P; }
    int q() const { return current.Q; }
    int r() const { return current.R; }
    int s() const { return current.S; }
    int end_of_PK() const { return current.end_of_PK; }

    AOIntegralsIterator integrals_iterator();
};

class SOShellCombinationsIterator
{
private:
    struct ShellQuartet {
        int P;
        int Q;
        int R;
        int S;
        bool end_of_PK;
    };

    ShellQuartet current;
    int usi_arr[3], usj_arr[3], usk_arr[3], usl_arr[3];
    int usii, usjj, uskk, usll, upk;

    int num_unique_pk;

    bool done;

    boost::shared_ptr<SOBasisSet> bs1_;
    boost::shared_ptr<SOBasisSet> bs2_;
    boost::shared_ptr<SOBasisSet> bs3_;
    boost::shared_ptr<SOBasisSet> bs4_;

public:
    SOShellCombinationsIterator(boost::shared_ptr<SOBasisSet>bs1, boost::shared_ptr<SOBasisSet>bs2,
                              boost::shared_ptr<SOBasisSet>bs3, boost::shared_ptr<SOBasisSet>bs4);
    SOShellCombinationsIterator();
    void init(boost::shared_ptr<SOBasisSet>bs1, boost::shared_ptr<SOBasisSet>bs2,
            boost::shared_ptr<SOBasisSet>bs3, boost::shared_ptr<SOBasisSet>bs4);

    void first();
    void next();
    bool is_done() { return done; }

    int p() const { return current.P; }
    int q() const { return current.Q; }
    int r() const { return current.R; }
    int s() const { return current.S; }
    int end_of_PK() const { return current.end_of_PK; }
};

class SO_PQ_Iterator
{
private:
    struct PQ_Pair {
        int P;
        int Q;
    };

    PQ_Pair current;
    int ii, jj;

    bool done;

    boost::shared_ptr<SOBasisSet> bs1_;

public:
    SO_PQ_Iterator(boost::shared_ptr<SOBasisSet>bs1);
    SO_PQ_Iterator();

    void first();
    void next();
    bool is_done() { return done; }

    int p() const { return current.P; }
    int q() const { return current.Q; }
};

class SO_RS_Iterator
{
private:
    struct RS_Pair {
        int P;
        int Q;
        int R;
        int S;
    };

    RS_Pair current;
    int usi_arr[3], usj_arr[3], usk_arr[3], usl_arr[3];
    int usii, usjj, uskk, usll, upk;

    int num_unique_pk;

    bool done;

    boost::shared_ptr<SOBasisSet> bs1_;
    boost::shared_ptr<SOBasisSet> bs2_;
    boost::shared_ptr<SOBasisSet> bs3_;
    boost::shared_ptr<SOBasisSet> bs4_;

public:
    SO_RS_Iterator(const int &P, const int &Q,
                   boost::shared_ptr<SOBasisSet>bs1, boost::shared_ptr<SOBasisSet>bs2,
                   boost::shared_ptr<SOBasisSet>bs3, boost::shared_ptr<SOBasisSet>bs4);
    SO_RS_Iterator(boost::shared_ptr<SOBasisSet>bs1, boost::shared_ptr<SOBasisSet>bs2,
                   boost::shared_ptr<SOBasisSet>bs3, boost::shared_ptr<SOBasisSet>bs4);

    SO_RS_Iterator();

    void first();
    void next();
    bool is_done() { return done; }

    int p() const { return current.P; }
    int q() const { return current.Q; }
    int r() const { return current.R; }
    int s() const { return current.S; }
};


/*! \ingroup MINTS */
class IntegralFactory
{
protected:
    /// Center 1 basis set
    boost::shared_ptr<BasisSet> bs1_;
    /// Center 2 basis set
    boost::shared_ptr<BasisSet> bs2_;
    /// Center 3 basis set
    boost::shared_ptr<BasisSet> bs3_;
    /// Center 4 basis set
    boost::shared_ptr<BasisSet> bs4_;

    /// Provides ability to transform to sphericals (d=0, f=1, g=2)
    std::vector<SphericalTransform> spherical_transforms_;
    /// Provides ability to transform from sphericals (d=0, f=1, g=2)
    std::vector<ISphericalTransform> ispherical_transforms_;

public:
    /** Initialize IntegralFactory object given a BasisSet for each center. */
    IntegralFactory(boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2,
                    boost::shared_ptr<BasisSet> bs3, boost::shared_ptr<BasisSet> bs4);
    /** Initialize IntegralFactory object given a BasisSet for two centers. Becomes (bs1 bs2 | bs1 bs2). */
    IntegralFactory(boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2);
    /** Initialize IntegralFactory object given a BasisSet for two centers. Becomes (bs1 bs1 | bs1 bs1). */
    IntegralFactory(boost::shared_ptr<BasisSet> bs1);

    virtual ~IntegralFactory();

    /// Return the basis set on center 1.
    boost::shared_ptr<BasisSet> basis1() const;
    /// Return the basis set on center 2.
    boost::shared_ptr<BasisSet> basis2() const;
    /// Return the basis set on center 3.
    boost::shared_ptr<BasisSet> basis3() const;
    /// Return the basis set on center 4.
    boost::shared_ptr<BasisSet> basis4() const;

    /// Set the basis set for each center.
    virtual void set_basis(boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2,
        boost::shared_ptr<BasisSet> bs3, boost::shared_ptr<BasisSet> bs4);

    /// Returns an OneBodyInt that computes the overlap integral.
    virtual OneBodyAOInt* ao_overlap(int deriv=0);

    /// Returns an OneBodyInt that computes the overlap integral.
    virtual OneBodySOInt* so_overlap(int deriv=0);

    /// Returns a ThreeCenterOverlapINt that computes the overlap between three centers
    virtual ThreeCenterOverlapInt* overlap_3c();

    /// Returns an OneBodyInt that computes the kinetic energy integral.
    virtual OneBodyAOInt* ao_kinetic(int deriv=0);
    virtual OneBodySOInt* so_kinetic(int deriv=0);

    /// Returns an OneBodyInt that computes the nuclear attraction integral.
    virtual OneBodyAOInt* ao_potential(int deriv=0);
    virtual OneBodySOInt* so_potential(int deriv=0);

    /// Returns the OneBodyInt that computes the pseudospectral grid integrals
    virtual OneBodyAOInt* ao_pseudospectral(int deriv = 0);
    virtual OneBodySOInt* so_pseudospectral(int deriv = 0);

    /// Returns an OneBodyInt that computes the dipole integral.
    virtual OneBodyAOInt* ao_dipole(int deriv=0);
    virtual OneBodySOInt* so_dipole(int deriv=0);

    /// Returns an OneBodyInt that computes the quadrupole integral.
    virtual OneBodyAOInt* ao_quadrupole();
    virtual OneBodySOInt* so_quadrupole();

    /// Returns an OneBodyInt that computes the traceless quadrupole integral.
    virtual OneBodyAOInt* ao_traceless_quadrupole();
    virtual OneBodySOInt* so_traceless_quadrupole();

    /// Returns an OneBodyInt that computes the nabla integral.
    virtual OneBodyAOInt* ao_nabla(int deriv=0);
    virtual OneBodySOInt* so_nabla(int deriv=0);

    /// Returns an OneBodyInt that computes the nabla integral.
    virtual OneBodyAOInt* ao_angular_momentum(int deriv=0);
    virtual OneBodySOInt* so_angular_momentum(int deriv=0);

    /// Returns an OneBodyInt that computes the electric field
    virtual OneBodyAOInt *electric_field();

    /// Returns an OneBodyInt that computes the point electrostatic potential
    virtual OneBodyAOInt *electrostatic();

    /// Returns an ERI integral object
    virtual TwoBodyAOInt* eri(int deriv=0, double schwarz = 0.0);

    /// Returns an erf ERI integral object (omega integral)
    virtual TwoBodyAOInt* erf_eri(double omega, int deriv=0, double schwarz = 0.0);

    /// Returns an erf complement ERI integral object (omega integral)
    virtual TwoBodyAOInt* erf_complement_eri(double omega, int deriv=0, double schwarz = 0.0);

    /// Returns an F12 integral object
    virtual TwoBodyAOInt* f12(boost::shared_ptr<CorrelationFactor> cf, int deriv=0, double schwarz = 0.0);

    /// Returns an F12 squared integral object
    virtual TwoBodyAOInt* f12_squared(boost::shared_ptr<CorrelationFactor> cf, int deriv=0, double schwarz = 0.0);

    /// Returns an F12G12 integral object
    virtual TwoBodyAOInt* f12g12(boost::shared_ptr<CorrelationFactor> cf, int deriv=0, double schwarz = 0.0);

    /// Returns an F12 double commutator integral object
    virtual TwoBodyAOInt* f12_double_commutator(boost::shared_ptr<CorrelationFactor> cf, int deriv=0, double schwarz = 0.0);

    /// Returns a general ERI iterator object for any (P Q | R S) in shells
    AOIntegralsIterator integrals_iterator(int p, int q, int r, int s);

    /// Returns an ERI iterator object, only coded for standard ERIs
    AOShellCombinationsIterator shells_iterator();

    /// Initializes spherical harmonic transformations
    virtual void init_spherical_harmonics(int max_am);

    /// Return spherical transform object for am
    const SphericalTransform* spherical_transform(int am) const { return &(spherical_transforms_[am]); }

    // Return spherical transform object for am
    std::vector<SphericalTransform> spherical_transform() { return spherical_transforms_; }

    /// Return a spherical transform iterator object for am
    SphericalTransformIter* spherical_transform_iter(int am, int inv=0, int subl=-1) const;

    /// Return a new Cartesian iterator
    CartesianIter* cartesian_iter(int l) const;
    /// Return a new rudundant Cartesian iterator
    RedundantCartesianIter* redundant_cartesian_iter(int l) const;
    /// Return a new rudundant Cartesian sub iterator
    RedundantCartesianSubIter* redundant_cartesian_sub_iter(int l) const;
    /** Return the ShellRotation object for a shell of a given angular
        momentum. Pass nonzero to pure to do solid harmonics. */
    ShellRotation shell_rotation(int am, SymmetryOperation&, int pure=0) const;
};

}

#endif
