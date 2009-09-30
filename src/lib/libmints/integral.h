#ifndef _psi_src_lib_libmints_integral_h_
#define _psi_src_lib_libmints_integral_h_

/*!
    \file libmints/integral.h
    \ingroup MINTS
*/

#include <libutil/ref.h>
#include <libmints/basisset.h>
#include <vector>

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

namespace psi {
    
class BasisSet;
class OneBodyInt;
class TwoBodyInt;
class Symmetry;

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

class SphericalTransform
{
protected:
    std::vector<SphericalTransformComponent> components_;
    int l_; // The angular momentum this transform is for.
    
    SphericalTransform();
public:
    SphericalTransform(int l);
    virtual ~SphericalTransform() {};
        
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

class SphericalTransformIter
{
private:
    SphericalTransform& trans_;
    int i_;
    
public:
    SphericalTransformIter(SphericalTransform& trans) : trans_(trans) { i_ = 0; }
    
    void first() { i_ = 0; }
    void next()  { i_++;   }
    bool is_done() { return i_ < trans_.n() ? true : false; }
    
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
};

class IntegralsIterator
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
    shared_ptr<GaussianShell> usi, usj, usk, usl;

    bool done;

    //std::vector<Integral> unique_integrals_;
    int ii, iimax, jj, jjmax, kk, kkmax, ll, llmax;
    int ni, nj, nk, nl, fii, fij, fik, fil;
        
public:
    IntegralsIterator(shared_ptr<GaussianShell> s1, shared_ptr<GaussianShell> s2,
                      shared_ptr<GaussianShell> s3, shared_ptr<GaussianShell> s4) {
        done = false;
        usi = s1;
        usj = s2;
        usk = s3;
        usl = s4;
        ni =usi->nfunction(0);
        nj =usj->nfunction(0);
        nk =usk->nfunction(0);
        nl =usl->nfunction(0);

        fii = usi->function_index();
        fij = usj->function_index();
        fik = usk->function_index();
        fil = usl->function_index();
        
        iimax = ni - 1;
        if (usi == usj && usk == usl && usi == usk) {
            kkmax = 0;
            llmax = 0;
            jjmax = 0;
        }
        else if(usi == usk && usj == usl){
            kkmax = 0;
            llmax = 0;
            jjmax = nj - 1;
        }
        else{
            kkmax = nk - 1;
            jjmax = (usi == usj) ? 0 : nj - 1;
            llmax = (usk == usl) ? 0 : nl - 1;
        }

        ii = 0;
        jj = 0;
        kk = 0;
        ll = 0;
        
    }
    
    void first();
    void next();
    bool is_done() { return done; }
    
    int i() const { return current.i; }
    int j() const { return current.j; }
    int k() const { return current.k; }
    int l() const { return current.l; }
    int index() const { return current.index;}
};

class ShellCombinationsIterator
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
    
    shared_ptr<BasisSet> bs1_;
    shared_ptr<BasisSet> bs2_;
    shared_ptr<BasisSet> bs3_;
    shared_ptr<BasisSet> bs4_;
    
//    void generate_combinations(BasisSet*bs1, BasisSet*bs2,
//        BasisSet*bs3, BasisSet*bs4);
        
public:
    ShellCombinationsIterator(shared_ptr<BasisSet>bs1, shared_ptr<BasisSet>bs2,
                              shared_ptr<BasisSet>bs3, shared_ptr<BasisSet>bs4) :
                       bs1_(bs1), bs2_(bs2), bs3_(bs3), bs4_(bs4) {    }
                 
    void first();
    void next();
    bool is_done() { return done; }
    
    int p() const { return current.P; }
    int q() const { return current.Q; }
    int r() const { return current.R; }
    int s() const { return current.S; }
    int end_of_PK() const { return current.end_of_PK; }
    
    IntegralsIterator integrals_iterator();
};

class IntegralFactory
{
protected:
    /// Center 1 basis set
    shared_ptr<BasisSet> bs1_;
    /// Center 2 basis set
    shared_ptr<BasisSet> bs2_;
    /// Center 3 basis set
    shared_ptr<BasisSet> bs3_;
    /// Center 4 basis set
    shared_ptr<BasisSet> bs4_;
    
    /// Provides ability to transform to and from sphericals (d=0, f=1, g=2)
    std::vector<SphericalTransform> spherical_transforms_;
    
public:
    /** Initialize IntegralFactory object given a GaussianBasisSet for each center. */
    IntegralFactory(shared_ptr<BasisSet> bs1, shared_ptr<BasisSet> bs2,
                    shared_ptr<BasisSet> bs3, shared_ptr<BasisSet> bs4);
    
    virtual ~IntegralFactory();
    
    /// Set the basis set for each center.
    virtual void set_basis(shared_ptr<BasisSet> bs1, shared_ptr<BasisSet> bs2,
        shared_ptr<BasisSet> bs3, shared_ptr<BasisSet> bs4);
        
    /// Returns an OneBodyInt that computes the overlap integral.
    virtual OneBodyInt* overlap(int deriv=0);
    
    /// Returns an OneBodyInt that computes the kinetic energy integral.
    virtual OneBodyInt* kinetic(int deriv=0);
    
    /// Returns an OneBodyInt that computes the nuclear attraction integral.
    virtual OneBodyInt* potential(int deriv=0);

    /// Returns an OneBodyInt that computes the dipole integral.
    virtual OneBodyInt* dipole(int deriv=0);
    
    /// Returns an OneBodyInt that computes the quadrupole integral.
    virtual OneBodyInt* quadrupole();
    
    /// Returns an ERI integral object
    virtual TwoBodyInt* eri(int deriv=0);

    /// Returns an ERI iterator object, only coded for standard ERIs
    ShellCombinationsIterator shells_iterator();
    
    /// Initializes spherical harmonic transformations
    virtual void init_spherical_harmonics(int max_am);
    
    // Return spherical transform object for am
    SphericalTransform* spherical_transform(int am) { return &(spherical_transforms_[am]); }
};

}

#endif
