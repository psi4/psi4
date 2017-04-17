/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#ifndef _psi_src_lib_libmints_integral_h_
#define _psi_src_lib_libmints_integral_h_

 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
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
    const GaussianShell& usi;
    const GaussianShell& usj;
    const GaussianShell& usk;
    const GaussianShell& usl;

    bool done;

    int ii, iimax, jj, jjmax, kk, kkmax, ll, llmax;
    int ni, nj, nk, nl, fii, fij, fik, fil;

public:
    AOIntegralsIterator(const GaussianShell& s1, const GaussianShell& s2,
                      const GaussianShell& s3, const GaussianShell& s4);

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

    std::shared_ptr<BasisSet> bs1_;
    std::shared_ptr<BasisSet> bs2_;
    std::shared_ptr<BasisSet> bs3_;
    std::shared_ptr<BasisSet> bs4_;

public:
    AOShellCombinationsIterator(std::shared_ptr<BasisSet>bs1, std::shared_ptr<BasisSet>bs2,
                              std::shared_ptr<BasisSet>bs3, std::shared_ptr<BasisSet>bs4);
    AOShellCombinationsIterator();
    void init(std::shared_ptr<BasisSet>bs1, std::shared_ptr<BasisSet>bs2,
            std::shared_ptr<BasisSet>bs3, std::shared_ptr<BasisSet>bs4);

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

    std::shared_ptr<SOBasisSet> bs1_;
    std::shared_ptr<SOBasisSet> bs2_;
    std::shared_ptr<SOBasisSet> bs3_;
    std::shared_ptr<SOBasisSet> bs4_;

public:
    SOShellCombinationsIterator(std::shared_ptr<SOBasisSet>bs1, std::shared_ptr<SOBasisSet>bs2,
                              std::shared_ptr<SOBasisSet>bs3, std::shared_ptr<SOBasisSet>bs4);
    SOShellCombinationsIterator();
    void init(std::shared_ptr<SOBasisSet>bs1, std::shared_ptr<SOBasisSet>bs2,
            std::shared_ptr<SOBasisSet>bs3, std::shared_ptr<SOBasisSet>bs4);

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

    std::shared_ptr<SOBasisSet> bs1_;

public:
    SO_PQ_Iterator(std::shared_ptr<SOBasisSet>bs1);
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

    std::shared_ptr<SOBasisSet> bs1_;
    std::shared_ptr<SOBasisSet> bs2_;
    std::shared_ptr<SOBasisSet> bs3_;
    std::shared_ptr<SOBasisSet> bs4_;

public:
    SO_RS_Iterator(const int &P, const int &Q,
                   std::shared_ptr<SOBasisSet>bs1, std::shared_ptr<SOBasisSet>bs2,
                   std::shared_ptr<SOBasisSet>bs3, std::shared_ptr<SOBasisSet>bs4);
    SO_RS_Iterator(std::shared_ptr<SOBasisSet>bs1, std::shared_ptr<SOBasisSet>bs2,
                   std::shared_ptr<SOBasisSet>bs3, std::shared_ptr<SOBasisSet>bs4);

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
    std::shared_ptr<BasisSet> bs1_;
    /// Center 2 basis set
    std::shared_ptr<BasisSet> bs2_;
    /// Center 3 basis set
    std::shared_ptr<BasisSet> bs3_;
    /// Center 4 basis set
    std::shared_ptr<BasisSet> bs4_;
    /// ECP basis set
    std::shared_ptr<BasisSet> bs_ecp_;

    /// Provides ability to transform to sphericals (d=0, f=1, g=2)
    std::vector<SphericalTransform> spherical_transforms_;
    /// Provides ability to transform from sphericals (d=0, f=1, g=2)
    std::vector<ISphericalTransform> ispherical_transforms_;

public:
    /** Initialize IntegralFactory object given a BasisSet for each center, and optionally an ECPBasisSet. */
    IntegralFactory(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                    std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4);
    IntegralFactory(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                    std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4,
                    std::shared_ptr<BasisSet> bsecp);
    /** Initialize IntegralFactory object given a BasisSet for two centers, and optionally an ECPBasisSet. Becomes (bs1 bs1 | bs1 bs1). */
    IntegralFactory(std::shared_ptr<BasisSet> bs1);
    IntegralFactory(std::shared_ptr<BasisSet> bs1,  std::shared_ptr<BasisSet> bsecp);

    virtual ~IntegralFactory();

    /// Return the basis set on center 1.
    std::shared_ptr<BasisSet> basis1() const;
    /// Return the basis set on center 2.
    std::shared_ptr<BasisSet> basis2() const;
    /// Return the basis set on center 3.
    std::shared_ptr<BasisSet> basis3() const;
    /// Return the basis set on center 4.
    std::shared_ptr<BasisSet> basis4() const;
    /// Return the ECP basis set
    std::shared_ptr<BasisSet> basisECP() const;
    /// Return whether ECPs are present
    bool hasECP() const; 

    /// Set the basis set for each center.
    virtual void set_basis(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
        std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4,  std::shared_ptr<BasisSet> bsecp);

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
    
    /// Returns an OneBodyInt that computes the ECP integral.
    virtual OneBodyAOInt* ao_ecp(int deriv=0);
    virtual OneBodySOInt* so_ecp(int deriv=0);

    /// Returns an OneBodyInt that computes the relativistic nuclear attraction integral.
    virtual OneBodyAOInt* ao_rel_potential(int deriv=0);
    virtual OneBodySOInt* so_rel_potential(int deriv=0);

    /// Returns the OneBodyInt that computes the pseudospectral grid integrals
    virtual OneBodyAOInt* ao_pseudospectral(int deriv = 0);
    virtual OneBodySOInt* so_pseudospectral(int deriv = 0);

    /// Returns an OneBodyInt that computes the dipole integral.
    virtual OneBodyAOInt* ao_dipole(int deriv=0);
    virtual OneBodySOInt* so_dipole(int deriv=0);

    /// Returns an OneBodyInt that computes the quadrupole integral.
    virtual OneBodyAOInt* ao_quadrupole();
    virtual OneBodySOInt* so_quadrupole();

    /// Returns an OneBodyInt that computes arbitrary-order multipole integrals.
    virtual OneBodyAOInt* ao_multipoles(int order);
    virtual OneBodySOInt* so_multipoles(int order);

    /// Returns an OneBodyInt that computes the traceless quadrupole integral.
    virtual OneBodyAOInt* ao_traceless_quadrupole();
    virtual OneBodySOInt* so_traceless_quadrupole();

    /// Returns an OneBodyInt that computes the nabla integral.
    virtual OneBodyAOInt* ao_nabla(int deriv=0);
    virtual OneBodySOInt* so_nabla(int deriv=0);

    /// Returns an OneBodyInt that computes the nabla integral.
    virtual OneBodyAOInt* ao_angular_momentum(int deriv=0);
    virtual OneBodySOInt* so_angular_momentum(int deriv=0);

    /// Returns a OneBodyInt that computes the multipole potential integrals for EFP
    virtual OneBodyAOInt* ao_efp_multipole_potential(int deriv=0);
    virtual OneBodySOInt* so_efp_multipole_potential(int deriv=0);

    /// Returns an OneBodyInt that computes the electric field
    virtual OneBodyAOInt *electric_field();

    /// Returns an OneBodyInt that computes the point electrostatic potential
    virtual OneBodyAOInt *electrostatic();

    /// Returns an OneBodyInt that computes the electrostatic potential at desired points
    /// Want to change the name of this after the PCM dust settles
    virtual OneBodyAOInt *pcm_potentialint();

    /// Returns an ERI integral object
    virtual TwoBodyAOInt* eri(int deriv=0, bool use_shell_pairs=true);

    /// Returns an ERD ERI integral object, if available.  Otherwise returns a libint integral object
    virtual TwoBodyAOInt* erd_eri(int deriv=0, bool use_shell_pairs=true);

    /// Returns an erf ERI integral object (omega integral)
    virtual TwoBodyAOInt* erf_eri(double omega, int deriv=0, bool use_shell_pairs=true);

    /// Returns an erf complement ERI integral object (omega integral)
    virtual TwoBodyAOInt* erf_complement_eri(double omega, int deriv=0, bool use_shell_pairs=true);

    /// Returns an F12 integral object
    virtual TwoBodyAOInt* f12(std::shared_ptr<CorrelationFactor> cf, int deriv=0, bool use_shell_pairs=true);

    /// Returns an F12Scaled integral object
    virtual TwoBodyAOInt* f12_scaled(std::shared_ptr<CorrelationFactor> cf, int deriv=0, bool use_shell_pairs=true);

    /// Returns an F12 squared integral object
    virtual TwoBodyAOInt* f12_squared(std::shared_ptr<CorrelationFactor> cf, int deriv=0, bool use_shell_pairs=true);

    /// Returns an F12G12 integral object
    virtual TwoBodyAOInt* f12g12(std::shared_ptr<CorrelationFactor> cf, int deriv=0, bool use_shell_pairs=true);

    /// Returns an F12 double commutator integral object
    virtual TwoBodyAOInt* f12_double_commutator(std::shared_ptr<CorrelationFactor> cf, int deriv=0, bool use_shell_pairs=true);

    /// Returns a general ERI iterator object for any (P Q | R S) in shells
    AOIntegralsIterator integrals_iterator(int p, int q, int r, int s);

    /// Returns an ERI iterator object, only coded for standard ERIs
    AOShellCombinationsIterator shells_iterator();
    AOShellCombinationsIterator* shells_iterator_ptr();

    /// Initializes spherical harmonic transformations
    virtual void init_spherical_harmonics(int max_am);

    /// Return spherical transform object for am
    const SphericalTransform* spherical_transform(int am) const { return &(spherical_transforms_[am]); }

    // Return spherical transform object for am
    std::vector<SphericalTransform>& spherical_transform() { return spherical_transforms_; }

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
