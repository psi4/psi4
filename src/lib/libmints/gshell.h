#ifndef _psi_src_lib_libmints_gshell_h_
#define _psi_src_lib_libmints_gshell_h_

/*!
    \file libmints/gshell.h
    \ingroup MINTS
*/

#include <cstdio>
#include <libmints/vector3.h>

namespace psi {
    
/// A Gaussian orbital shell.
class GaussianShell
{
public:
    enum PrimitiveType { Normalized, Unnormalized };
    enum GaussianType { Cartesian, Pure };

private:
    /// number of primitives used in this contraction
    int nprimitives_;
    /// number of contractions
    int ncontractions_;
    /// Angular momentum for each contraction (length ncontractions_)
    int *l_;
    /// Flag for pure angular momentum for each contraction (length ncontractions_)
    int *puream_;
    /// Exponents (of length nprimitives_)
    double *exp_;
    /// Contraction coefficients (ncontractions_ x nprimitives_)
    double **coef_;

    /// Atom number this shell goes to. Needed when indexing integral derivatives.
    int nc_;
    /// Atomic center number in the Molecule
    Vector3 center_;
    int start_;
    
    /// How many cartesian functions? (1=s, 3=p, 6=d, ...)
    int ncartesians_;
    /// How many functions? (1=s, 3=p, 6=d, ...)
    int nfunctions_;
    /// Worry about pure angular momentum?
    bool has_pure_;

    int max_am_;
    int min_am_;
    
    void init_data();
    void copy_data(int *l, double *exp, double **coef);

    double primitive_normalization(int);
    void contraction_normalization(int);
    void normalize_shell();
    
    static const char *amtypes;
    static const char *AMTYPES;

public:
    GaussianShell() {};
    ~GaussianShell();
    
    void init(int ncn, int nprm, double* e, int* am, GaussianType pure,
        double** c, int nc, Vector3& center, int start, PrimitiveType pt = GaussianShell::Normalized);

    /// The number of primitive Gaussians
    int nprimitive() const          { return nprimitives_; }
    /// The number of contractions formed from the primitives
    int ncontraction() const        { return ncontractions_; }
    /// The number of basis functions
    int nfunction(int) const;
    /// Total number of basis functions
    int nfunction() const           { return nfunctions_; }
    /// Total number of functions if this shell was Cartesian
    int ncartesian() const          { return ncartesians_; }
    /// The number of Cartesian functions for the given contraction
    int ncartesian(int c) const     { return ((l_[c]+2)*(l_[c]+1))>>1; }
    /// The angular momentum of the given contraction
    int am(int con) const           { return l_[con]; }
    /// The minimum angular momentum
    int min_am() const              { return min_am_; }
    /// The maximum angular momentum
    int max_am() const              { return max_am_; }
    /// The character symbol for the angular momentum of the given contraction
    char amchar(int con) const      { return amtypes[l_[con]]; }
    /// Returns true if contraction is Cartesian
    bool is_cartesian(int c) const  { return !puream_[c]; }
    /// Returns true if contraction is pure
    bool is_pure(int c) const       { return puream_[c]; }
    /// Return true if any contraction is pure
    bool has_pure() const           { return has_pure_; }
    
    /// Returns the center of the Molecule this shell is on
    Vector3 center() const          { return center_; }
    /// Returns the atom number this shell is on. Used by integral derivatives for indexing.
    int ncenter() const             { return nc_; }
    
    /// Returns the exponent of the given primitive
    double exp(int prim) const { return exp_[prim]; }
    /// Return coefficient of pi'th primitive and ci'th contraction
    double coef(int ci, int pi) const { return coef_[ci][pi]; }
    
    /// Print out the shell
    void print(FILE *out) const;
    
    /// Normalize the angular momentum component
    double normalize(int l, int m, int n);
    
    /// Basis function index where this shell starts.
    int function_index() const { return start_; }
};

}

#endif
