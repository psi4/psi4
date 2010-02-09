#ifndef _psi_src_lib_libmints_gshell_h_
#define _psi_src_lib_libmints_gshell_h_

#include <cstdio>
#include <libmints/vector3.h>

namespace psi {

/*! \ingroup MINTS
 *  \class GaussianShell
 *  \brief Gaussian orbital shell.
 */
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

    /** Initializes some basic data about the GaussianShell
     *
     *  Determines maximum and minimum AM and totall number of Cartesian and functions.
     */
    void init_data();
    /** Makes a copy of the data for the Gaussian shell.
     *  @param l Array of angular momentum
     *  @param exp Array of exponents
     *  @param coef Matrix of contraction coefficients. nprimitive x ncontraction
     */
    void copy_data(int *l, double *exp, double **coef);

    /** Normalizes a single primitive.
     *  @param p The primitive index to normalize.
     *  @return Normalization constant to be applied to the primitive.
     */
    double primitive_normalization(int p);
    /** Normalizes an entire contraction set. Applies the normalization to the coefficients
     *  @param gs The contraction set to normalize.
     */
    void contraction_normalization(int gs);
    /** Handles calling primitive_normalization and contraction_normalization for you. */
    void normalize_shell();

    /** Lookup array that when you index the angular momentum it returns the lowercase letter corresponding to it. */
    static const char *amtypes;
    /** Lookup array that when you index the angular momentum it returns the uppercase letter corresponding to it. */
    static const char *AMTYPES;

public:
    /** Constructor, does nothing. */
    GaussianShell() {};
    ~GaussianShell();

    /** Initializes the GaussianShell with the data provided.
     *  @param ncn The number of contractions for the shell (most likely you'll say 1).
     *  @param nprm The number of primitives in each contraction.
     *  @param e An array of exponent values.
     *  @param am An array of angular momentum of each contraction (most likely of length 1).
     *  @param pure Pure spherical harmonics, or Cartesian.
     *  @param c A matrix of contraction coefficients (most likely of dimensions nprm x 1).
     *  @param nc The atomic center that this shell is located on. Must map back to the correct atom in the owning BasisSet molecule_. Used in integral derivatives for indexing.
     *  @param center The x, y, z position of the shell. This is passed to reduce the number of calls to the molecule.
     *  @param start The starting index of the first function this shell provides. Used to provide starting positions in matrices.
     *  @param pt Is the shell already normalized?
     */
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
