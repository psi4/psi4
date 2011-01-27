#ifndef _psi_src_lib_libmints_gshell_h_
#define _psi_src_lib_libmints_gshell_h_

#include <cstdio>
#include <libmints/vector3.h>

namespace psi {

class Vector3;

enum PrimitiveType { Normalized, Unnormalized };
enum GaussianType { Cartesian, Pure };

/*! \ingroup MINTS
 *  \class GaussianShell
 *  \brief Gaussian orbital shell.
 */
class GaussianShell
{
private:
    /// number of primitives used in this contraction
    int nprimitives_;
    /// Angular momentum
    int l_;
    /// Flag for pure angular momentum
    int puream_;
    /// Exponents (of length nprimitives_)
    double *exp_;
    /// Contraction coefficients (of length nprimitives_)
    double *coef_;

    /// Atom number this shell goes to. Needed when indexing integral derivatives.
    int nc_;
    /// Atomic center number in the Molecule
    Vector3 center_;
    int start_;

    /// How many cartesian functions? (1=s, 3=p, 6=d, ...)
    int ncartesians_;
    /// How many functions? (1=s, 3=p, 6=d, ...)
    int nfunctions_;

    /** Initializes some basic data about the GaussianShell
     *
     *  Determines maximum and minimum AM and totall number of Cartesian and functions.
     */
    void init_data();

    /** Makes a copy of the data for the Gaussian shell.
     *  @param l Angular momentum
     *  @param exp Array of exponents
     *  @param coef Array of contraction coefficients. nprimitive
     */
    void copy_data(int l, double *exp, double *coef);

    /** Normalizes a single primitive.
     *  @param p The primitive index to normalize.
     *  @return Normalization constant to be applied to the primitive.
     */
    double primitive_normalization(int p);
    /** Normalizes an entire contraction set. Applies the normalization to the coefficients
     *  @param gs The contraction set to normalize.
     */
    void contraction_normalization();

    /** Lookup array that when you index the angular momentum it returns the lowercase letter corresponding to it. */
    static const char *amtypes;
    /** Lookup array that when you index the angular momentum it returns the uppercase letter corresponding to it. */
    static const char *AMTYPES;

public:
    /** Constructor, does nothing. */
    GaussianShell() {};
    ~GaussianShell();

    /** Handles calling primitive_normalization and contraction_normalization for you. */
    void normalize_shell();
    /** Initializes the GaussianShell with the data provided.
     *  @param nprm The number of primitives.
     *  @param e An array of exponent values.
     *  @param am Angular momentum.
     *  @param pure Pure spherical harmonics, or Cartesian.
     *  @param c An array of contraction coefficients.
     *  @param nc The atomic center that this shell is located on. Must map back to the correct atom in the owning BasisSet molecule_. Used in integral derivatives for indexing.
     *  @param center The x, y, z position of the shell. This is passed to reduce the number of calls to the molecule.
     *  @param start The starting index of the first function this shell provides. Used to provide starting positions in matrices.
     *  @param pt Is the shell already normalized?
     */
    void init(int nprm,
              double* e,
              int am,
              GaussianType pure,
              double* c,
              int nc,
              Vector3& center,
              int start,
              PrimitiveType pt = Normalized);

    /// The number of primitive Gaussians
    int nprimitive() const          { return nprimitives_; }
    /// Total number of basis functions
    int nfunction() const;
    /// Total number of functions if this shell was Cartesian
    int ncartesian() const          { return ncartesians_; }
    /// The angular momentum of the given contraction
    int am() const           { return l_; }
    /// The character symbol for the angular momentum of the given contraction
    char amchar() const      { return amtypes[l_]; }
    /// The character symbol for the angular momentum of the given contraction (upper case)
    char AMCHAR() const      { return AMTYPES[l_]; }
    /// Returns true if contraction is Cartesian
    bool is_cartesian() const  { return !puream_; }
    /// Returns true if contraction is pure
    bool is_pure() const       { return l_ < 2 ? false : puream_; }
//    bool is_pure() const       { puream_; }

    /// Returns the center of the Molecule this shell is on
    Vector3 center() const;
    /// Returns the atom number this shell is on. Used by integral derivatives for indexing.
    int ncenter() const             { return nc_; }

    /// Returns the exponent of the given primitive
    double exp(int prim) const { return exp_[prim]; }
    /// Return coefficient of pi'th primitive and ci'th contraction
    double coef(int pi) const { return coef_[pi]; }
    /// Returns the exponent of the given primitive
    double* exps() const { return exp_; }
    /// Return coefficient of pi'th primitive and ci'th contraction
    double* coefs() const { return coef_; }

    /// Print out the shell
    void print(FILE *out) const;

    /// Normalize the angular momentum component
    double normalize(int l, int m, int n);

    /// Basis function index where this shell starts.
    int function_index() const { return start_; }
    void set_function_index(int i) { start_ = i; }
    int ncartesians() const { return ncartesians_; }
    int nfunctions() const { return nfunctions_; }
};

}

#endif
