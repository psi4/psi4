#ifndef _psi_src_lib_libmints_gshell_h_
#define _psi_src_lib_libmints_gshell_h_

#include <cstdio>
#include <libmints/vector3.h>
#include <vector>

namespace psi {

class Vector3;

enum PrimitiveType {
    Normalized,
    Unnormalized
};

enum GaussianType {
    Cartesian = 0,
    Pure = 1
};

/*! \ingroup MINTS
 *  \class GaussianShell
 *  \brief Gaussian orbital shell.
 */
class GaussianShell
{
private:
    /// Angular momentum
    int l_;
    /// Flag for pure angular momentum
    int puream_;
    /// Exponents (of length nprimitives_)
    std::vector<double> exp_;
    /// Contraction coefficients (of length nprimitives_)
    std::vector<double> coef_;

    /// Atom number this shell goes to. Needed when indexing integral derivatives.
    int nc_;
    /// Atomic center number in the Molecule
    Vector3 center_;
    int start_;

    /// How many cartesian functions? (1=s, 3=p, 6=d, ...)
    int ncartesian_;
    /** How many functions? (1=s, 3=p, 5/6=d, ...)
     * Dependent on the value of puream_
     */
    int nfunction_;

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
    /** Constructor.
     *  @param e An array of exponent values.
     *  @param am Angular momentum.
     *  @param pure Pure spherical harmonics, or Cartesian.
     *  @param c An array of contraction coefficients.
     *  @param nc The atomic center that this shell is located on. Must map back to the correct atom in the owning BasisSet molecule_. Used in integral derivatives for indexing.
     *  @param center The x, y, z position of the shell. This is passed to reduce the number of calls to the molecule.
     *  @param start The starting index of the first function this shell provides. Used to provide starting positions in matrices.
     *  @param pt Is the shell already normalized?
     */
    GaussianShell(int am,
                  const std::vector<double>& c,
                  const std::vector<double>& e,
                  GaussianType pure,
                  int nc,
                  const Vector3& center,
                  int start,
                  PrimitiveType pt = Normalized);

    /** Handles calling primitive_normalization and contraction_normalization for you. */
    void normalize_shell();

    /// Make a copy of the GaussianShell.
    GaussianShell copy();

    /// Make a copy of the GaussianShell.
    GaussianShell copy(int nc, const Vector3& c);

    /// The number of primitive Gaussians
    int nprimitive() const;
    /// Total number of basis functions
    int nfunction() const;
    /// Total number of functions if this shell was Cartesian
    int ncartesian() const          { return ncartesian_; }
    /// The angular momentum of the given contraction
    int am() const                  { return l_; }
    /// The character symbol for the angular momentum of the given contraction
    char amchar() const             { return amtypes[l_]; }
    /// The character symbol for the angular momentum of the given contraction (upper case)
    char AMCHAR() const             { return AMTYPES[l_]; }
    /// Returns true if contraction is Cartesian
    bool is_cartesian() const       { return !puream_; }
    /// Returns true if contraction is pure
    bool is_pure() const            { return puream_; }

    /// Returns the center of the Molecule this shell is on
    const Vector3& center() const;
    /// Returns the atom number this shell is on. Used by integral derivatives for indexing.
    int ncenter() const             { return nc_; }

    /// Returns the exponent of the given primitive
    double exp(int prim) const      { return exp_[prim]; }
    /// Return coefficient of pi'th primitive and ci'th contraction
    double coef(int pi) const       { return coef_[pi]; }
    /// Returns the exponent of the given primitive
    const std::vector<double>& exps() const { return exp_; }
    /// Return coefficient of pi'th primitive and ci'th contraction
    const std::vector<double>& coefs() const { return coef_; }

    /// Print out the shell
    void print(FILE *out) const;

    /// Normalize the angular momentum component
    static double normalize(int l, int m, int n);

    /// Basis function index where this shell starts.
    int function_index() const      { return start_; }
    void set_function_index(int i)  { start_ = i; }
};

}

#endif
