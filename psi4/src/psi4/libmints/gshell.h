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

#ifndef _psi_src_lib_libmints_gshell_h_
#define _psi_src_lib_libmints_gshell_h_

#include <cstdio>
#include <vector>
#include "psi4/libmints/vector3.h"

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

enum ShellType {
    Gaussian = 0,
    ECPType1 = 1,
    ECPType2 = 2
};

/*! \ingroup MINTS
 *  \class ShellInfo
 *  \brief This class has the same behavior as GaussianShell, but implements everything using
 *         slower data structures, which are easier to construct. These are used to build the
 *         basis set, which builds more efficient pointer-based GaussianShell objects.
 */
class ShellInfo
{
protected:
    /// Angular momentum
    int l_;
    /// Flag for pure angular momentum
    int puream_;
    /// Exponents (of length nprimitives_)
    std::vector<double> exp_;
    /// Contraction coefficients (of length nprimitives_)
    std::vector<double> coef_;
    /// R exponents, for ECP basis sets (of length nprimitives_)
    std::vector<int> n_;
    /// ERD normalized contraction coefficients (of length nprimitives_)
    std::vector<double> erd_coef_;
    /// Original (un-normalized) contraction coefficients (of length nprimitives)
    /// Only used in printing.
    std::vector<double> original_coef_;
    /// The type of shell this object encodes
    ShellType shelltype_;
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
     */
    void contraction_normalization();

    /** Lookup array that when you index the angular momentum it returns the lowercase letter corresponding to it. */
    static const char *amtypes;
    /** Lookup array that when you index the angular momentum it returns the uppercase letter corresponding to it. */
    static const char *AMTYPES;

public:
    /** Constructor; use this version for ECP basis sets.
     *  @param am Angular momentum.
     *  @param e An array of exponent values.
     *  @param n An array of r exponents for ECPs.
     *  @param c An array of contraction coefficients.
     *  @param nc The atomic center that this shell is located on. Must map back to the correct atom in the owning BasisSet molecule_. Used in integral derivatives for indexing.
     *  @param center The x, y, z position of the shell. This is passed to reduce the number of calls to the molecule.
     *  @param start The starting index of the first function this shell provides. Used to provide starting positions in matrices.
     */
    ShellInfo(int am,
              const std::vector<double>& c,
              const std::vector<double>& e,
              const std::vector<int>& n,
              int nc,
              const Vector3& center,
              int start);

    /** Constructor; use this version for regular Gaussian basis sets.
     *  @param e An array of exponent values.
     *  @param am Angular momentum.
     *  @param pure Pure spherical harmonics, or Cartesian.
     *  @param c An array of contraction coefficients.
     *  @param nc The atomic center that this shell is located on. Must map back to the correct atom in the owning BasisSet molecule_. Used in integral derivatives for indexing.
     *  @param center The x, y, z position of the shell. This is passed to reduce the number of calls to the molecule.
     *  @param start The starting index of the first function this shell provides. Used to provide starting positions in matrices.
     *  @param pt Is the shell already normalized?
     */
    ShellInfo(int am,
              const std::vector<double>& c,
              const std::vector<double>& e,
              GaussianType pure,
              int nc,
              const Vector3& center,
              int start,
              PrimitiveType pt = Normalized);

    /** Handles calling primitive_normalization and contraction_normalization for you. */
    void normalize_shell();
    void erd_normalize_shell();

    /// Make a copy of the ShellInfo.
    ShellInfo copy();

    /// Make a copy of the ShellInfo.
    ShellInfo copy(int nc, const Vector3& c);

    /// The type of shell this object encodes information for.
    ShellType shell_type() const { return shelltype_; }
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
    /// Return coefficient of pi'th primitive
    double coef(int pi) const       { return coef_[pi]; }
    /// Return r exponent for pi'th ECP primitive
    int nval(int pi) const          { return n_[pi]; }
    /// Return ERD normalized coefficient of pi'th primitive
    double erd_coef(int pi) const       { return erd_coef_[pi]; }
    /// Return unnormalized coefficient of pi'th primitive
    double original_coef(int pi) const { return original_coef_[pi]; }
    /// Returns the exponent of the given primitive
    const std::vector<double>& exps() const { return exp_; }
    /// Return coefficient of pi'th primitive and ci'th contraction
    const std::vector<double>& coefs() const { return coef_; }
    /// Return unnormalized coefficient of pi'th primitive and ci'th contraction
    const std::vector<double>& original_coefs() const { return original_coef_; }

    /// Print out the shell
    void print(std::string out) const;

    /// Normalize the angular momentum component
    static double normalize(int l, int m, int n);

    /// Basis function index where this shell starts.
    int function_index() const      { return start_; }
    void set_function_index(int i)  { start_ = i; }

    ///True if two ShellInfo instances are exactly equal
    bool operator==(const ShellInfo& RHS)const;

    ///True if any part of the two ShellInfo instances is unequal
    bool operator!=(const ShellInfo& RHS)const
    {
        return !(*this==RHS);
    }
};


/*! \ingroup MINTS
 *  \class GaussianShell
 *  \brief Gaussian orbital shell.
 */
class GaussianShell
{
protected:
    /// Angular momentum
    int l_;
    /// Flag for pure angular momentum
    int puream_;
    /// Exponents (of length nprimitives_)
    const double* exp_;
    /// Original (un-normalized) contraction coefficients (of length nprimitives)
    const double* original_coef_;
    /// Contraction coefficients (of length nprimitives_)
    const double* coef_;
    /// Contraction coefficients normalized for the ERD integral package (of length nprimitives_)
    const double* erd_coef_;
    /// R exponents for ECPs (of length nprimitives_)
    const int* n_;

    /// The type of shell this structure encodes.
    ShellType shelltype_;
    /// Atom number this shell goes to. Needed when indexing integral derivatives.
    int nc_;
    /// Atomic coordinates of this center
    const double *center_;
    /// First basis function in this shell
    int start_;

    /// How many cartesian functions? (1=s, 3=p, 6=d, ...)
    int ncartesian_;
    /// The number of primitives in this shell
    int nprimitive_;
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
     */
    void contraction_normalization();

    /** Lookup array that when you index the angular momentum it returns the lowercase letter corresponding to it. */
    static const char *amtypes;
    /** Lookup array that when you index the angular momentum it returns the uppercase letter corresponding to it. */
    static const char *AMTYPES;

public:
    /** Constructor; Use this version for regular Gaussian basis sets.
     *  @param shelltype The type of shell this structure describes
     *  @param am Angular momentum.
     *  @param pure Pure spherical harmonics, or Cartesian.
     *  @param oc An array of contraction coefficients.
     *  @param c An array of normalized contraction coefficients.
     *  @param ec An array of ERD normalized contraction coefficients.
     *  @param e An array of exponent values.
     *  @param nc The atomic center that this shell is located on. Must map back to the correct atom in the owning BasisSet molecule_. Used in integral derivatives for indexing.
     *  @param center The x, y, z position of the shell. This is passed to reduce the number of calls to the molecule.
     *  @param start The starting index of the first function this shell provides. Used to provide starting positions in matrices.
     */
    GaussianShell(ShellType shelltype,
                  int am,
                  int nprimitive,
                  const double *oc,
                  const double *c, const double *ec,
                  const double *e,
                  GaussianType pure,
                  int nc,
                  const double* center,
                  int start);

    /** Constructor; use this version for ECPs.
     *  @param shelltype The type of shell this structure describes
     *  @param am Angular momentum.
     *  @param pure Pure spherical harmonics, or Cartesian.
     *  @param oc An array of contraction coefficients.
     *  @param e An array of exponent values.
     *  @param n An array of r exponent values for ECPs.
     *  @param nc The atomic center that this shell is located on. Must map back to the correct atom in the owning BasisSet molecule_. Used in integral derivatives for indexing.
     *  @param center The x, y, z position of the shell. This is passed to reduce the number of calls to the molecule.
     *  @param start The starting index of the first function this shell provides. Used to provide starting positions in matrices.
     */
    GaussianShell(ShellType shelltype,
                  int am,
                  int nprimitive,
                  const double *oc,
                  const double *e,
                  const int *n,
                  GaussianType pure,
                  int nc,
                  const double* center,
                  int start);
    ///Builds and empty GShell
    GaussianShell() {}

    /// The type of shell this object encodes information for.
    ShellType shell_type() const { return shelltype_; }
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
    const double* center() const;
    /// Returns the atom number this shell is on. Used by integral derivatives for indexing.
    int ncenter() const             { return nc_; }
    /// Returns the first basis function in this shell
    int start() const { return start_; }

    /// Returns the exponent of the given primitive
    double exp(int prim) const      { return exp_[prim]; }
    /// Return coefficient of pi'th primitive
    double coef(int pi) const       { return coef_[pi]; }
    /// Return unnormalized coefficient of pi'th primitive
    double original_coef(int pi) const { return original_coef_[pi]; }
    /// Return the pi'th r exponent for ECPs
    int nval(int pi) const          { return n_[pi]; }
    /// Return unnormalized coefficient of pi'th primitive
    double erd_coef(int pi) const { return erd_coef_[pi]; }
    /// Returns the exponents
    const double* exps() const { return exp_; }
    /// Return coefficients
    const double* coefs() const { return coef_; }
    /// Return unnormalized coefficients
    const double* original_coefs() const { return original_coef_; }
    /// Return ERD normalized coefficients
    const double* erd_coefs() const { return erd_coef_; }

    /// Print out the shell
    void print(std::string out) const;

    /// Basis function index where this shell starts.
    int function_index() const      { return start_; }
    void set_function_index(int i)  { start_ = i; }
    double evaluate(double r, int l) const;
};

}

#endif
