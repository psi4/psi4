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

#ifndef _psi_src_lib_libmints_basisset_h_
#define _psi_src_lib_libmints_basisset_h_

#include <cstdio>
#include <string>
#include <vector>

#include "gshell.h"
#include "molecule.h"

 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
#include "psi4/psi4-dec.h"
namespace psi {



class Molecule;
class GaussianShell;

class Chkpt;
class BasisSetParser;
class DealiasBasisSet;
class SOTransformShell;
class SphericalTransform;
class SOTransform;
class Matrix;
class Vector3;
class SOBasisSet;
class IntegralFactory;

/*! \ingroup MINTS */

//! Basis set container class
/*! Reads the basis set from a checkpoint file object. Also reads the molecule
    from the checkpoint file storing the information in an internal Molecule class
    which can be accessed using molecule().
*/
class BasisSet
{
protected:
    friend class BasisSetParser;

    //! The name of this basis set (e.g. "BASIS", "RI BASIS")
    std::string name_;

    // Key and target information from constructing
    std::string key_;
    std::string target_;

    //! Array of gaussian shells
    GaussianShell *shells_;

    //! vector of shells numbers sorted in acending AM order.
    std::vector<int> sorted_ao_shell_list_;

    //! Molecule object.
    std::shared_ptr<Molecule> molecule_;

    // Has static information been initialized?
    static bool initialized_shared_;

    /*
     * Scalars
     */
    /// Number of atomic orbitals (Cartesian)
    int nao_;
    /// Number of basis functions (either cartesian or spherical)
    int nbf_;
    /// The number of unique primitives
    int n_uprimitive_;
    /// The number of shells
    int n_shells_;
    /// The number of primitives
    int nprimitive_;
    /// The maximum angular momentum
    int max_am_;
    /// The maximum number of primitives in a shell
    int max_nprimitive_;
    /// Whether the basis set is uses spherical basis functions or not
    bool puream_;

    /*
     * Arrays
     */
    /// The number of primitives (and exponents) in each shell
    int *n_prim_per_shell_;
    /// The first (Cartesian) atomic orbital in each shell
    int *shell_first_ao_;
    /// The first (Cartesian / spherical) basis function in each shell
    int *shell_first_basis_function_;
    /// Shell number to atomic center.
    int *shell_center_;
    /// Which shell does a given (Cartesian / spherical) function belong to?
    int *function_to_shell_;
    /// Which shell does a given Cartesian function belong to?
    int *ao_to_shell_;
    /// Which center is a given function on?
    int *function_center_;
    /// How many shells are there on each center?
    int *center_to_nshell_;
    /// What's the first shell on each center?
    int *center_to_shell_;

    /// The flattened lists of unique exponents
    double *uexponents_;
    /// The flattened lists of unique contraction coefficients (normalized)
    double *ucoefficients_;
    /// The flattened lists of unique contraction coefficients (as provided by the user)
    double *uoriginal_coefficients_;
    /// The flattened list of r exponenets for ECP calculations
    int *uns_;
    /// The flattened lists of ERD normalized contraction coefficients
    double *uerd_coefficients_;
    /// The flattened list of Cartesian coordinates for each atom
    double *xyz_;




public:
    BasisSet();

    BasisSet(const std::string &basistype, SharedMolecule mol,
             std::map<std::string, std::map<std::string, std::vector<ShellInfo> > > &shell_map);

    /** Builder factory method
     * @param molecule the molecule to build the BasisSet around
     * @param shells array of *atom-numbered* GaussianShells to build the BasisSet from
     * @return BasisSet corresponding to this molecule and set of shells
     */
    static std::shared_ptr<BasisSet> build(std::shared_ptr<Molecule> molecule,
                                             const std::vector<ShellInfo> &shells);

    /** Initialize singleton values that are shared by all basis set objects. */
    static void initialize_singletons();

    /** Number of primitives.
     *  @return The total number of primitives in all contractions.
     */
    int nprimitive() const             { return nprimitive_; }
    /** Maximum number of primitives in a shell.
     *  Examines each shell and find the shell with the maximum number of primitives returns that
     *  number of primitives.
     *  @return Maximum number of primitives.
     */
    int max_nprimitive() const         { return max_nprimitive_; }
    /** Number of shells.
     *  @return Number of shells.
     */
    int nshell() const                 { return n_shells_;  }
    /** Number of atomic orbitals (Cartesian).
     * @return The number of atomic orbitals (Cartesian orbitals, always).
     */
    int nao() const                    { return nao_;         }
    /** Number of basis functions (Spherical).
     *  @return The number of basis functions (Spherical, if has_puream() == true).
     */
    int nbf() const                    { return nbf_;         }
    /** Maximum angular momentum used in the basis set.
     *  @return Maximum angular momentum.
     */
    int max_am() const                 { return max_am_;      }
    /** Spherical harmonics?
     *  @return true if using spherical harmonics
     */
    bool has_puream() const            { return puream_;      }
    /** Compute the maximum number of basis functions contained in a shell.
     *  @return The max number of basis functions in a shell.
     */
    int max_function_per_shell() const { return (puream_) ? 2*max_am_+1 : (max_am_+1)*(max_am_+2)/2; }
    /** Molecule this basis is for.
     *  @return Shared pointer to the molecule for this basis set.
     */
    std::shared_ptr<Molecule> molecule() const;
    /** Given a shell what is its first AO function
     *  @param i Shell number
     *  @return The function number for the first function for the i'th shell.
     */
    int shell_to_ao_function(int i) const { return shell_first_ao_[i]; }
    /** Given a shell what is its atomic center
     *  @param i Shell number
     *  @return The atomic center for the i'th shell.
     */
    int shell_to_center(int i) const { return shell_center_[i]; }
    /** Given a shell what is its first basis function (spherical) function
     *  @param i Shell number
     *  @return The function number for the first function for the i'th shell.
     */
    int shell_to_basis_function(int i) const { return shell_first_basis_function_[i]; }

    /** Given a function number what shell does it correspond to. */
    int function_to_shell(int i) const { return function_to_shell_[i]; }
    /** Given a function what is its atomic center
     *  @param i Function number
     *  @return The atomic center for the i'th function.
     */
    int function_to_center(int i) const { return function_center_[i]; }

    /** Given a Cartesian function (AO) number what shell does it correspond to. */
    int ao_to_shell(int i) const { return ao_to_shell_[i]; }

    /** Return the si'th Gaussian shell
     *  @param si Shell number
     *  @return A shared pointer to the GaussianShell object for the i'th shell.
     */
    const GaussianShell& shell(int si) const;

    /** Return the i'th Gaussian shell on center
     *  @param center atomic center
     *  @param si Shell number
     *  @return A shared pointer to the GaussianShell object for the i'th shell.
     */
    const GaussianShell& shell(int center, int si) const;

    /** @{
     *  Print the basis set.
     *  @param out The file stream to use for printing. Defaults to outfile.
     */
    void print(std::string out) const;
    void print() const { print("outfile"); }
    /// @}

    /// Returns the name of this basis set
    const std::string & name() const { return name_; }
    void set_name(const std::string str) {name_ = str;}

    /// Return the construction key and target information
    const std::string & key() const { return key_; }
    const std::string & target() const { return target_; }

    /** Print basis set information according to the level of detail in print_level
     *  @param out The file stream to use for printing. Defaults to outfile.
     *  @param print_level  < 1: Nothing
                              1: Brief summary
                              2: Summary and contraction details
                            > 2: Full details
                            Defaults to 2
     */
    void print_by_level(std::string out = "outfile", int print_level = 2) const;
    /** Prints a short string summarizing the basis set
     *  @param out The file stream to use for printing. Defaults to outfile.
     */
    void print_summary(std::string out = "outfile") const;

    /** Prints a detailed PSI3-style summary of the basis (per-atom)
     *  @param out The file stream to use for printing. Defaults to outfile.
     */
    void print_detail(std::string out) const;
    void print_detail() const { print_detail("outfile"); }

    /** Returns a string in CFOUR-style of the basis (per-atom)
     *  Format from http://slater.chemie.uni-mainz.de/cfour/index.php?n=Main.OldFormatOfAnEntryInTheGENBASFile
     */
    std::string print_detail_cfour() const;

    /** Refresh internal basis set data. Useful if someone has pushed to shells_.
     *  Pushing to shells_ happens in the BasisSetParsers, so the parsers will
     *  call refresh().
     */
    void refresh();

    /// Return the number of shells on a given center.
    int nshell_on_center(int i) const { return center_to_nshell_[i]; }
    /// Return the overall shell number
    int shell_on_center(int center, int shell) const { return center_to_shell_[center] + shell; }


    /** Returns an empty basis set object.
     *
     *  Returns a BasisSet object that actually has a single s-function
     *  at the origin with an exponent of 0.0 and contraction of 1.0.
     *  @return A new empty BasisSet object.
     */
    static std::shared_ptr<BasisSet> zero_ao_basis_set();

    /** Returns an empty SO basis set object.
     *
     *  Returns an SOBasis object that actually has a single s-function
     *  at the origin with an exponent of 0.0 and contraction of 1.0.
     *  @return A new empty SOBasis object.
     */
    static std::shared_ptr<SOBasisSet> zero_so_basis_set(const std::shared_ptr<IntegralFactory>& factory);

    /** Returns a shell-labeled test basis set object
     *
     * @param max_am maximum angular momentum to build
     * @return pair containing shell labels and four-center
     * test basis for use in benchmarking
     * See libmints/benchmark.cc for details
     */
    static std::pair<std::vector<std::string>, std::shared_ptr<BasisSet> > test_basis_set(int max_am);


    /** Returns a new basis set object
     * Constructs a basis set from the parsed information
     *
     * @param mol           Psi4 molecule
     * @param py::dict      Python dictionary containing the basis information
     * @param forced_puream Force puream or not
    **/
    static std::shared_ptr<BasisSet> construct_from_pydict(const std::shared_ptr <Molecule> &mol, py::dict pybs, const int forced_puream);

    /** Returns a new basis set object
     * Constructs an ECP basis set from the parsed information
     *
     * @param mol           Psi4 molecule.  WARNING: The nuclear charges are modified by this routine
     * @param py::dict      Python dictionary containing the basis information
     * @param forced_puream Force puream or not
    **/

    static std::shared_ptr<BasisSet> construct_ecp_from_pydict(std::shared_ptr <Molecule> mol, py::dict pybs, const int forced_puream);

    /** Converts basis set name to a compatible filename.
     * @param basisname Basis name
     * @return Compatible file name.
     */
    static std::string make_filename(const std::string& basisname);

    /// Global arrays of x, y, z exponents
    static std::vector<Vector3> exp_ao[];

    //! Returns the value of the sorted shell list.
    int get_ao_sorted_shell(const int &i) { return sorted_ao_shell_list_[i]; }
    //! Returns the vector of sorted shell list.
    std::vector<int> get_ao_sorted_list() { return sorted_ao_shell_list_; }

    // Returns the values of the basis functions at a point
    void compute_phi(double *phi_ao, double x, double y, double z);

    // BasisSet friends
    friend class Gaussian94BasisSetParser;
};

}

#endif
