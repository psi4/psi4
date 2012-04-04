#ifndef _psi_src_lib_libmints_basisset_h_
#define _psi_src_lib_libmints_basisset_h_

#include <cstdio>
#include <string>
#include <vector>

// Need libint for maximum angular momentum
#include <libint/libint.h>

#include "gshell.h"
#include "molecule.h"

#include <boost/shared_ptr.hpp>

namespace boost {
class once_flag;
}

namespace psi {

    extern FILE *outfile;

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
    friend class BasisSetParser;
    //! Number of primitives.
    int nprimitive_;
    //! Number of atomic orbitals.
    int nao_;
    //! Number of basis functions (either cartesian or spherical)
    int nbf_;
    //! Maximum angular momentum
    int max_am_;
    //! Maximum number of primitives.
    int max_nprimitive_;
    //! Shell number to first basis function index.
    std::vector<int> shell_first_basis_function_;           // Is this used?
    //! Shell number to first atomic function index.
    std::vector<int> shell_first_ao_;
    //! Shell number to atomic center.
    std::vector<int> shell_center_;
    //! Function number to atomic center.
    std::vector<int> function_center_;

    //! Map function number to shell
    std::vector<int> function_to_shell_;
    //! Map Cartesian function number to shell
    std::vector<int> ao_to_shell_;

    //! Does the loaded basis set contain pure angular momentum functions?
    bool puream_;

    //! The name of this basis set (e.g. "BASIS", "RI BASIS")
    std::string name_;
    //! Number of shells per center
    std::vector<int> center_to_nshell_;
    //! For a given center, its first shell.
    std::vector<int> center_to_shell_;

    //! Array of gaussian shells
    std::vector<GaussianShell> shells_;

    //! vector of shells numbers sorted in acending AM order.
    std::vector<int> sorted_ao_shell_list_;

    //! Molecule object.
    boost::shared_ptr<Molecule> molecule_;

    // Has static information been initialized?
    static boost::once_flag initialized_shared_;

public:
    BasisSet();

    /** Builder factory method
     * @param molecule the molecule to build the BasisSet around
     * @param shells array of *atom-numbered* GaussianShells to build the BasisSet from
     * @return BasisSet corresponding to this molecule and set of shells
     */
    static boost::shared_ptr<BasisSet> build(boost::shared_ptr<Molecule> molecule,
                                             const std::vector<GaussianShell>& shells);

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
    int nshell() const                 { return shells_.size();  }
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
    boost::shared_ptr<Molecule> molecule() const;
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
     *  @param i Shell number
     *  @return A shared pointer to the GaussianShell object for the i'th shell.
     */
    const GaussianShell& shell(int si) const;

    /** Return the i'th Gaussian shell on center
     *  @param i Shell number
     *  @return A shared pointer to the GaussianShell object for the i'th shell.
     */
    const GaussianShell& shell(int center, int si) const;

    /** @{
     *  Print the basis set.
     *  @param out The file stream to use for printing. Defaults to outfile.
     */
    void print(FILE *out) const;
    void print() const { print(outfile); }
    /// @}

    /// Returns the name of this basis set
    const std::string & name() const { return name_; }
    void set_name(const std::string str) {name_ = str;}

    /** Print basis set information according to the level of detail in print_level
     *  @param out The file stream to use for printing. Defaults to outfile.
     *  @param print_level: < 1: Nothing
                              1: Brief summary
                              2: Summary and contraction details
                            > 2: Full details
                            Defaults to 2
     */
    void print_by_level(FILE* out = outfile, int print_level = 2) const;
    /** Prints a short string summarizing the basis set
     *  @param out The file stream to use for printing. Defaults to outfile.
     */
    void print_summary(FILE *out = outfile) const;

    /** Prints a detailed PSI3-style summary of the basis (per-atom)
     *  @param out The file stream to use for printing. Defaults to outfile.
     */
    void print_detail(FILE *out) const;
    void print_detail() const { print_detail(outfile); }

    /** Refresh internal basis set data. Useful if someone has pushed to shells_.
     *  Pushing to shells_ happens in the BasisSetParsers, so the parsers will
     *  call refresh().
     */
    void refresh();

    /// Return the number of shells on a given center.
    int nshell_on_center(int i) const { return center_to_nshell_[i]; }
    /// Return the overall shell number
    int shell_on_center(int center, int shell) const { return center_to_shell_[center] + shell; }

    /** Return a BasisSet object containing all shells at center i
     *
     * Used for Atomic HF computations for SAD Guesses
     *
     * @param center Atomic center to provide a basis object for.
     * @returns A new basis set object for the atomic center.
     */
    boost::shared_ptr<BasisSet> atomic_basis_set(int center);

    /** Returns an empty basis set object.
     *
     *  Returns a BasisSet object that actually has a single s-function
     *  at the origin with an exponent of 0.0 and contraction of 1.0.
     *  @return A new empty BasisSet object.
     */
    static boost::shared_ptr<BasisSet> zero_ao_basis_set();

    /** Returns an empty SO basis set object.
     *
     *  Returns an SOBasis object that actually has a single s-function
     *  at the origin with an exponent of 0.0 and contraction of 1.0.
     *  @return A new empty SOBasis object.
     */
//    static boost::shared_ptr<SOBasisSet> zero_so_basis_set(const boost::shared_ptr<IntegralFactory>& factory);

    /** Returns a shell-labeled test basis set object
     *
     * @param max_am maximum angular momentum to build
     * @return pair containing shell labels and four-center
     * test basis for use in benchmarking
     * See libmints/benchmark.cc for details
     */
    static std::pair<std::vector<std::string>, boost::shared_ptr<BasisSet> > test_basis_set(int max_am);

    /** Returns a new BasisSet object.
     *
     * Returns a new BasisSet object configured with the provided Molecule object.
     * @param parser The basis set parser object that will be used to interpret the basis set file.
     * @param mol Molecule to construct basis set for.
     * @param basisnames Name of the basis set for each atom in molecule to search for in pbasis.dat
     * @return A new basis set object constructed from the information passed in.
     */
    static boost::shared_ptr<BasisSet> construct(const boost::shared_ptr<BasisSetParser>& parser,
        const boost::shared_ptr<Molecule>& mol,
        const std::string& type);

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

    // BasisSet friends
    friend class Gaussian94BasisSetParser;
    friend BasisSet operator +(const BasisSet& a, const BasisSet& b);
    friend boost::shared_ptr<BasisSet> operator +(const boost::shared_ptr<BasisSet>& a, const boost::shared_ptr<BasisSet>& b);
};

inline
bool shell_sorter_ncenter(const GaussianShell& d1, const GaussianShell& d2)
{
    return d1.ncenter() < d2.ncenter();
}

inline
bool shell_sorter_am(const GaussianShell& d1, const GaussianShell& d2)
{
    return d1.am() < d2.am();
}

inline
BasisSet operator +(const BasisSet& a, const BasisSet& b) {
    if (a.molecule() != b.molecule()) {
        fprintf(stderr, "BasisSet::operator+ : Unable to add basis sets from different molecules.");
        return BasisSet();
    }
    BasisSet temp;

    temp.name_ = a.name_ + " + " + b.name_;
    temp.molecule_ = a.molecule();

    // Copy a's shells to temp
    temp.shells_ = a.shells_;

    // Append b's shells to temp
    temp.shells_.insert(temp.shells_.end(), b.shells_.begin(), b.shells_.end());

    // Sort by center number
    std::sort(temp.shells_.begin(), temp.shells_.end(), shell_sorter_ncenter);

    // Call refresh to regenerate center_to_shell and center_to_nshell
    temp.refresh();

    // Sort by AM in each center
    for (int atom=0; atom < temp.molecule_->natom(); ++atom) {
        std::sort(temp.shells_.begin()+temp.center_to_shell_[atom],
                  temp.shells_.begin()+temp.center_to_shell_[atom]+temp.center_to_nshell_[atom],
                  shell_sorter_am);
    }

    temp.refresh();

    return temp;
}

inline
boost::shared_ptr<BasisSet> operator +(const boost::shared_ptr<BasisSet>& a, const boost::shared_ptr<BasisSet>& b) {
    return boost::shared_ptr<BasisSet>(new BasisSet(*a.get() + *b.get()));
}

}

#endif
