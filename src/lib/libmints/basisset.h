#ifndef _psi_src_lib_libmints_basisset_h_
#define _psi_src_lib_libmints_basisset_h_

#include <cstdio>
#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

// Need libint for maximum angular momentum

// These probably should be remove and anything we need to be forward declared.
#include <libint/libint.h>

#include <boost/thread/once.hpp>

namespace psi {

    extern FILE *outfile;

    class Molecule;
    class GaussianShell;

    class Chkpt;
    class BasisSetParser;
    class SOTransformShell;
    class SphericalTransform;
    class SOTransform;
    class Matrix;
    class SimpleMatrix;
    class Vector3;

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
    //! Number of shells.
    int nshell_;
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

    //! Map function number to shell
    std::vector<int> function_to_shell_;

    //! Does the loaded basis set contain pure angular momentum functions?
    bool puream_;

    //! Number of shells per center
    std::vector<int> center_to_nshell_;
    //! For a given center, what is its first shell.
    std::vector<int> center_to_shell_;

    //! Array of gaussian shells
    std::vector<boost::shared_ptr<GaussianShell> > shells_;

    //! Molecule object.
    boost::shared_ptr<Molecule> molecule_;
    /** Symmetry orbital transformation (used in one-electron integrals)
     *  NOTE: This will likely change.
     */
    boost::shared_ptr<SOTransform> sotransform_;
    //! Spherical transfromation (used in two-electron integrals)
    std::vector<SphericalTransform> sphericaltransforms_;

    //! No default constructor
    BasisSet();
    //! No assignment
    BasisSet& operator=(const BasisSet&);

    //! Initialize shells based on information found in checkpoint
    /*! Reads in information from the checkpoint file constructing GaussianShells
        as it goes. If set, basiskey is passed along to libchkpt when reading to
        read in the non-default basis set information.
        @param chkpt Checkpoint library object to read from.
        @param basiskey If reading non-default basis set information then this is set to the suffix of the TOC entries.
      */
//    void initialize_shells(boost::shared_ptr<psi::Chkpt> chkpt, std::string& basiskey);

    // Has static information been initialized?
    static boost::once_flag initialized_shared_;
    // Global arrays of x, y, z exponents
    static std::vector<Vector3> exp_ao[];

public:

    /** Constructor, reads in the basis set from the checkpoint file using basiskey
     *  @param chkpt Checkpoint library object that contains the basis set information.
     *  @param basiskey To load the default basis set leave this parameter empty.
     *                  If an RI-basis is wanted pass "DF"
     */
//    BasisSet(boost::shared_ptr<psi::Chkpt> chkpt, std::string basiskey = "");

    /** Copy constructor, currently errors if used. */
    BasisSet(const BasisSet&);
    /// Destructor
    ~BasisSet();

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
    int nshell() const                 { return nshell_;     }
    /** Number of atomic orbitals (Cartesian).
     * @return The number of atomic orbitals (Cartesian orbitals).
     */
    int nao() const                    { return nao_;         }
    /** Number of basis functions (Spherical).
     *  @return The number of basis functions (Spherical).
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
    /** Molecule this basis is for.
     *  @return Shared pointer to the molecule for this basis set.
     */
    boost::shared_ptr<Molecule> molecule() const         { return molecule_;    }
    /** Given a shell what is its first AO function
     *  @param i Shell number
     *  @return The function number for the first function for the i'th shell.
     */
    int shell_to_function(int i) const { return shell_first_ao_[i]; }
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

    /** Return the si'th Gaussian shell
     *  @param i Shell number
     *  @return A shared pointer to the GaussianShell object for the i'th shell.
     */
    boost::shared_ptr<GaussianShell> shell(int si) const;

    /** Return the i'th Gaussian shell on center
     *  @param i Shell number
     *  @return A shared pointer to the GaussianShell object for the i'th shell.
     */
    boost::shared_ptr<GaussianShell> shell(int center, int si) const;

    /** Returns the transformation object for a given angular momentum. Used in ERIs.
     *  @param am Angular momentum
     *  @return A SphericalTransform object that details how to transfrom from AO to BF.
     */
    SphericalTransform& spherical_transform(int am);

    /** Print the basis set.
     *  @param out The file stream to use for printing. Defaults to outfile.
     */
    void print(FILE *out = outfile) const;

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
    static boost::shared_ptr<BasisSet> zero_basis_set();

    /** Returns a new BasisSet object.
     *
     * Returns a new BasisSet object configured with the provided Molecule object.
     * @param parser The basis set parser object that will be used to interpret the basis set file.
     * @param mol Molecule to construct basis set for.
     * @param basisname Name of the basis set to search for in pbasis.dat
     * @return A new basis set object constructed from the information passed in.
     */
    static boost::shared_ptr<BasisSet> construct(const boost::shared_ptr<BasisSetParser>& parser,
        const boost::shared_ptr<Molecule>& mol,
        const std::string &basisname);

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
        const std::vector<std::string> &basisnames);

    friend class Gaussian94BasisSetParser;
};

}

#endif
