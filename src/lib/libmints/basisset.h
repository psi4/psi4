#ifndef _psi_src_lib_libmints_basisset_h_
#define _psi_src_lib_libmints_basisset_h_

#include <cstdio>
#include <libchkpt/chkpt.hpp>

#include <libmints/molecule.h>
#include <libmints/gshell.h>
#include <libmints/sobasis.h>
#include <libmints/integral.h>
// Need libint for maximum angular momentum
#include <libint/libint.h>
#include <psi4-dec.h>

#include <boost/thread/once.hpp>

namespace psi {

    class BasisSetParser;
    class SOTransformShell;
    class SphericalTransform;
    class SOTransform;
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
    int nprimitives_;
    //! Number of shells.
    int nshells_;
    //! Number of atomic orbitals.
    int nao_;
    //! Number of basis functions (either cartesian or spherical)
    int nbf_;
    //! Maximum angular momentum
    int max_am_;
    //! Maximum number of primitives.
    int max_nprimitives_;
    //! Shell number to first basis function index.
    int *shell_first_basis_function_;
    //! Shell number to first atomic function index.
    int *shell_first_ao_;
    //! Shell number to atomic center.
    int *shell_center_;
    //! Not used, yet.
    int max_stability_index_;
    //! Unique symmetry orbitals to atomic orbitals.
    double **uso2ao_;
    shared_ptr<SimpleMatrix> simple_mat_uso2ao_;
    double **uso2bf_;
    shared_ptr<SimpleMatrix> simple_mat_uso2bf_;

    //! Does the loaded basis set contain pure angular momentum functions?
    bool puream_;

    //! Array of gaussian shells
    std::vector<shared_ptr<GaussianShell> > shells_;
    //! Molecule object.
    shared_ptr<Molecule> molecule_;
    //! Symmetry orbital transformation (used in one-electron integrals)
    shared_ptr<SOTransform> sotransform_;
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
    void initialize_shells(shared_ptr<psi::Chkpt> chkpt, std::string& basiskey);

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
    BasisSet(shared_ptr<psi::Chkpt> chkpt, std::string basiskey = "");

    /** Copy constructor, currently errors if used. */
    BasisSet(const BasisSet&);
    /// Destructor
    ~BasisSet();

    /** Initialize singleton values that are shared by all basis set objects. */
    static void initialize_singletons();

    /** Number of primitives.
     *  @return The total number of primitives in all contractions.
     */
    int nprimitive() const             { return nprimitives_; }
    /** Maximum number of primitives in a shell.
     *  Examines each shell and find the shell with the maximum number of primitives returns that
     *  number of primitives.
     *  @return Maximum number of primitives.
     */
    int max_nprimitive() const         { return max_nprimitives_; }
    /** Number of shells.
     *  @return Number of shells.
     */
    int nshell() const                 { return nshells_;     }
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
    shared_ptr<Molecule> molecule() const         { return molecule_;    }
    /// Maximum stabilizer index
    int max_stability_index() const    { return max_stability_index_; }
    /** Given a shell what is its first AO function
     *  @param i Shell number
     *  @return The function number for the first function for the i'th shell.
     */
    int shell_to_function(int i) const { return shell_first_ao_[i]; }
    /** Given a shell what is its first basis function (spherical) function
     *  @param i Shell number
     *  @return The function number for the first function for the i'th shell.
     */
    int shell_to_basis_function(int i) const { return shell_first_basis_function_[i]; }

    /** Return the si'th Gaussian shell
     *  @param i Shell number
     *  @return A shared pointer to the GaussianShell object for the i'th shell.
     */
    shared_ptr<GaussianShell> shell(int si) const;

    /** Returns i'th shell's transform object.
     *  @param i Shell number
     *  @return A SOTransformShell object that details how to transform from AO to SO.
     */
    SOTransformShell* so_transform(int i) { return sotransform_->aoshell(i); }

    /** Returns the transformation object for a given angular momentum. Used in ERIs.
     *  @param am Angular momentum
     *  @return A SphericalTransform object that details how to transfrom from AO to BF.
     */
    SphericalTransform& spherical_transform(int am) { return sphericaltransforms_[am]; }

    /** Print the basis set.
     *  @param out The file stream to use for printing. Defaults to outfile.
     */
    void print(FILE *out = outfile) const;

    /** Returns the uso2ao_ matrix.
     *  @return The transformation matrix for USO to AO.
     */
    const shared_ptr<SimpleMatrix> uso_to_ao() const { return simple_mat_uso2ao_; }

    /** Returns the uso2bf_ matrix.
     *  @return The transformation matrix for USO to BF.
     */
    const shared_ptr<SimpleMatrix> uso_to_bf() const { return simple_mat_uso2bf_; }

    /** Refresh internal basis set data. Useful if someone has pushed to shells_.
     *  Pushing to shells_ happens in the BasisSetParsers, so the parsers will
     *  call refresh().
     */
    void refresh();

    /** Return a BasisSet object containing all shells at center i
    *
    * Used for Atomic HF computations for SAD Guesses
    */
    shared_ptr<BasisSet> atomic_basis_set(int center);    

    /** Returns an empty basis set object.
     *
     *  Returns a BasisSet object that actually has a single s-function
     *  at the origin with an exponent of 0.0 and contraction of 1.0.
     *  @return A new empty BasisSet object.
     */
    static shared_ptr<BasisSet> zero_basis_set();

    /** Returns a new BasisSet object.
     *
     * Returns a new BasisSet object configured with the provided Molecule object.
     * @param parser The basis set parser object that will be used to interpret the basis set file.
     * @param mol Molecule to construct basis set for.
     * @param basisname Name of the basis set to search for in pbasis.dat
     * @return A new basis set object constructed from the information passed in.
     */
    static shared_ptr<BasisSet> construct(const shared_ptr<BasisSetParser>& parser,
        const shared_ptr<Molecule>& mol,
        const std::string &basisname);

    /** Returns a new BasisSet object.
     *
     * Returns a new BasisSet object configured with the provided Molecule object.
     * @param parser The basis set parser object that will be used to interpret the basis set file.
     * @param mol Molecule to construct basis set for.
     * @param basisnames Name of the basis set for each atom in molecule to search for in pbasis.dat
     * @return A new basis set object constructed from the information passed in.
     */
    static shared_ptr<BasisSet> construct(const shared_ptr<BasisSetParser>& parser,
        const shared_ptr<Molecule>& mol,
        const std::vector<std::string> &basisnames);

    friend class Gaussian94BasisSetParser;
};

}

#endif
