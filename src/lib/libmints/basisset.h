#ifndef _psi_src_lib_libmints_basisset_h_
#define _psi_src_lib_libmints_basisset_h_

/*!
    \defgroup MINTS libmints: Integral library
    \file libmints/basisset.h
    \ingroup MINTS
*/

#include <cstdio>
#include <libchkpt/chkpt.hpp>

#include <libmints/ref.h>
#include <libmints/molecule.h>
#include <libmints/gshell.h>
#include <libmints/sobasis.h>
#include <libmints/integral.h>

extern FILE *outfile;

namespace psi {
    
//! Basis set container class
/*!
    Reads the basis set from a checkpoint file object. Also reads the molecule
    from the checkpoint file storing the information in an internal Molecule class
    which can be accessed using molecule().
*/
class BasisSet
{
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
    SimpleMatrix *simple_mat_uso2ao_;
    double **uso2bf_;
    SimpleMatrix *simple_mat_uso2bf_;
    
    //! Does the loaded basis set contain pure angular momentum functions?
    bool puream_;
    
    //! Array of gaussian shells
    GaussianShell** shells_;
    //! Molecule object.
    Molecule *molecule_;
    //! Symmetry orbital transformation (used in one-electron integrals)
    SOTransform *sotransform_;
    //! Spherical transfromation (used in two-electron integrals)
    std::vector<SphericalTransform> sphericaltransforms_;
    
    //! No default constructor
    BasisSet();
    //! No assignment
    BasisSet& operator=(const BasisSet&);
    
    //! Initialize shells based on information found in checkpoint
    void initialize_shells(psi::Chkpt* chkpt, std::string& basiskey);
    
    //! Initialize shells based on information found in GENBAS file
    void initialize_shells_via_genbas(std::string& genbas_filename, std::string& genbas_basis);
    
public:

    /// Constructor, reads in the basis set from the checkpoint file using basiskey
    BasisSet(psi::Chkpt* chkpt, std::string basiskey = "");
    
    /// Constructor, reads in the basis set from the checkpoint file
    BasisSet(psi::Chkpt* chkpt, std::string genbas_filename, std::string genbas_basis);
    /// Copy constructor, currently errors if used
    BasisSet(const BasisSet&);
    /// Destructor
    ~BasisSet();
    
    /// Total number of primitives
    int nprimitive() const             { return nprimitives_; }
    /// Maximum number of primitives in a shell
    int max_nprimitive() const         { return max_nprimitives_; }
    /// Number of shells
    int nshell() const                 { return nshells_;     }
    /// Number of atomic orbitals
    int nao() const                    { return nao_;         }
    /// Number of basis functions
    int nbf() const                    { return nbf_;         }
    /// Maximum angular momentum
    int max_am() const                 { return max_am_;      }
    /// Spherical harmonics?
    bool has_puream() const            { return puream_;      }
    /// Molecule this basis is for
    Molecule* molecule() const         { return molecule_;    }
    /// Maximum stabilizer index
    int max_stability_index() const    { return max_stability_index_; }
    /// Given a shell what is its first AO function
    int shell_to_function(int i) const { return shell_first_ao_[i]; }
    int shell_to_basis_function(int i) const { return shell_first_basis_function_[i]; }
    
    /// Return the si'th Gaussian shell
    GaussianShell* shell(int si) const;
    
    /// Returns i'th shell's transform
    SOTransformShell* so_transform(int i) { return sotransform_->aoshell(i); }
    
    /// Returns the transformation object for a given angular momentum. Used in ERIs.
    SphericalTransform* spherical_transform(int am) { return &sphericaltransforms_[am]; }
    
    /// Print the basis set
    void print(FILE *out = outfile) const;
    
    /// Returns the uso2ao_ matrix.
    const SimpleMatrix* uso_to_ao() const { return simple_mat_uso2ao_; }

    /// Returns the uso2bf_ matrix.
    const SimpleMatrix* uso_to_bf() const { return simple_mat_uso2bf_; }
    
    /// Returns an empty basis set object
    static BasisSet* zero_basis_set();
};

}

#endif
