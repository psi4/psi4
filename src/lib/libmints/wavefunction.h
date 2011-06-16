#ifndef _psi_src_lib_libmints_wavefunction_h
#define _psi_src_lib_libmints_wavefunction_h

#include <stddef.h>

#define MAX_IOFF 30000
extern size_t ioff[MAX_IOFF];

#define MAX_DF 500
extern double df[MAX_DF];

#define MAX_BC 20
extern double bc[MAX_BC][MAX_BC];

#define MAX_FAC 100
extern double fac[MAX_FAC];

#define INDEX2(i, j) ( (i) >= (j) ? ioff[(i)] + (j) : ioff[(j)] + (i) )
#define INDEX4(i, j, k, l) ( INDEX2( INDEX2((i), (j)), INDEX2((k), (l)) ) )

#include "factory.h"

#include <vector>

#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif

namespace psi {

class Molecule;
class BasisSet;
class MatrixFactory;
class Options;
class SOBasisSet;

/*! \ingroup MINTS
 *  \class Wavefunction
 *  \brief Simple wavefunction base class.
 */
class Wavefunction {
protected:

    /// Primary basis set for AO integrals
    boost::shared_ptr<BasisSet> basisset_;

    /// Primary basis set for SO integrals
    boost::shared_ptr<SOBasisSet> sobasisset_;

    /// Molecule that this wavefunction is run on
    boost::shared_ptr<Molecule> molecule_;

    /// Options object
    Options & options_;

    // PSI file access variables
    boost::shared_ptr<PSIO> psio_;
    boost::shared_ptr<Chkpt> chkpt_;

    /// Matrix factory for creating standard sized matrices
    boost::shared_ptr<MatrixFactory> factory_;

    boost::shared_ptr<Wavefunction> reference_wavefunction_;

    /// How much memory you have access to.
    long int memory_;

    /// Debug flag
    unsigned int debug_;

    /// Energy convergence threshold
    double energy_threshold_;

    /// Density convergence threshold
    double density_threshold_;

    /// Total alpha and beta electrons per irrep
    int nalpha_, nbeta_;

    /// Number of doubly occupied per irrep
    int doccpi_[8];
    /// Number of singly occupied per irrep
    int soccpi_[8];
    /// Number of frozen core per irrep
    int frzcpi_[8];
    /// Number of frozen virtuals per irrep
    int frzvpi_[8];
    /// Number of alpha electrons per irrep
    int nalphapi_[8];
    /// Number of beta electrons per irrep
    int nbetapi_[8];

    /// Number of so per irrep
    int nsopi_[8];
    /// Number of mo per irrep
    int nmopi_[8];

    /// The energy associated with this wavefunction
    double energy_;

    /// Total number of SOs
    int nso_;
    /// Total number of MOs
    int nmo_;
    /// Number of irreps
    int nirrep_;

    /// Alpha MO coefficients
    SharedMatrix Ca_;
    /// Beta MO coefficients
    SharedMatrix Cb_;

    SharedMatrix Da_;
    SharedMatrix Db_;

    /// Alpha Fock matrix
    SharedMatrix Fa_;
    /// Beta Fock matrix
    SharedMatrix Fb_;

    /// Alpha orbital eneriges
    SharedVector epsilon_a_;
    /// Beta orbital energies
    SharedVector epsilon_b_;

    // Callback routines to Python
    std::vector<void*> precallbacks_;
    std::vector<void*> postcallbacks_;

    /// If a gradient is available it will be here:
    SharedMatrix gradient_;

private:
    // Wavefunction() {}
    void common_init();

public:
    /// Set the PSIO object.
    Wavefunction(Options & options, boost::shared_ptr<PSIO> psio);
    Wavefunction(Options & options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);

    virtual ~Wavefunction();

    /// Compute energy. Subclasses override this function to compute its energy.
    virtual double compute_energy() = 0;

    /// Initialize internal variables from checkpoint file.
    void init_with_chkpt();

    /// Is this a restricted wavefunction?
    virtual bool restricted() const { return true; }

    /// Returns the molecule object that pertains to this wavefunction.
    boost::shared_ptr<Molecule> molecule() const;
    boost::shared_ptr<PSIO> psio() const;
    Options& options() const;

    /// Returns the basis set object that pertains to this wavefunction.
    boost::shared_ptr<BasisSet> basisset() const;
    /// Returns the SO basis set object that pertains to this wavefunction.
    boost::shared_ptr<SOBasisSet> sobasisset() const;
    /// Returns the MatrixFactory object that pertains to this wavefunction
    boost::shared_ptr<MatrixFactory> matrix_factory() const;
    /// Returns the reference wavefunction
    boost::shared_ptr<Wavefunction> reference_wavefunction() const;
    /// Returns the reference wavefunction
    void set_reference_wavefunction(const boost::shared_ptr<Wavefunction> wfn);

    static void initialize_singletons();

    /// Returns the DOCC per irrep array. You DO NOT own this array.
    int* doccpi() const { return (int*)doccpi_; }
    /// Returns the SOCC per irrep array. You DO NOT own this array.
    int* soccpi() const { return (int*)soccpi_; }
    /// Returns the number of SOs per irrep array. You DO NOT own this array.
    int* nsopi() const { return (int*)nsopi_; }
    /// Returns the number of MOs per irrep array. You DO NOT own this array.
    int* nmopi() const { return (int*)nmopi_; }
    /// Returns the number of alpha electrons per irrep array. You DO NOT own this array.
    int* nalphapi() const { return (int*)nalphapi_; }
    /// Returns the number of beta electrons per irrep array. You DO NOT own this array.
    int* nbetapi() const { return (int*)nbetapi_; }
    /// Returns the frozen core orbitals per irrep array. You DO NOT own this array.
    int* frzcpi() const { return (int*)frzcpi_; }
    /// Returns the frozen virtual orbitals per irrep array. You DO NOT own this array.
    int* frzvpi() const { return (int*)frzvpi_; }
    /// Returns the number of SOs
    int nso() const { return nso_; }
    /// Returns the number of MOs
    int nmo() const { return nmo_; }
    /// Returns the number of irreps
    int nirrep() const { return nirrep_; }
    /// Returns the reference energy
    double reference_energy () const { return energy_; }

    /// Returns the alpha electrons MO coefficients
    SharedMatrix Ca() const { return Ca_; }
    /// Returns the beta electrons MO coefficients
    SharedMatrix Cb() const { return Cb_; }
    /// Returns the alpha Fock matrix
    SharedMatrix Fa() const { return Fa_; }
    /// Returns the beta Fock matrix
    SharedMatrix Fb() const { return Fb_; }
    /// Returns the alpha orbital energies
    SharedVector epsilon_a() const { return epsilon_a_; }
    /// Returns the beta orbital energies
    SharedVector epsilon_b() const { return epsilon_b_; }

    /// Returns the alpha OPDM for the wavefunction
    SharedMatrix Da() const { return Da_; }
    /// Returns the beta OPDM for the wavefunction
    SharedMatrix Db() const { return Db_; }

    /// Adds a pre iteration Python callback function
    void add_preiteration_callback(PyObject*);
    /// Adds a post iteration Python callback function
    void add_postiteration_callback(PyObject*);

    /// Call pre iteration callbacks
    void call_preiteration_callbacks();
    /// Call post iteration callbacks
    void call_postiteration_callbacks();

    /// Returns the gradient
    SharedMatrix gradient() const { return gradient_; }
    /// Set the gradient for the wavefunction
    void set_gradient(SharedMatrix& grad) { gradient_ = grad; }
};

}

#endif
