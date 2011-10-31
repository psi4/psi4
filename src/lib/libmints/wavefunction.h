#ifndef _psi_src_lib_libmints_wavefunction_h
#define _psi_src_lib_libmints_wavefunction_h

#include <stddef.h>

#include "typedefs.h"
#include "dimension.h"
#include <libparallel/parallel.h>

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

#include <vector>

#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif

namespace psi {

class Molecule;
class BasisSet;
class IntegralFactory;
class Matrix;
class Vector;
class MatrixFactory;
class Options;
class SOBasisSet;
class PSIO;
class Chkpt;

/*! \ingroup MINTS
 *  \class Wavefunction
 *  \brief Simple wavefunction base class.
 */
class Wavefunction
{
protected:
    /// Name of the wavefunction
    std::string name_;

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

    /// Integral factory
    boost::shared_ptr<IntegralFactory> integral_;

    /// Matrix factory for creating standard sized matrices
    boost::shared_ptr<MatrixFactory> factory_;

    boost::shared_ptr<Wavefunction> reference_wavefunction_;

    /// How much memory you have access to.
    long int memory_;

    /// Debug flag
    unsigned int debug_;

    /// Whether this wavefunction was obtained using density fitting
    bool density_fitted_;

    /// Energy convergence threshold
    double energy_threshold_;

    /// Density convergence threshold
    double density_threshold_;

    /// Total alpha and beta electrons
    int nalpha_, nbeta_;

    /// Number of doubly occupied per irrep
    Dimension doccpi_;
    /// Number of singly occupied per irrep
    Dimension soccpi_;
    /// Number of frozen core per irrep
    Dimension frzcpi_;
    /// Number of frozen virtuals per irrep
    Dimension frzvpi_;
    /// Number of alpha electrons per irrep
    Dimension nalphapi_;
    /// Number of beta electrons per irrep
    Dimension nbetapi_;

    /// Number of so per irrep
    Dimension nsopi_;
    /// Number of mo per irrep
    Dimension nmopi_;

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

    /// Alpha density matrix
    SharedMatrix Da_;
    /// Beta density matrix
    SharedMatrix Db_;

    /// Lagrangian matrix
    SharedMatrix Lagrangian_;

    /// Alpha Fock matrix
    SharedMatrix Fa_;
    /// Beta Fock matrix
    SharedMatrix Fb_;

    /// Alpha orbital eneriges
    boost::shared_ptr<Vector> epsilon_a_;
    /// Beta orbital energies
    boost::shared_ptr<Vector> epsilon_b_;

    // Callback routines to Python
    std::vector<void*> precallbacks_;
    std::vector<void*> postcallbacks_;

    /// If a gradient is available it will be here:
    SharedMatrix gradient_;

    /// The TPDM contribution to the gradient
    boost::shared_ptr<Matrix> tpdm_gradient_contribution_;

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

    /// An integral factory with basisset() on each center.
    boost::shared_ptr<IntegralFactory> integral() const;
    /// Returns the basis set object that pertains to this wavefunction.
    boost::shared_ptr<BasisSet> basisset() const;
    /// Returns the SO basis set object that pertains to this wavefunction.
    boost::shared_ptr<SOBasisSet> sobasisset() const;
    /// Returns the MatrixFactory object that pertains to this wavefunction
    boost::shared_ptr<MatrixFactory> matrix_factory() const;
    /// Returns the reference wavefunction
    boost::shared_ptr<Wavefunction> reference_wavefunction() const;
    /// Sets the reference wavefunction
    void set_reference_wavefunction(const boost::shared_ptr<Wavefunction> wfn);

    /// Returns whether this wavefunction was obtained using density fitting or not
    bool density_fitted() const { return density_fitted_; }

    static void initialize_singletons();

    /// Returns the DOCC per irrep array. You DO NOT own this array.
    const Dimension& doccpi() const { return doccpi_; }
    /// Returns the SOCC per irrep array. You DO NOT own this array.
    const Dimension& soccpi() const { return soccpi_; }
    /// Returns the number of SOs per irrep array. You DO NOT own this array.
    const Dimension& nsopi() const { return nsopi_; }
    /// Returns the number of MOs per irrep array. You DO NOT own this array.
    const Dimension& nmopi() const { return nmopi_; }
    /// Returns the number of alpha electrons per irrep array. You DO NOT own this array.
    const Dimension& nalphapi() const { return nalphapi_; }
    /// Returns the number of beta electrons per irrep array. You DO NOT own this array.
    const Dimension& nbetapi() const { return nbetapi_; }
    /// Returns the frozen core orbitals per irrep array. You DO NOT own this array.
    const Dimension& frzcpi() const { return frzcpi_; }
    /// Returns the frozen virtual orbitals per irrep array. You DO NOT own this array.
    const Dimension& frzvpi() const { return frzvpi_; }
    /// Returns the number of SOs
    int nso() const { return nso_; }
    /// Returns the number of MOs
    int nmo() const { return nmo_; }
    /// Returns the number of irreps
    int nirrep() const { return nirrep_; }
    /// Returns the reference energy
    double reference_energy () const { return energy_; }

    /// Returns the alpha electrons MO coefficients
    SharedMatrix Ca() const;
    /// Returns the beta electrons MO coefficients
    SharedMatrix Cb() const;
    /// Returns the alpha Fock matrix
    SharedMatrix Fa() const;
    /// Returns the beta Fock matrix
    SharedMatrix Fb() const;
    /// Returns the alpha orbital energies
    boost::shared_ptr<Vector> epsilon_a() const;
    /// Returns the beta orbital energies
    boost::shared_ptr<Vector> epsilon_b() const;
    /// Returns the SO basis Lagrangian
    boost::shared_ptr<Matrix> Lagrangian() const;
    /// The two particle density matrix contribution to the gradient
    virtual boost::shared_ptr<Matrix> tpdm_gradient_contribution() const;

    /// Returns the alpha OPDM for the wavefunction
    const SharedMatrix Da() const;
    /// Returns the beta OPDM for the wavefunction
    SharedMatrix Db() const;

    /// Returns the Lagrangian in SO basis for the wavefunction
    SharedMatrix X() const;

    /// Adds a pre iteration Python callback function
    void add_preiteration_callback(PyObject*);
    /// Adds a post iteration Python callback function
    void add_postiteration_callback(PyObject*);

    /// Call pre iteration callbacks
    void call_preiteration_callbacks();
    /// Call post iteration callbacks
    void call_postiteration_callbacks();

    /// Returns the gradient
    SharedMatrix gradient() const;
    /// Set the gradient for the wavefunction
    void set_gradient(SharedMatrix& grad);

    /// Set the wavefunction name (e.g. "RHF", "ROHF", "UHF", "CCEnergyWavefunction")
    void set_name(const std::string& name) { name_ = name; }

    /// Returns the wavefunction name
    const std::string& name() const { return name_; }

    /// Save the wavefunction to checkpoint
    virtual void save() const;
};

}

#endif
