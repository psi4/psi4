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

#ifndef _psi_src_lib_libmints_wavefunction_h
#define _psi_src_lib_libmints_wavefunction_h

#include "typedefs.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/libmints/dimension.h"

#include "psi4/pybind11.h"
#include <stddef.h>
#include <vector>
#include <memory>

#define MAX_IOFF 30000
extern size_t ioff[MAX_IOFF];

#define MAX_DF 500
extern double df[MAX_DF];

#define MAX_BC 20
extern double bc[MAX_BC][MAX_BC];

#define MAX_FAC 100
extern double fac[MAX_FAC];

#if !defined( EXPLICIT_IOFF )
#   define EXPLICIT_IOFF(i) ( (i) * ((i) + 1) / 2 )
#endif

#if !defined( INDEX2 )
#   define INDEX2(i, j) ( (i) >= (j) ? EXPLICIT_IOFF(i) + (j) : EXPLICIT_IOFF(j) + (i) )
#endif

#if !defined( INDEX4 )
#   define INDEX4(i, j, k, l) ( INDEX2( INDEX2((i), (j)), INDEX2((k), (l)) ) )
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
class OrbitalSpace;
class OEProp;

/*! \ingroup MINTS
 *  \class Wavefunction
 *  \brief Simple wavefunction base class.
 */
class Wavefunction : public std::enable_shared_from_this<Wavefunction>
{
protected:
    /// Name of the wavefunction
    std::string name_;

    /// DF/RI/F12/etc basis sets
    std::map<std::string, std::shared_ptr<BasisSet>> basissets_;

    /// The ORBITAL basis
    std::shared_ptr<BasisSet> basisset_;

    /// The ECP basis set
    std::shared_ptr<BasisSet> ecpbasisset_;

    /// Primary basis set for SO integrals
    std::shared_ptr<SOBasisSet> sobasisset_;

    /// AO2SO conversion matrix (AO in rows, SO in cols)
    SharedMatrix AO2SO_;

    /// Molecule that this wavefunction is run on
    std::shared_ptr<Molecule> molecule_;

    /// Options object
    Options & options_;

    // PSI file access variables
    std::shared_ptr<PSIO> psio_;

    /// Integral factory
    std::shared_ptr<IntegralFactory> integral_;

    /// Matrix factory for creating standard sized matrices
    std::shared_ptr<MatrixFactory> factory_;

    std::shared_ptr<Wavefunction> reference_wavefunction_;

    std::shared_ptr<OEProp> oeprop_;

    /// How much memory you have access to.
    long int memory_;

    /// Debug flag
    unsigned int debug_;
    /// Print flag
    unsigned int print_;

    /// Total alpha and beta electrons
    int nalpha_, nbeta_;

    /// Total frozen core orbitals
    int nfrzc_;

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

    /// Whether this wavefunction was obtained using density fitting
    bool density_fitted_;

    /// The energy associated with this wavefunction
    double energy_;

    /// Frozen-core energy associated with this wavefunction
    double efzc_;

    /// Total number of SOs
    int nso_;
    /// Total number of MOs
    int nmo_;
    /// Number of irreps
    int nirrep_;

    /// Overlap matrix
    SharedMatrix S_;

    /// Core Hamiltonian matrix
    SharedMatrix H_;
    SharedMatrix Horig_;

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
    std::shared_ptr<Vector> epsilon_a_;
    /// Beta orbital energies
    std::shared_ptr<Vector> epsilon_b_;

    // Callback routines to Python
    std::vector<void*> precallbacks_;
    std::vector<void*> postcallbacks_;

    /// If a gradient is available it will be here:
    SharedMatrix gradient_;

    /// If a Hessian is available it will be here:
    SharedMatrix hessian_;

    /// The TPDM contribution to the gradient
    std::shared_ptr<Matrix> tpdm_gradient_contribution_;

    /// Helpers for C/D/epsilon transformers
    SharedMatrix C_subset_helper(SharedMatrix C, const Dimension& noccpi, SharedVector epsilon, const std::string& basis, const std::string& subset);
    SharedMatrix F_subset_helper(SharedMatrix F, SharedMatrix C, const std::string& basis);
    SharedVector epsilon_subset_helper(SharedVector epsilon, const Dimension& noccpi, const std::string& basis, const std::string& subset);
    std::vector<std::vector<int> > subset_occupation(const Dimension& noccpi, const std::string& subset);

    /// If atomic point charges are available they will be here
    std::shared_ptr<std::vector<double>> atomic_point_charges_;

    /// If frequencies are available, they will be here:
    std::shared_ptr<Vector> frequencies_;

    /// If normal modes are available, they will be here:
    std::shared_ptr<Vector> normalmodes_;

    /// Same orbs or dens
    bool same_a_b_dens_;
    bool same_a_b_orbs_;

    // Collection of variables
    std::map<std::string, double> variables_;
    std::map<std::string, SharedMatrix> arrays_;

private:
    // Wavefunction() {}
    void common_init();

public:

    /// Constructor for an entirely new wavefunction with an existing basis
    Wavefunction(std::shared_ptr<Molecule> molecule,
                 std::shared_ptr<BasisSet> basis,
                 Options& options);

    /// Constructor for an entirely new wavefunction with an existing basis
    Wavefunction(std::shared_ptr<Molecule> molecule,
                 std::shared_ptr<BasisSet> basis,
                 std::shared_ptr<BasisSet> ecpbasis);

    /// Constructor for an entirely new wavefunction with an existing basis and global options
    Wavefunction(std::shared_ptr<Molecule> molecule,
                 std::shared_ptr<BasisSet> basis);

    /// Blank constructor for derived classes
    Wavefunction(Options & options);

    /**
    * Copy the contents of another Wavefunction into this one.
    * Useful at the beginning of correlated wavefunction computations.
    * -Does not set options, callbacks, or reference_wavefunction_
    * -Matrices and Vectors (Ca,Da,Fa,epsilon_a, etc) are copied by reference,
    *  so if you change these, you must reallocate to avoid compromising the
    *  reference wavefunction's data.
    **/
    void shallow_copy(SharedWavefunction other);
    void shallow_copy(const Wavefunction* other);

    /**
    * Copy the contents of another Wavefunction into this one.
    * Useful at the beginning of correlated wavefunction computations.
    * -Does not set options or callbacks
    * -reference_wavefunction_ is set to other
    * -Matrices and Vectors (Ca,Da,Fa,epsilon_a, etc) are deep copied.
    **/
    void deep_copy(SharedWavefunction other);
    void deep_copy(const Wavefunction* other);

    virtual ~Wavefunction();

    /// Compute energy. Subclasses override this function to compute its energy.
    virtual double compute_energy() {throw PSIEXCEPTION("Compute energy has not been implemented for this wavefunction.");}

    /// Compute gradient.  Subclasses override this function to compute the gradient.
    virtual SharedMatrix compute_gradient() {throw PSIEXCEPTION("Analytic gradients are not available for this wavefunction.");}

    /// Compute Hessian.  Subclasses override this function to compute the Hessian.
    virtual SharedMatrix compute_hessian() {throw PSIEXCEPTION("Analytic Hessians are not available for this wavefunction.");}

    /// Is this a restricted wavefunction?
    bool same_a_b_orbs() const { return same_a_b_orbs_; }
    bool same_a_b_dens() const { return same_a_b_dens_; }

    /// Takes a Dimension object (e.g. DOCC) and returns a new Dimension object
    /// with occupations mapped to the current point group
    Dimension map_irreps(const Dimension &dimpi);

    /// Returns the molecule object that pertains to this wavefunction.
    std::shared_ptr<Molecule> molecule() const;
    std::shared_ptr<PSIO> psio() const;
    Options& options() const;

    /// An integral factory with basisset() on each center.
    std::shared_ptr<IntegralFactory> integral() const;
    /// Returns the basis set object that pertains to this wavefunction.
    std::shared_ptr<BasisSet> basisset() const;
    /// Returns this wavefunction's ECP basisset
    std::shared_ptr<BasisSet> ecpbasisset() const;
    /// Returns the SO basis set object that pertains to this wavefunction.
    std::shared_ptr<SOBasisSet> sobasisset() const;

    /// Getters and setters for other basis sets
    std::shared_ptr<BasisSet> get_basisset(std::string label);
    void set_basisset(std::string label, std::shared_ptr<BasisSet> basis);
    bool basisset_exists(std::string label);


    /// Returns the MatrixFactory object that pertains to this wavefunction
    std::shared_ptr<MatrixFactory> matrix_factory() const;
    /// Returns the reference wavefunction
    std::shared_ptr<Wavefunction> reference_wavefunction() const;
    /// Sets the reference wavefunction
    void set_reference_wavefunction(const std::shared_ptr<Wavefunction> wfn);

    /// Returns whether this wavefunction was obtained using density fitting or not
    bool density_fitted() const { return density_fitted_; }

    static void initialize_singletons();

    /// Returns the DOCC per irrep array.
    const Dimension& doccpi() const { return doccpi_; }
    /// Returns the SOCC per irrep array.
    const Dimension& soccpi() const { return soccpi_; }
    /// Returns the number of SOs per irrep array.
    const Dimension& nsopi() const { return nsopi_; }
    /// Returns the number of MOs per irrep array.
    const Dimension& nmopi() const { return nmopi_; }
    /// Returns the number of alpha electrons per irrep array.
    const Dimension& nalphapi() const { return nalphapi_; }
    /// Returns the number of beta electrons per irrep array.
    const Dimension& nbetapi() const { return nbetapi_; }
    /// Returns the frozen core orbitals per irrep array.
    const Dimension& frzcpi() const { return frzcpi_; }
    /// Returns the frozen virtual orbitals per irrep array.
    const Dimension& frzvpi() const { return frzvpi_; }

    void set_doccpi(const Dimension& doccpi);
    void set_soccpi(const Dimension& soccpi);

    /// Sets the frozen virtual orbitals per irrep array.
    void set_frzvpi(const Dimension& frzvpi) { for(int h=0; h < nirrep_; h++) frzvpi_[h] = frzvpi[h]; }

    /// Return the number of frozen core orbitals
    int nfrzc() const { return nfrzc_; }
    /// Return the number of alpha electrons
    int nalpha() const { return nalpha_; }
    /// Return the number of beta electrons
    int nbeta() const { return nbeta_; }
    /// Returns the number of SOs
    int nso() const { return nso_; }
    /// Returns the number of MOs
    int nmo() const { return nmo_; }
    /// Returns the number of irreps
    int nirrep() const { return nirrep_; }
    /// Returns the reference energy
    double reference_energy () const { return energy_; }
    /// Returns the frozen-core energy
    double efzc() const { return efzc_; }
    /// Sets the frozen-core energy
    void set_efzc(double efzc) { efzc_ = efzc; }

    /// Returns the overlap matrix
    SharedMatrix S() const { return S_; }

    /// Returns the core Hamiltonian matrix
    SharedMatrix H() const { return H_; }

    /// Returns the alpha electrons MO coefficients
    SharedMatrix Ca() const;
    /// Returns the beta electrons MO coefficients
    SharedMatrix Cb() const;
    /// Returns the (SO basis) alpha Fock matrix
    SharedMatrix Fa() const;
    /// Returns the (SO basis) beta Fock matrix
    SharedMatrix Fb() const;
    /// Returns the alpha orbital energies
    std::shared_ptr<Vector> epsilon_a() const;
    /// Returns the beta orbital energies
    std::shared_ptr<Vector> epsilon_b() const;
    /// Returns the SO basis Lagrangian
    std::shared_ptr<Matrix> Lagrangian() const;
    /// The two particle density matrix contribution to the gradient
    virtual std::shared_ptr<Matrix> tpdm_gradient_contribution() const;

    SharedMatrix aotoso() const { return AO2SO_; }

    std::shared_ptr<OEProp> get_oeprop() const { return oeprop_; }
    void set_oeprop( std::shared_ptr<OEProp> oeprop ) { oeprop_ = oeprop; }

    /// Returns the alpha OPDM for the wavefunction
    const SharedMatrix Da() const;
    /// Returns the beta OPDM for the wavefunction
    SharedMatrix Db() const;

    /**
    * Return a subset of the Ca matrix in a desired basis
    * @param basis the symmetry basis to use
    *  AO, SO
    * @param subset the subset of orbitals to return
    *  ALL, ACTIVE, FROZEN, OCC, VIR, FROZEN_OCC, ACTIVE_OCC, ACTIVE_VIR, FROZEN_VIR
    * @return the matrix in Pitzer order in the desired basis
    **/
    SharedMatrix Ca_subset(const std::string& basis = "SO", const std::string& subset = "ALL");

    /**
    * Return a subset of the Cb matrix in a desired basis
    * @param basis the symmetry basis to use
    *  AO, SO
    * @param subset the subset of orbitals to return
    *  ALL, ACTIVE, FROZEN, OCC, VIR, FROZEN_OCC, ACTIVE_OCC, ACTIVE_VIR, FROZEN_VIR
    * @return the matrix in Pitzer order in the desired basis
    **/
    SharedMatrix Cb_subset(const std::string& basis = "SO", const std::string& subset = "ALL");

    /**
     * @brief Creates an OrbitalSpace object containing information about the request alpha orbital space.
     * @param id unique name for the orbital space
     * @param basis the symmetry basis to use
     *  AO, SO
     * @param subset the subset of orbitals to return
     *  ALL, ACTIVE, FROZEN, OCC, VIR, FROZEN_OCC, ACTIVE_OCC, ACTIVE_VIR, FROZEN_VIR
     * @return OrbitalSpace object containing data for the requested space.
     */
    OrbitalSpace alpha_orbital_space(const std::string& id, const std::string& basis = "SO", const std::string& subset = "ALL");
    /**
     * @brief Creates an OrbitalSpace object containing information about the request beta orbital space.
     * @param id unique name for the orbital space
     * @param basis the symmetry basis to use
     *  AO, SO
     * @param subset the subset of orbitals to return
     *  ALL, ACTIVE, FROZEN, OCC, VIR, FROZEN_OCC, ACTIVE_OCC, ACTIVE_VIR, FROZEN_VIR
     * @return OrbitalSpace object containing data for the requested space.
     */
    OrbitalSpace beta_orbital_space(const std::string& id, const std::string& basis = "SO", const std::string& subset = "ALL");

    /**
    * Return the Da matrix in the desired basis
    * @param basis the symmetry basis to use
    *  AO, SO, MO
    * @return the matrix in the desired basis
    **/
    SharedMatrix Da_subset(const std::string& basis = "SO");

    /**
    * Return the Db matrix in the desired basis
    * @param basis the symmetry basis to use
    *  AO, SO, MO
    * @return the matrix in the desired basis
    **/
    SharedMatrix Db_subset(const std::string& basis = "SO");

    /**
    * Return the D matrix in the desired basis
    * @param D matrix in the SO basis to transform
    * @param C matrix in the SO basis to use as a transformer
    * @param basis the symmetry basis to use
    *  AO, SO, MO, CartAO
    * @return the D matrix in the desired basis
    **/
    SharedMatrix D_subset_helper(SharedMatrix D, SharedMatrix C, const std::string& basis);

    /**
    * Return the alpha orbital eigenvalues in the desired basis
    * @param basis the symmetry basis to use
    *  AO, SO, MO (SO and MO return the same thing)
    * @param subset the subset of orbitals to return
    *  ALL, ACTIVE, FROZEN, OCC, VIR, FROZEN_OCC, ACTIVE_OCC, ACTIVE_VIR, FROZEN_VIR
    */
    SharedVector epsilon_a_subset(const std::string& basis = "SO", const std::string& subset = "ALL");

    /**
    * Return the beta orbital eigenvalues in the desired basis
    * @param basis the symmetry basis to use
    *  AO, SO, MO (SO and MO return the same thing)
    * @param subset the subset of orbitals to return
    *  ALL, ACTIVE, FROZEN, OCC, VIR, FROZEN_OCC, ACTIVE_OCC, ACTIVE_VIR, FROZEN_VIR
    */
    SharedVector epsilon_b_subset(const std::string& basis = "SO", const std::string& subset = "ALL");

    /**
     * Projects the given orbitals from the old to the new basis
     * @param  Cold      Orbitals to project
     * @param  noccpi    Number of eigenvectors of importance
     * @param  old_basis The old basis set
     * @param  new_basis The new basis set
     * @return           The projected basis (nso x noccpi)
     */
    SharedMatrix basis_projection(SharedMatrix Cold, Dimension noccpi,
                                  std::shared_ptr<BasisSet> old_basis,
                                  std::shared_ptr<BasisSet> new_basis);

    /// Returns the Lagrangian in SO basis for the wavefunction
    SharedMatrix X() const;

    /// Returns the gradient
    SharedMatrix gradient() const;
    /// Set the gradient for the wavefunction
    void set_gradient(SharedMatrix& grad);

    /// Returns the Hessian
    SharedMatrix hessian() const;
    /// Set the Hessian for the wavefunction
    void set_hessian(SharedMatrix& hess);

    /// Returns the atomic point charges
    std::shared_ptr<std::vector<double>> atomic_point_charges()const{
       return atomic_point_charges_;
    }
    /// Returns the atomic point charges in Vector form for python output.
    std::shared_ptr<Vector> get_atomic_point_charges() const;

    /// Sets the atomic point charges
    void set_atomic_point_charges(const std::shared_ptr<std::vector<double>>& apcs){
       atomic_point_charges_=apcs;
    }

    /// Returns the frequencies
    std::shared_ptr<Vector> frequencies() const;
    /// Set the frequencies for the wavefunction
    void set_frequencies(std::shared_ptr<Vector>& freqs);

    /// Returns the normalmodes
    std::shared_ptr<Vector> normalmodes() const;
    /// Set the normalmodes for the wavefunction
    void set_normalmodes(std::shared_ptr<Vector>& norms);

    /// Set the wavefunction name (e.g. "RHF", "ROHF", "UHF", "CCEnergyWavefunction")
    void set_name(const std::string& name) { name_ = name; }

    /// Returns the wavefunction name
    const std::string& name() const { return name_; }

    // Set the print flag level
    void set_print(unsigned int print) { print_ = print; }

    // Set the debug flag level
    void set_debug(unsigned int debug) { debug_ = debug; }

    /// Save the wavefunction to checkpoint
    virtual void save() const;

    /// Get and set variables dictionary
    double get_variable(const std::string key);
    void set_variable(const std::string key, double value) { variables_[key] = value; }
    std::map<std::string, double> variables(void) { return variables_; }

    /// Get and set arrays dictionary
    SharedMatrix get_array(const std::string key);
    void set_array(const std::string key, SharedMatrix value) { arrays_[key] = value; }
    std::map<std::string, SharedMatrix> arrays(void) { return arrays_; }
};

}

#endif
