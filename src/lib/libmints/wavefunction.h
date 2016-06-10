/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#include <stddef.h>
#include <vector>

#include "typedefs.h"
#include "exception.h"
#include <libparallel/parallel.h>

#include <boost/enable_shared_from_this.hpp>

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

#if !defined( PyObject_HEAD )
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
class OrbitalSpace;

/*! \ingroup MINTS
 *  \class Wavefunction
 *  \brief Simple wavefunction base class.
 */
class Wavefunction : public boost::enable_shared_from_this<Wavefunction>
{
protected:
    /// Name of the wavefunction
    std::string name_;

    /// Primary basis set for AO integrals
    boost::shared_ptr<BasisSet> basisset_;

    /// Primary basis set for SO integrals
    boost::shared_ptr<SOBasisSet> sobasisset_;

    /// AO2SO conversion matrix (AO in rows, SO in cols)
    SharedMatrix AO2SO_;

    /// Molecule that this wavefunction is run on
    boost::shared_ptr<Molecule> molecule_;

    /// Options object
    Options & options_;

    // PSI file access variables
    boost::shared_ptr<PSIO> psio_;

    /// Integral factory
    boost::shared_ptr<IntegralFactory> integral_;

    /// Matrix factory for creating standard sized matrices
    boost::shared_ptr<MatrixFactory> factory_;

    boost::shared_ptr<Wavefunction> reference_wavefunction_;

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
    boost::shared_ptr<Vector> epsilon_a_;
    /// Beta orbital energies
    boost::shared_ptr<Vector> epsilon_b_;

    // Callback routines to Python
    std::vector<void*> precallbacks_;
    std::vector<void*> postcallbacks_;

    /// If a gradient is available it will be here:
    SharedMatrix gradient_;

    /// If a Hessian is available it will be here:
    SharedMatrix hessian_;

    /// The TPDM contribution to the gradient
    boost::shared_ptr<Matrix> tpdm_gradient_contribution_;

    /// Helpers for C/D/epsilon transformers
    SharedMatrix C_subset_helper(SharedMatrix C, const Dimension& noccpi, SharedVector epsilon, const std::string& basis, const std::string& subset);
    SharedMatrix F_subset_helper(SharedMatrix F, SharedMatrix C, const std::string& basis);
    SharedVector epsilon_subset_helper(SharedVector epsilon, const Dimension& noccpi, const std::string& basis, const std::string& subset);
    std::vector<std::vector<int> > subset_occupation(const Dimension& noccpi, const std::string& subset);

    /// If atomic point charges are available they will be here
    boost::shared_ptr<double[]> atomic_point_charges_;

    /// If frequencies are available, they will be here:
    boost::shared_ptr<Vector> frequencies_;

    /// If normal modes are available, they will be here:
    boost::shared_ptr<Vector> normalmodes_;

    /// Same orbs or dens
    bool same_a_b_dens_;
    bool same_a_b_orbs_;

private:
    // Wavefunction() {}
    void common_init();

public:

    /// Constructor for an entirely new wavefunction
    Wavefunction(boost::shared_ptr<Molecule> molecule, const std::string& basis,
                 Options & options);

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

    /// Takes an irrep-by-irrep array (e.g. DOCC) and maps it into the current point group
    void map_irreps(std::vector<int*> &arrays);
    /// A wrapper to the Vector version of this function
    void map_irreps(int* &array);
    /// A wrapper to the Vector version of this function
    void map_irreps(Dimension &array);

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
    /// Sets the frozen virtual orbitals per irrep array.
    void set_frzvpi(const Dimension& frzvpi) { for(int h=0; h < nirrep_; h++) frzvpi_[h] = frzvpi[0]; }
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
    boost::shared_ptr<Vector> epsilon_a() const;
    /// Returns the beta orbital energies
    boost::shared_ptr<Vector> epsilon_b() const;
    /// Returns the SO basis Lagrangian
    boost::shared_ptr<Matrix> Lagrangian() const;
    /// The two particle density matrix contribution to the gradient
    virtual boost::shared_ptr<Matrix> tpdm_gradient_contribution() const;

    SharedMatrix aotoso() const { return AO2SO_; }

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

    /// Returns the Hessian
    SharedMatrix hessian() const;
    /// Set the Hessian for the wavefunction
    void set_hessian(SharedMatrix& hess);

    /// Returns the atomic point charges
    boost::shared_ptr<double[]> atomic_point_charges()const{
       return atomic_point_charges_;
    }
    /// Returns the atomic point charges in Vector form for python output.
    boost::shared_ptr<Vector> get_atomic_point_charges() const;

    /// Sets the atomic point charges
    void set_atomic_point_charges(const boost::shared_ptr<double[]>& apcs){
       atomic_point_charges_=apcs;
    }

    /// Returns the frequencies
    boost::shared_ptr<Vector> frequencies() const;
    /// Set the frequencies for the wavefunction
    void set_frequencies(boost::shared_ptr<Vector>& freqs);

    /// Returns the normalmodes
    boost::shared_ptr<Vector> normalmodes() const;
    /// Set the normalmodes for the wavefunction
    void set_normalmodes(boost::shared_ptr<Vector>& norms);

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
};

}

#endif
