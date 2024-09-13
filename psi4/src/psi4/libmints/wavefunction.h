/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi_src_lib_libmints_wavefunction_h
#define _psi_src_lib_libmints_wavefunction_h

#include "typedefs.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/dimension.h"

#include <cstddef>
#include <vector>
#include <array>
#include <memory>
#include <map>

#define MAX_IOFF 30000
extern size_t ioff[MAX_IOFF];

#define MAX_DF 500
extern double df[MAX_DF];

#define MAX_BC 20
extern double bc[MAX_BC][MAX_BC];

#define MAX_FAC 100
extern double fac[MAX_FAC];

#if !defined(EXPLICIT_IOFF)
#define EXPLICIT_IOFF(i) ((i) * ((i) + 1) / 2)
#endif

#if !defined(INDEX2)
#define INDEX2(i, j) ((i) >= (j) ? EXPLICIT_IOFF(i) + (j) : EXPLICIT_IOFF(j) + (i))
#endif

#if !defined(INDEX4)
#define INDEX4(i, j, k, l) (INDEX2(INDEX2((i), (j)), INDEX2((k), (l))))
#endif

namespace psi {

class Molecule;
class BasisSet;
class IntegralFactory;
class MintsHelper;
class Matrix;
class Vector;
class MatrixFactory;
class Options;
class SOBasisSet;
class PCM;
class PSIO;
class OrbitalSpace;
class ExternalPotential;

/*! \ingroup MINTS
 *  \class Wavefunction
 *  \brief Simple wavefunction base class.
 */
class PSI_API Wavefunction : public std::enable_shared_from_this<Wavefunction> {
   protected:
    /// Name of the wavefunction
    std::string name_;

    /// Module name for CURRENT ENERGY
    std::string module_;

    /// The ORBITAL basis
    std::shared_ptr<BasisSet> basisset_;

    /// Primary basis set for SO integrals
    std::shared_ptr<SOBasisSet> sobasisset_;

    /// AO2SO conversion matrix (AO in rows, SO in cols)
    SharedMatrix AO2SO_;

    /// Molecule that this wavefunction is run on
    std::shared_ptr<Molecule> molecule_;

    /// Options object
    Options& options_;

    // PSI file access variables
    std::shared_ptr<PSIO> psio_;

    /// Integral factory
    std::shared_ptr<IntegralFactory> integral_;

    /// MintsHelper
    std::shared_ptr<MintsHelper> mintshelper_;

    /// Matrix factory for creating standard sized matrices
    std::shared_ptr<MatrixFactory> factory_;

    std::shared_ptr<Wavefunction> reference_wavefunction_;

    /// How much memory you have access to.
    long int memory_;

    /// Perturb the Hamiltonian?
    int perturb_h_;
    /// With what...
    enum FieldType { nothing, dipole_x, dipole_y, dipole_z, dipole, embpot, dx, sphere };
    FieldType dipole_field_type_;
    /// How big of a field perturbation to apply
    std::array<double, 3> dipole_field_strength_;

    /// Debug flag
    size_t debug_;
    /// Print flag
    size_t print_;

    /// Total alpha and beta electrons
    int nalpha_, nbeta_;

    /// Total frozen core orbitals
    int nfrzc_;

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
    SharedVector epsilon_a_;
    /// Beta orbital energies
    SharedVector epsilon_b_;

    /// gradient, if available, as natom_ x 3 SharedMatrix
    SharedMatrix gradient_;

    /// Hessian, if available, as natom_*3 x natom_*3 SharedMatrix (NOT mass-weighted!)
    SharedMatrix hessian_;

    /// Helpers for C/D/epsilon transformers
    SharedMatrix C_subset_helper(SharedMatrix C, const Dimension& noccpi, SharedVector epsilon,
                                 const std::string& basis, const std::string& subset) const;
    // Return the desired subset of orbital energies.
    // @param epsilon The vector of orbital energies
    // @param noccpi The dimension of "occupied" orbitals for the case of interest.
    //               Usual use case: nalphapi_ or nbetapi_?
    // @param basis "AO", "SO", or "MO" - should the return vector use symmetry
    // @param subset A label appended to the return vector name, "Epsilon {basis} {subset}"
    SharedVector epsilon_subset_helper(SharedVector epsilon, const Dimension& noccpi, const std::string& basis,
                                       const std::string& subset) const;
    // Helper function needed by the helper functions.
    // Return the desired MO indices, per irrep.
    // @param noccpi The dimension of occupied indices used as the source-of-
    //               truth for the dimension of occupied orbitals.
    // @param subset "FROZEN_OCC", "FROZEN_VIR", "ACTIVE_OCC", "ACTIVE_VIR"
    //               "ACTIVE", "OCC", "VIR", or "ALL". The space, the indices
    //                of which are to be returned.
    std::vector<std::vector<int>> subset_occupation(const Dimension& noccpi, const std::string& subset) const;

    /// Should nuclear electrostatic potentials be available, they will be here
    std::shared_ptr<std::vector<double>> esp_at_nuclei_;

    /// Should molecular orbital extents be available, they will be here
    std::vector<SharedVector> mo_extents_;

    /// If atomic point charges are available they will be here
    std::shared_ptr<std::vector<double>> atomic_point_charges_;

    /// Should natural orbital occupations be available, they will be here
    std::vector<std::vector<std::tuple<double, int, int>>> no_occupations_;

    /// Same orbs or dens
    bool same_a_b_dens_;
    bool same_a_b_orbs_;

    // The external potential for the current wave function
    std::shared_ptr<ExternalPotential> external_pot_;

    // Collection of scalar variables
    std::map<std::string, double> variables_;

    // Collection of Matrix variables
    // * any '<mtd> GRADIENT' is an energy derivative w.r.t. nuclear perturbations (a.u.) as a (nat, 3) Matrix
    // * any '<mtd> DIPOLE GRADIENT' is a dipole derivative w.r.t. nuclear perturbations (a.u.) as a degree-of-freedom
    //   by dipole component (3 * nat, 3) Matrix
    std::map<std::string, SharedMatrix> arrays_;

    // Collection of external potentials; this member variable is provisional and might be removed in the future.
    // This member variable is currently used for passing ExternalPotential objects to the F/I-SAPT code
    // The above defined external_pot_ member variable contains the total external potential defined for the current
    // wave function. For F/I-SAPT, we need a set of external potential that can be assigned to either the interacting
    // fragments or to the environment
    // For F/I-SAPT the keys can be A, B, or C (all optionals), where A and B signify the interacting subsystem
    // and C signify the envirnoment
    std::map<std::string, std::shared_ptr<ExternalPotential>> potentials_;

    // Polarizable continuum model
    bool PCM_enabled_;
    std::shared_ptr<PCM> PCM_;

   private:
    // Wavefunction() {}
    void common_init();

   public:
    /// Constructor for an entirely new wavefunction with an existing basis
    Wavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis, Options& options);

    /// Constructor for an entirely new wavefunction with an existing basis and global options
    Wavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis);

    /// Constructor for a wavefunction deserialized from a file and initialized in the form of maps to all member
    /// variables
    Wavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basisset,
                 std::map<std::string, std::shared_ptr<Matrix>> matrices,
                 std::map<std::string, std::shared_ptr<Vector>> vectors, std::map<std::string, Dimension> dimensions,
                 std::map<std::string, int> ints, std::map<std::string, std::string> strings,
                 std::map<std::string, bool> booleans, std::map<std::string, double> floats);

    /// Blank constructor for derived classes
    Wavefunction(SharedWavefunction reference_wavefunction, Options& options);

    /// Blank constructor for derived classes
    Wavefunction(Options& options);

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

    /**
     * Creates a new wavefunction in C1-symmetry format from the current
     * Wavefunction that may be in a higher point group symmetry format.
     *
     * @param basis A C1-symmetry basis set object (we don't yet have
     *        the ability to copy this straight from the symmetric Wavefunction)
     **/
    std::shared_ptr<Wavefunction> c1_deep_copy(std::shared_ptr<BasisSet> basis);

    virtual ~Wavefunction();

    /// Compute energy. Subclasses override this function to compute its energy.
    virtual double compute_energy() {
        throw PSIEXCEPTION("Compute energy has not been implemented for this wavefunction.");
    }

    /// Compute gradient.  Subclasses override this function to compute the gradient.
    virtual SharedMatrix compute_gradient() {
        throw PSIEXCEPTION("Analytic gradients are not available for this wavefunction.");
    }

    /// Compute Hessian.  Subclasses override this function to compute the Hessian.
    virtual SharedMatrix compute_hessian() {
        throw PSIEXCEPTION("Analytic Hessians are not available for this wavefunction.");
    }

    /// Is this a restricted wavefunction?
    bool same_a_b_orbs() const { return same_a_b_orbs_; }
    bool same_a_b_dens() const { return same_a_b_dens_; }

    /// Takes a Dimension object (e.g. DOCC) and returns a new Dimension object
    /// with occupations mapped to the current point group
    Dimension map_irreps(const Dimension& dimpi);

    /// Returns the molecule object that pertains to this wavefunction.
    std::shared_ptr<Molecule> molecule() const;
    std::shared_ptr<PSIO> psio() const;
    Options& options() const;

    /// An integral factory with basisset() on each center.
    std::shared_ptr<IntegralFactory> integral() const;
    /// An molecular integrals helper with basisset() on each center.
    std::shared_ptr<MintsHelper> mintshelper() const;
    /// Returns the basis set object that pertains to this wavefunction.
    std::shared_ptr<BasisSet> basisset() const;
    /// Returns the SO basis set object that pertains to this wavefunction.
    std::shared_ptr<SOBasisSet> sobasisset() const;

    /// Getters and setters for other basis sets
    std::map<std::string, std::shared_ptr<BasisSet>> basissets() const;
    std::shared_ptr<BasisSet> get_basisset(std::string label);
    void set_basisset(std::string label, std::shared_ptr<BasisSet> basis);
    bool basisset_exists(std::string label);

    /// Returns the MatrixFactory object that pertains to this wavefunction
    std::shared_ptr<MatrixFactory> matrix_factory() const;
    /// Returns the reference wavefunction
    std::shared_ptr<Wavefunction> reference_wavefunction() const;
    /// Sets the reference wavefunction
    void set_reference_wavefunction(const std::shared_ptr<Wavefunction> wfn);

    /// Returns the print level
    int get_print() const { return print_; }
    static void initialize_singletons();

    /// Returns the DOCC per irrep array. Not recommended for unrestricted code.
    /// Flag `warn_on_beta_socc` triggers warning on singly occupied beta orbitals,
    /// which break assumptions made in pre-1.7 DOCC/SOCC handling.
    const Dimension doccpi(bool warn_on_beta_socc = true) const;
    /// Returns the SOCC per irrep array. Not recommended for unrestricted code.
    /// Flag `warn_on_beta_socc` triggers warning on singly occupied beta orbitals,
    /// which break assumptions made in pre-1.7 DOCC/SOCC handling.
    const Dimension soccpi(bool warn_on_beta_socc = true) const;
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

    /* Return the magnitude of the dipole perturbation strength in the x,y,z direction */
    std::array<double, 3> get_dipole_field_strength() const;
    FieldType get_dipole_perturbation_type() const;

    /**
     * @brief Expert specialized use only. Sets the number of doubly and singly occupied orbitals per irrep. Results in an
     * inconsistent Wavefunction object for SCF purposes, so caution is advised.
     * @param doccpi the new list of doubly occupied orbitals per irrep
     */
    void force_occpi(const Dimension& input_doccpi, const Dimension& input_soccpi);

    /// Sets the frozen virtual orbitals per irrep array.
    void set_frzvpi(const Dimension& frzvpi);

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
    /// Returns the energy
    PSI_DEPRECATED(
        "Using `Wavefunction.reference_energy` instead of `Wavefunction.energy` is deprecated, and as soon as 1.4 it will "
        "stop working")
    double reference_energy() const { return energy_; }
    double energy() const { return energy_; }
    /// Sets the energy
    void set_energy(double ene);
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
    SharedVector epsilon_a() const;
    /// Returns the beta orbital energies
    SharedVector epsilon_b() const;

    SharedMatrix aotoso() const { return AO2SO_; }

    /// Returns the alpha OPDM for the wavefunction
    const SharedMatrix Da() const;
    /// Returns the beta OPDM for the wavefunction
    const SharedMatrix Db() const;

    /**
     * Return a subset of the Ca matrix in a desired basis
     * @param basis the symmetry basis to use
     *  AO, SO
     * @param subset the subset of orbitals to return
     *  ALL, ACTIVE, FROZEN, OCC, VIR, FROZEN_OCC, ACTIVE_OCC, ACTIVE_VIR, FROZEN_VIR
     * @return the matrix in Pitzer order in the desired basis
     *  Pitzer ordering is in c1 symmetry if AO is selected
     **/
    SharedMatrix Ca_subset(const std::string& basis = "SO", const std::string& subset = "ALL") const;

    /**
     * Return a subset of the Cb matrix in a desired basis
     * @param basis the symmetry basis to use
     *  AO, SO
     * @param subset the subset of orbitals to return
     *  ALL, ACTIVE, FROZEN, OCC, VIR, FROZEN_OCC, ACTIVE_OCC, ACTIVE_VIR, FROZEN_VIR
     * @return the matrix in Pitzer order in the desired basis
     *  Pitzer ordering is in c1 symmetry if AO is selected
     **/
    SharedMatrix Cb_subset(const std::string& basis = "SO", const std::string& subset = "ALL") const;

    /**
     * @brief Creates an OrbitalSpace object containing information about the request alpha orbital space.
     * @param id unique name for the orbital space
     * @param basis the symmetry basis to use
     *  AO, SO
     * @param subset the subset of orbitals to return
     *  ALL, ACTIVE, FROZEN, OCC, VIR, FROZEN_OCC, ACTIVE_OCC, ACTIVE_VIR, FROZEN_VIR
     * @return OrbitalSpace object containing data for the requested space.
     */
    OrbitalSpace alpha_orbital_space(const std::string& id, const std::string& basis = "SO",
                                     const std::string& subset = "ALL");
    /**
     * @brief Creates an OrbitalSpace object containing information about the request beta orbital space.
     * @param id unique name for the orbital space
     * @param basis the symmetry basis to use
     *  AO, SO
     * @param subset the subset of orbitals to return
     *  ALL, ACTIVE, FROZEN, OCC, VIR, FROZEN_OCC, ACTIVE_OCC, ACTIVE_VIR, FROZEN_VIR
     * @return OrbitalSpace object containing data for the requested space.
     */
    OrbitalSpace beta_orbital_space(const std::string& id, const std::string& basis = "SO",
                                    const std::string& subset = "ALL");

    /**
     * Return the Da matrix in the desired basis
     * @param basis the symmetry basis to use
     *  AO, SO, MO
     * @return the matrix in the desired basis
     **/
    SharedMatrix Da_subset(const std::string& basis = "SO") const;

    /**
     * Return the Db matrix in the desired basis
     * @param basis the symmetry basis to use
     *  AO, SO, MO
     * @return the matrix in the desired basis
     **/
    SharedMatrix Db_subset(const std::string& basis = "SO") const;

    /**
     * Return the Fa matrix in the desired basis
     * @param basis the symmetry basis to use
     *  AO, SO, MO
     * @return the matrix in the desired basis
     **/
    SharedMatrix Fa_subset(const std::string& basis = "SO") const;

    /**
     * Return the Fb matrix in the desired basis
     * @param basis the symmetry basis to use
     *  AO, SO, MO
     * @return the matrix in the desired basis
     **/
    SharedMatrix Fb_subset(const std::string& basis = "SO") const;

    /**
     * Transform a matrix M into the desired basis
     * @param M matrix in the SO basis to transform
     * @param C matrix in the SO basis to use for transforms to MO basis
     * @param basis the symmetry basis to use
     *  AO, SO, MO, CartAO
     * @param MO_as_overlap whether the AO-to-MO transformation requires the overlap.
     *   Only used when MO basis requested.
     * @return the matrix M in the desired basis
     **/
    SharedMatrix matrix_subset_helper(SharedMatrix M, SharedMatrix C, const std::string& basis,
                                      const std::string matrix_basename, bool MO_as_overlap = true) const;

    /**
     * Return the alpha orbital eigenvalues in the desired basis
     * @param basis the symmetry basis to use
     *  AO, SO, MO (SO and MO return the same thing)
     * @param subset the subset of orbitals to return
     *  ALL, ACTIVE, FROZEN, OCC, VIR, FROZEN_OCC, ACTIVE_OCC, ACTIVE_VIR, FROZEN_VIR
     */
    SharedVector epsilon_a_subset(const std::string& basis = "SO", const std::string& subset = "ALL") const;

    /**
     * Return the beta orbital eigenvalues in the desired basis
     * @param basis the symmetry basis to use
     *  AO, SO, MO (SO and MO return the same thing)
     * @param subset the subset of orbitals to return
     *  ALL, ACTIVE, FROZEN, OCC, VIR, FROZEN_OCC, ACTIVE_OCC, ACTIVE_VIR, FROZEN_VIR
     */
    SharedVector epsilon_b_subset(const std::string& basis = "SO", const std::string& subset = "ALL") const;

    /**
     * Projects the given orbitals from the old to the new basis
     * @param  Cold      Orbitals to project
     * @param  noccpi    Number of eigenvectors of importance
     * @param  old_basis The old basis set
     * @param  new_basis The new basis set
     * @return           The projected basis (nso x noccpi)
     */
    SharedMatrix basis_projection(SharedMatrix Cold, Dimension noccpi, std::shared_ptr<BasisSet> old_basis,
                                  std::shared_ptr<BasisSet> new_basis);

    /// Returns the SO basis Lagrangian
    SharedMatrix lagrangian() const;
    /// Set Lagrangian matrix in SO basis
    void set_lagrangian(SharedMatrix X);
    /// Returns the SO basis Lagrangian
    PSI_DEPRECATED(
        "Using `Wavefunction.Lagrangian` instead of `Wavefunction.lagrangian` is deprecated,"
        " and as soon as 1.5 it will stop working")
    SharedMatrix Lagrangian() const { return lagrangian(); }
    /// Returns the SO basis Lagrangian (duplicated one)
    PSI_DEPRECATED(
        "Using `Wavefunction.X` instead of `Wavefunction.lagrangian` is deprecated,"
        " and as soon as 1.5 it will stop working")
    SharedMatrix X() const { return lagrangian(); }

    /// Returns the gradient
    SharedMatrix gradient() const;
    /// Set the gradient for the wavefunction
    void set_gradient(SharedMatrix grad);

    /// Returns the Hessian
    SharedMatrix hessian() const;
    /// Set the Hessian for the wavefunction
    void set_hessian(SharedMatrix hess);

    /// Returns electrostatic potentials at nuclei
    std::shared_ptr<std::vector<double>> esp_at_nuclei() const { return esp_at_nuclei_; }

    /// Returns electrostatic potentials at nuclei in Vector form for python output
    std::shared_ptr<Vector> get_esp_at_nuclei() const;

    /// Sets the electrostatic potentials at nuclei
    void set_esp_at_nuclei(const std::shared_ptr<std::vector<double>>& nesps) { esp_at_nuclei_ = nesps; }

    /// Returns Molecular orbital extents
    std::vector<SharedVector> mo_extents() const { return mo_extents_; }

    /// Returns Molecular orbital extents in Vector form for python output.
    std::vector<SharedVector> get_mo_extents() const;

    /// Sets molecular orbital extents
    void set_mo_extents(const std::vector<SharedVector> mo_es) { mo_extents_ = mo_es; }

    /// Returns the atomic point charges
    std::shared_ptr<std::vector<double>> atomic_point_charges() const { return atomic_point_charges_; }
    /// Returns the atomic point charges in Vector form for python output.
    SharedVector get_atomic_point_charges() const;

    /// Sets the atomic point charges
    void set_atomic_point_charges(const std::shared_ptr<std::vector<double>>& apcs) { atomic_point_charges_ = apcs; }

    /// Returns NO occupations
    std::vector<std::vector<std::tuple<double, int, int>>> no_occupations() const { return no_occupations_; }

    /// Returns the NO occupations in vector form for python output
    std::vector<std::vector<std::tuple<double, int, int>>> get_no_occupations() const;

    /// Sets the NO occupations
    void set_no_occupations(const std::vector<std::vector<std::tuple<double, int, int>>> no_ocs) {
        no_occupations_ = no_ocs;
    }

    /// Set the wavefunction name (e.g. "RHF", "ROHF", "UHF", "CCEnergyWavefunction")
    void set_name(const std::string& name) { name_ = name; }

    /// Returns the wavefunction name
    const std::string& name() const { return name_; }

    /// Set the module name (e.g. "OCC", "CCENERGY", "CCT3")
    void set_module(const std::string& module) { module_ = module; }

    /// Returns the module name
    const std::string& module() const { return module_; }

    // Set the print flag level
    void set_print(size_t print) { print_ = print; }

    // Set the debug flag level
    void set_debug(size_t debug) { debug_ = debug; }

    /// Save the wavefunction to checkpoint
    virtual void save() const;

    // Get the external potential
    std::shared_ptr<ExternalPotential> external_pot() const;

    // Set the external potential
    void set_external_potential(std::shared_ptr<ExternalPotential> external) { external_pot_ = external; }

    /// Get and set variables, arrays, and potentials dictionaries
    bool has_scalar_variable(const std::string& key);
    bool has_array_variable(const std::string& key);
    // The function below is provisional and might be removed in the future
    bool has_potential_variable(const std::string& key);
    double scalar_variable(const std::string& key);
    SharedMatrix array_variable(const std::string& key);
    // The function below is provisional and might be removed in the future
    std::shared_ptr<ExternalPotential> potential_variable(const std::string& key);
    void set_scalar_variable(const std::string& key, double value);
    void set_array_variable(const std::string& key, SharedMatrix value);
    // The function below is provisional and might be removed in the future
    void set_potential_variable(const std::string& key, std::shared_ptr<ExternalPotential> value);
    int del_scalar_variable(const std::string& key);
    int del_array_variable(const std::string& key);
    // The function below is provisional and might be removed in the future
    int del_potential_variable(const std::string& key);
    std::map<std::string, double> scalar_variables();
    std::map<std::string, SharedMatrix> array_variables();
    // The function below is provisional and might be removed in the future
    std::map<std::string, std::shared_ptr<ExternalPotential>> potential_variables();

    PSI_DEPRECATED(
        "Using `Wavefunction.get_variable` instead of `Wavefunction.scalar_variable` is deprecated, and as soon as 1.4 it will "
        "stop working")
    double get_variable(const std::string& key);
    PSI_DEPRECATED(
        "Using `Wavefunction.set_variable` instead of `Wavefunction.set_scalar_variable` is deprecated, and as soon as 1.4 it "
        "will stop working")
    void set_variable(const std::string& key, double value);
    PSI_DEPRECATED(
        "Using `Wavefunction.variables` instead of `Wavefunction.scalar_variables` is deprecated, and as soon as 1.4 it will "
        "stop working")
    std::map<std::string, double> variables();
    PSI_DEPRECATED(
        "Using `Wavefunction.get_array` instead of `Wavefunction.array_variable` is deprecated, and as soon as 1.4 it will "
        "stop working")
    SharedMatrix get_array(const std::string& key);
    PSI_DEPRECATED(
        "Using `Wavefunction.set_array` instead of `Wavefunction.set_array_variable` is deprecated, and as soon as 1.4 it will "
        "stop working")
    void set_array(const std::string& key, SharedMatrix value);
    PSI_DEPRECATED(
        "Using `Wavefunction.arrays` instead of `Wavefunction.array_variables` is deprecated, and as soon as 1.4 it will stop "
        "working")
    std::map<std::string, SharedMatrix> arrays();

    /// Set PCM object
    void set_PCM(const std::shared_ptr<PCM>& pcm);
    /// Get PCM object
    std::shared_ptr<PCM> get_PCM() const;
    bool PCM_enabled() const { return PCM_enabled_; }

    /// The below members are experimental and are designed to hold densities when the
    /// "current density" is ambiguous, e.g., non-orbital optimized methods and multi-
    /// stage methods. ~ JPM - Apr. '22
    /// This is public because `ccdensity` doesn't subclass wfn like it should, so we need
    /// SOME way to let it get/set.
    /// Vector of density matrices
    std::map<std::string, SharedMatrix> density_map_;

};

}  // namespace psi

#endif
