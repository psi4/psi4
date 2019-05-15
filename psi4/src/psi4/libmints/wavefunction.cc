/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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
#include "psi4/libmints/wavefunction.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/psifiles.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/orbitalspace.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/vector3.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/corrtab.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/process.h"
#ifdef USING_PCMSolver
#include "psi4/libpsipcm/psipcm.h"
#endif

#include <typeinfo>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <tuple>

using namespace psi;

// Globals (Seriously? This is where we instantiate these.)
size_t ioff[MAX_IOFF];
double df[MAX_DF];
double bc[MAX_BC][MAX_BC];
double fac[MAX_FAC];

Wavefunction::Wavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis, Options &options)
    : options_(options),
      basisset_(basis),
      molecule_(molecule),
      dipole_field_strength_{{0.0, 0.0, 0.0}},
      PCM_enabled_(false) {
    common_init();
}

Wavefunction::Wavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis)
    : options_(Process::environment.options),
      basisset_(basis),
      molecule_(molecule),
      dipole_field_strength_{{0.0, 0.0, 0.0}},
      PCM_enabled_(false) {
    common_init();
}

Wavefunction::Wavefunction(SharedWavefunction reference_wavefunction, Options &options)
    : options_(options), dipole_field_strength_{{0.0, 0.0, 0.0}}, PCM_enabled_(false) {
    // Copy the wavefuntion then update
    shallow_copy(reference_wavefunction);
    set_reference_wavefunction(reference_wavefunction);
}

// TODO: pass Options object to constructor instead of relying on globals
Wavefunction::Wavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basisset,
                           std::map<std::string, std::shared_ptr<Matrix>> matrices,
                           std::map<std::string, std::shared_ptr<Vector>> vectors,
                           std::map<std::string, Dimension> dimensions, std::map<std::string, int> ints,
                           std::map<std::string, std::string> strings, std::map<std::string, bool> booleans,
                           std::map<std::string, double> floats)
    : options_(Process::environment.options), basisset_(basisset), molecule_(molecule) {
    // Check the point group of the molecule. If it is not set, set it.
    if (!molecule_->point_group()) {
        molecule_->set_point_group(molecule_->find_point_group());
    }

    // Create an SO basis...we need the point group for this part.
    integral_ = std::make_shared<IntegralFactory>(basisset_, basisset_, basisset_, basisset_);
    sobasisset_ = std::make_shared<SOBasisSet>(basisset_, integral_);

    // set matrices
    Ca_ = matrices["Ca"];
    Cb_ = matrices["Cb"];
    Da_ = matrices["Da"];
    Db_ = matrices["Db"];
    Fa_ = matrices["Fa"];
    Fb_ = matrices["Fb"];
    H_ = matrices["H"];
    S_ = matrices["S"];
    Lagrangian_ = matrices["X"];
    AO2SO_ = matrices["aotoso"];
    gradient_ = matrices["gradient"];
    hessian_ = matrices["hessian"];

    // set vectors
    epsilon_a_ = vectors["epsilon_a"];
    epsilon_b_ = vectors["epsilon_b"];
    frequencies_ = vectors["frequencies"];

    // set dimensions
    doccpi_ = dimensions["doccpi"];
    soccpi_ = dimensions["soccpi"];
    frzcpi_ = dimensions["frzcpi"];
    frzvpi_ = dimensions["frzvpi"];
    nalphapi_ = dimensions["nalphapi"];
    nbetapi_ = dimensions["nbetapi"];
    nmopi_ = dimensions["nmopi"];
    nsopi_ = dimensions["nsopi"];

    // set integers
    nalpha_ = ints["nalpha"];
    nbeta_ = ints["nbeta"];
    nfrzc_ = ints["nfrzc"];
    nirrep_ = ints["nirrep"];
    nmo_ = ints["nmo"];
    nso_ = ints["nso"];
    print_ = ints["print"];

    // set strings
    name_ = strings["name"];

    // set booleans
    PCM_enabled_ = booleans["PCM_enabled"];
    same_a_b_dens_ = booleans["same_a_b_dens"];
    same_a_b_orbs_ = booleans["same_a_b_orbs"];
    density_fitted_ = booleans["density_fitted"];

    // set floats
    energy_ = floats["energy"];
    efzc_ = floats["efzc"];
    dipole_field_strength_[0] = floats["dipole_field_x"];
    dipole_field_strength_[1] = floats["dipole_field_y"];
    dipole_field_strength_[2] = floats["dipole_field_z"];
}

Wavefunction::Wavefunction(Options &options)
    : options_(options), dipole_field_strength_{{0.0, 0.0, 0.0}}, PCM_enabled_(false) {}

Wavefunction::~Wavefunction() {}

void Wavefunction::shallow_copy(SharedWavefunction other) { shallow_copy(other.get()); }

void Wavefunction::shallow_copy(const Wavefunction *other) {
    name_ = other->name_;
    basisset_ = other->basisset_;
    basissets_ = other->basissets_;
    sobasisset_ = other->sobasisset_;
    AO2SO_ = other->AO2SO_;
    S_ = other->S_;
    molecule_ = other->molecule_;

    psio_ = other->psio_;
    integral_ = other->integral_;
    factory_ = other->factory_;
    memory_ = other->memory_;
    nalpha_ = other->nalpha_;
    nbeta_ = other->nbeta_;
    nfrzc_ = other->nfrzc_;

    print_ = other->print_;
    debug_ = other->debug_;
    density_fitted_ = other->density_fitted_;

    energy_ = other->energy_;
    efzc_ = other->efzc_;

    doccpi_ = other->doccpi_;
    soccpi_ = other->soccpi_;
    frzcpi_ = other->frzcpi_;
    frzvpi_ = other->frzvpi_;
    nalphapi_ = other->nalphapi_;
    nbetapi_ = other->nbetapi_;
    nsopi_ = other->nsopi_;
    nmopi_ = other->nmopi_;

    dipole_field_type_ = other->dipole_field_type_;
    perturb_h_ = other->perturb_h_;
    std::copy(other->dipole_field_strength_.begin(), other->dipole_field_strength_.end(),
              dipole_field_strength_.begin());

    nso_ = other->nso_;
    nmo_ = other->nmo_;
    nirrep_ = other->nirrep_;

    same_a_b_dens_ = other->same_a_b_dens_;
    same_a_b_orbs_ = other->same_a_b_orbs_;

    // Set by other classes
    H_ = other->H_;
    Ca_ = other->Ca_;
    Cb_ = other->Cb_;
    Da_ = other->Da_;
    Db_ = other->Db_;
    Fa_ = other->Fa_;
    Fb_ = other->Fb_;
    epsilon_a_ = other->epsilon_a_;
    epsilon_b_ = other->epsilon_b_;

    gradient_ = other->gradient_;
    hessian_ = other->hessian_;
    external_pot_ = other->external_pot_;

    variables_ = other->variables_;
    arrays_ = other->arrays_;
    PCM_enabled_ = other->PCM_enabled_;
#ifdef USING_PCMSolver
    if (PCM_enabled_) {
        PCM_ = other->PCM_;
    }
#endif
}

void Wavefunction::deep_copy(SharedWavefunction other) { deep_copy(other.get()); }

void Wavefunction::deep_copy(const Wavefunction *other) {
    if (!other->S_) {
        throw PSIEXCEPTION("Wavefunction::deep_copy must copy an initialized wavefunction.");
    }

    /// From typical constructor
    /// Some member data is not clone-able so we will copy
    name_ = other->name_;
    molecule_ = std::make_shared<Molecule>(other->molecule_->clone());
    basisset_ = other->basisset_;
    basissets_ = other->basissets_;  // Still cannot copy basissets
    integral_ = std::make_shared<IntegralFactory>(basisset_, basisset_, basisset_, basisset_);
    sobasisset_ = std::make_shared<SOBasisSet>(basisset_, integral_);
    factory_ = std::make_shared<MatrixFactory>();
    factory_->init_with(other->nsopi_, other->nsopi_);
    AO2SO_ = other->AO2SO_->clone();
    S_ = other->S_->clone();

    psio_ = other->psio_;  // We dont actually copy psio
    memory_ = other->memory_;
    nalpha_ = other->nalpha_;
    nbeta_ = other->nbeta_;
    nfrzc_ = other->nfrzc_;

    dipole_field_type_ = other->dipole_field_type_;
    perturb_h_ = other->perturb_h_;
    std::copy(other->dipole_field_strength_.begin(), other->dipole_field_strength_.end(),
              dipole_field_strength_.begin());

    print_ = other->print_;
    debug_ = other->debug_;
    density_fitted_ = other->density_fitted_;

    energy_ = other->energy_;
    efzc_ = other->efzc_;
    variables_ = other->variables_;

    doccpi_ = other->doccpi_;
    soccpi_ = other->soccpi_;
    frzcpi_ = other->frzcpi_;
    frzvpi_ = other->frzvpi_;
    nalphapi_ = other->nalphapi_;
    nbetapi_ = other->nbetapi_;
    nsopi_ = other->nsopi_;
    nmopi_ = other->nmopi_;

    nso_ = other->nso_;
    nmo_ = other->nmo_;
    nirrep_ = other->nirrep_;

    same_a_b_dens_ = other->same_a_b_dens_;
    same_a_b_orbs_ = other->same_a_b_orbs_;

    /// Below is not set in the typical constructor
    if (other->H_) H_ = other->H_->clone();
    if (other->Ca_) Ca_ = other->Ca_->clone();
    if (other->Cb_) Cb_ = other->Cb_->clone();
    if (other->Da_) Da_ = other->Da_->clone();
    if (other->Db_) Db_ = other->Db_->clone();
    if (other->Fa_) Fa_ = other->Fa_->clone();
    if (other->Fb_) Fb_ = other->Fb_->clone();
    if (other->epsilon_a_) epsilon_a_ = SharedVector(other->epsilon_a_->clone());
    if (other->epsilon_b_) epsilon_b_ = SharedVector(other->epsilon_b_->clone());

    if (other->gradient_) gradient_ = other->gradient_->clone();
    if (other->hessian_) hessian_ = other->hessian_->clone();

    // Not sure...
    external_pot_ = other->external_pot_;

    // Copy assignment
    variables_ = other->variables_;

    // Need to explicitly call copy
    for (auto const &kv : other->arrays_) {
        arrays_[kv.first] = kv.second->clone();
    }

    PCM_enabled_ = other->PCM_enabled_;
#ifdef USING_PCMSolver
    if (PCM_enabled_) {
        PCM_ = std::make_shared<PCM>(other->PCM_.get());
    }
#endif
}

std::shared_ptr<Wavefunction> Wavefunction::c1_deep_copy(std::shared_ptr<BasisSet> basis) {
    if (!S_) {
        throw PSIEXCEPTION("Wavefunction::c1_deep_copy must copy an initialized wavefunction.");
    }

    auto wfn = std::make_shared<Wavefunction>(basis->molecule(), basis, options_);

    /// From typical constructor
    /// Some member data is not clone-able so we will copy
    wfn->name_ = name_;
    wfn->integral_ = std::make_shared<IntegralFactory>(wfn->basisset_, wfn->basisset_, wfn->basisset_, wfn->basisset_);
    wfn->sobasisset_ = std::make_shared<SOBasisSet>(wfn->basisset_, wfn->integral_);
    wfn->factory_ = std::make_shared<MatrixFactory>();

    Dimension c1_nsopi = Dimension(1);
    c1_nsopi[0] = nsopi_.sum();
    wfn->factory_->init_with(c1_nsopi, c1_nsopi);

    // Need to re-generate AO2SO_ because the new SO basis is different
    // than the old one (lower symmetry)
    auto pet = std::make_shared<PetiteList>(wfn->basisset_, wfn->integral_);
    wfn->AO2SO_ = pet->aotoso();

    wfn->psio_ = psio_;  // We dont actually copy psio
    wfn->memory_ = memory_;
    wfn->nalpha_ = nalpha_;
    wfn->nbeta_ = nbeta_;
    wfn->nfrzc_ = nfrzc_;

    wfn->print_ = print_;
    wfn->debug_ = debug_;
    wfn->density_fitted_ = density_fitted_;

    wfn->energy_ = energy_;
    wfn->efzc_ = efzc_;

    // collapse all the Dimension objects down to one element
    wfn->doccpi_.init(1, doccpi_.name());
    wfn->doccpi_[0] = doccpi_.sum();
    wfn->soccpi_.init(1, soccpi_.name());
    wfn->soccpi_[0] = soccpi_.sum();
    wfn->frzcpi_.init(1, frzcpi_.name());
    wfn->frzcpi_[0] = frzcpi_.sum();
    wfn->frzvpi_.init(1, frzvpi_.name());
    wfn->frzvpi_[0] = frzvpi_.sum();
    wfn->nalphapi_.init(1, nalphapi_.name());
    wfn->nalphapi_[0] = nalphapi_.sum();
    wfn->nbetapi_.init(1, nbetapi_.name());
    wfn->nbetapi_[0] = nbetapi_.sum();
    wfn->nsopi_.init(1, nsopi_.name());
    wfn->nsopi_[0] = nsopi_.sum();
    wfn->nmopi_.init(1, nmopi_.name());
    wfn->nmopi_[0] = nmopi_.sum();

    wfn->nso_ = nso_;
    wfn->nmo_ = nmo_;
    wfn->nirrep_ = 1;

    wfn->same_a_b_dens_ = same_a_b_dens_;
    wfn->same_a_b_orbs_ = same_a_b_orbs_;

    /// Need the SO2AO matrix for remove_symmetry(), have the AO2SO matrix
    SharedMatrix SO2AO = aotoso()->transpose();

    wfn->S_ = wfn->factory_->create_shared_matrix("S");
    wfn->S_->remove_symmetry(S_, SO2AO);

    /// Below is not set in the typical constructor

    wfn->H_ = wfn->factory_->create_shared_matrix("One-electron Hamiltonian");
    wfn->H_->remove_symmetry(H_, SO2AO);

    /* This stuff we need to copy in the subclass functions, b/c
    ** constructors like RHF() just blow these away anyway
    if (Ca_) wfn->Ca_ = Ca_subset("AO", "ALL");
    if (Cb_) wfn->Cb_ = Cb_subset("AO", "ALL");
    if (Da_) wfn->Da_ = Da_subset("AO");
    if (Db_) wfn->Db_ = Db_subset("AO");
    if (Fa_) wfn->Fa_ = Fa_subset("AO");
    if (Fb_) wfn->Fb_ = Fb_subset("AO");
    if (epsilon_a_) wfn->epsilon_a_ =
        epsilon_subset_helper(epsilon_a_, nsopi_, "AO", "ALL");
    if (epsilon_b_) wfn->epsilon_b_ =
        epsilon_subset_helper(epsilon_b_, nsopi_, "AO", "ALL");
    */

    // these are simple SharedMatrices of size 3*natom_, etc., so should
    // not depend on symmetry ... can just copy them
    if (gradient_) wfn->gradient_ = gradient_->clone();
    if (hessian_) wfn->hessian_ = hessian_->clone();

    wfn->variables_ = variables_;
    // ok for deep copy?
    wfn->external_pot_ = external_pot_;

    // Need to explicitly call copy
    for (auto const &kv : arrays_) {
        wfn->arrays_[kv.first] = kv.second->clone();
    }

    wfn->PCM_enabled_ = PCM_enabled_;
#ifdef USING_PCMSolver
    if (wfn->PCM_enabled_) {
        wfn->PCM_ = std::make_shared<PCM>(PCM_.get());
    }
#endif

    return wfn;
}

void Wavefunction::common_init() {
    Wavefunction::initialize_singletons();
    if (!basisset_) {
        throw PSIEXCEPTION(
            "You can't initialize a Wavefunction that doesn't "
            "have a basis set");
    }

    // Check the point group of the molecule. If it is not set, set it.
    if (!molecule_->point_group()) {
        molecule_->set_point_group(molecule_->find_point_group());
    }

    // Create an SO basis...we need the point group for this part.
    integral_ = std::make_shared<IntegralFactory>(basisset_, basisset_, basisset_, basisset_);
    sobasisset_ = std::make_shared<SOBasisSet>(basisset_, integral_);

    auto pet = std::make_shared<PetiteList>(basisset_, integral_);
    AO2SO_ = pet->aotoso();

    // Obtain the dimension object to initialize the factory.
    nsopi_ = sobasisset_->dimension();
    nsopi_.set_name("SOs per irrep");

    factory_ = std::make_shared<MatrixFactory>();
    factory_->init_with(nsopi_, nsopi_);

    nirrep_ = nsopi_.n();

    S_ = factory_->create_shared_matrix("S");
    std::shared_ptr<OneBodySOInt> Sint(integral_->so_overlap());
    Sint->compute(S_);

    // Initialize array that hold dimensionality information
    nmopi_ = Dimension(nirrep_, "MOs per irrep");
    nalphapi_ = Dimension(nirrep_, "Alpha electrons per irrep");
    nbetapi_ = Dimension(nirrep_, "Beta electrons per irrep");
    doccpi_ = Dimension(nirrep_, "Doubly occupied orbitals per irrep");
    soccpi_ = Dimension(nirrep_, "Singly occupied orbitals per irrep");
    frzcpi_ = Dimension(nirrep_, "Frozen core orbitals per irrep");
    frzvpi_ = Dimension(nirrep_, "Frozen virtual orbitals per irrep");

    // Obtain memory amount from the environment
    memory_ = Process::environment.get_memory();

    nso_ = basisset_->nbf();
    nmo_ = basisset_->nbf();
    for (int k = 0; k < nirrep_; k++) {
        nmopi_[k] = 0;
        doccpi_[k] = 0;
        soccpi_[k] = 0;
        nalphapi_[k] = 0;
        nbetapi_[k] = 0;
    }

    density_fitted_ = false;
    energy_ = 0.0;
    efzc_ = 0.0;
    same_a_b_dens_ = true;
    same_a_b_orbs_ = false;

    // Read in the debug flag
    debug_ = options_.get_int("DEBUG");
    print_ = options_.get_int("PRINT");

    // Determine the number of electrons in the system
    int nelectron = 0;
    for (int i = 0; i < molecule_->natom(); ++i) {
        nelectron += (int)molecule_->Z(i);
    }
    nelectron -= molecule_->molecular_charge();

    // Make sure that the multiplicity is reasonable
    int multiplicity = molecule_->multiplicity();
    if (multiplicity - 1 > nelectron) {
        char *str = new char[200];
        sprintf(str,
                "There are not enough electrons for multiplicity = %d.\n"
                "Please check your input",
                multiplicity);
        throw SanityCheckError(str, __FILE__, __LINE__);
        delete[] str;
    }
    if (multiplicity % 2 == nelectron % 2) {
        char *str = new char[200];
        sprintf(str,
                "A multiplicity of %d with %d electrons is impossible.\n"
                "Please check your input",
                multiplicity, nelectron);
        throw SanityCheckError(str, __FILE__, __LINE__);
        delete[] str;
    }

    nbeta_ = (nelectron - multiplicity + 1) / 2;
    nalpha_ = nbeta_ + multiplicity - 1;

    // Setup dipole field perturbation information
    perturb_h_ = options_.get_bool("PERTURB_H");
    dipole_field_type_ = nothing;
    std::fill(dipole_field_strength_.begin(), dipole_field_strength_.end(), 0.0);
    if (perturb_h_) {
        std::string perturb_with;
        if (options_["PERTURB_WITH"].has_changed()) {
            perturb_with = options_.get_str("PERTURB_WITH");
            // Do checks to see what perturb_with is.
            if (perturb_with == "DIPOLE_X") {
                dipole_field_type_ = dipole_x;
                dipole_field_strength_[0] = options_.get_double("PERTURB_MAGNITUDE");
                outfile->Printf(
                    " WARNING: the DIPOLE_X and PERTURB_MAGNITUDE keywords are deprecated."
                    "  Use DIPOLE and the PERTURB_DIPOLE array instead.");
            } else if (perturb_with == "DIPOLE_Y") {
                dipole_field_type_ = dipole_y;
                dipole_field_strength_[1] = options_.get_double("PERTURB_MAGNITUDE");
                outfile->Printf(
                    " WARNING: the DIPOLE_Y and PERTURB_MAGNITUDE keywords are deprecated."
                    "  Use DIPOLE and the PERTURB_DIPOLE array instead.");
            } else if (perturb_with == "DIPOLE_Z") {
                dipole_field_type_ = dipole_z;
                dipole_field_strength_[2] = options_.get_double("PERTURB_MAGNITUDE");
                outfile->Printf(
                    " WARNING: the DIPOLE_Z and PERTURB_MAGNITUDE keywords are deprecated."
                    "  Use DIPOLE and the PERTURB_DIPOLE array instead.");
            } else if (perturb_with == "DIPOLE") {
                dipole_field_type_ = dipole;
                if (options_["PERTURB_DIPOLE"].size() != 3)
                    throw PSIEXCEPTION("The PERTURB dipole should have exactly three floating point numbers.");
                for (int n = 0; n < 3; ++n) dipole_field_strength_[n] = options_["PERTURB_DIPOLE"][n].to_double();
            } else if (perturb_with == "EMBPOT") {
                dipole_field_type_ = embpot;
            } else if (perturb_with == "DX") {
                dipole_field_type_ = dx;
            } else if (perturb_with == "SPHERE") {
                dipole_field_type_ = sphere;
            } else {
                outfile->Printf("Unknown PERTURB_WITH. Applying no perturbation.\n");
            }
        } else {
            outfile->Printf("PERTURB_H is true, but PERTURB_WITH not found, applying no perturbation.\n");
        }
    }
}

std::array<double, 3> Wavefunction::get_dipole_field_strength() const { return dipole_field_strength_; }

Wavefunction::FieldType Wavefunction::get_dipole_perturbation_type() const { return dipole_field_type_; }

Dimension Wavefunction::map_irreps(const Dimension &dimpi) {
    auto ps = options_.get_str("PARENT_SYMMETRY");

    // If the parent symmetry hasn't been set, no displacements have been made
    if (ps == "") return dimpi;

    auto full = std::make_shared<PointGroup> (ps);
    std::shared_ptr<PointGroup> sub = molecule_->point_group();

    // If the point group between the full and sub are the same return
    if (full->symbol() == sub->symbol()) return dimpi;

    // Build the correlation table between full, and subgroup
    CorrelationTable corrtab(full, sub);
    Dimension mapped_dimpi(sub->char_table().nirrep());
    for (int h = 0; h < full->char_table().nirrep(); ++h) {
        int target = corrtab.gamma(h, 0);
        mapped_dimpi[target] += dimpi[h];
    }

    return mapped_dimpi;
}

void Wavefunction::initialize_singletons() {
    static bool done = false;

    if (done) return;

    ioff[0] = 0;
    for (size_t i = 1; i < MAX_IOFF; ++i) {
        ioff[i] = ioff[i - 1] + i;
    }

    df[0] = 1.0;
    df[1] = 1.0;
    df[2] = 1.0;
    for (int i = 3; i < MAX_DF; ++i) {
        df[i] = (i - 1) * df[i - 2];
    }

    for (int i = 0; i < MAX_BC; ++i) {
        for (int j = 0; j <= i; ++j) {
            bc[i][j] = combinations(i, j);
        }
    }

    fac[0] = 1.0;
    for (int i = 1; i < MAX_FAC; ++i) {
        fac[i] = i * fac[i - 1];
    }

    done = true;
}

std::shared_ptr<Molecule> Wavefunction::molecule() const { return molecule_; }

std::shared_ptr<PSIO> Wavefunction::psio() const { return psio_; }

Options &Wavefunction::options() const { return options_; }

std::shared_ptr<IntegralFactory> Wavefunction::integral() const { return integral_; }

std::shared_ptr<BasisSet> Wavefunction::basisset() const { return basisset_; }

std::shared_ptr<BasisSet> Wavefunction::get_basisset(std::string label) {
    // This may be slightly confusing, but better than changing this in 800 other places
    if (label == "ORBITAL") {
        return basisset_;
    } else if (basissets_.count(label) == 0) {
        outfile->Printf("Could not find requested basisset (%s).", label.c_str());
        throw PSIEXCEPTION("Wavefunction::get_basisset: Requested basis set (" + label + ") was not set!\n");
    } else {
        return basissets_[label];
    }
}
void Wavefunction::set_basisset(std::string label, std::shared_ptr<BasisSet> basis) {
    if (label == "ORBITAL") {
        throw PSIEXCEPTION("Cannot set the ORBITAL basis after the Wavefunction is built!");
    } else {
        basissets_[label] = basis;
    }
}

bool Wavefunction::basisset_exists(std::string label) {
    if (basissets_.count(label) == 0)
        return false;
    else
        return true;
}

std::shared_ptr<SOBasisSet> Wavefunction::sobasisset() const { return sobasisset_; }

std::shared_ptr<MatrixFactory> Wavefunction::matrix_factory() const { return factory_; }

std::shared_ptr<Wavefunction> Wavefunction::reference_wavefunction() const { return reference_wavefunction_; }

void Wavefunction::set_reference_wavefunction(const std::shared_ptr<Wavefunction> wfn) {
    reference_wavefunction_ = wfn;
}

void Wavefunction::force_doccpi(const Dimension &doccpi) {
    for (int h = 0; h < nirrep_; h++) {
        if ((soccpi_[h] + doccpi[h]) > nmopi_[h]) {
            throw PSIEXCEPTION(
                "Wavefunction::force_doccpi: Number of doubly and singly occupied orbitals in an irrep cannot exceed "
                "the total number of molecular orbitals.");
        }
        doccpi_[h] = doccpi[h];
        nalphapi_[h] = doccpi_[h] + soccpi_[h];
        nbetapi_[h] = doccpi_[h];
    }
    nalpha_ = doccpi_.sum() + soccpi_.sum();
    nbeta_ = doccpi_.sum();
}

void Wavefunction::force_soccpi(const Dimension &soccpi) {
    if (same_a_b_dens_) {
        throw PSIEXCEPTION(
            "Wavefunction::force_soccpi: Cannot set soccpi since alpha and beta densities must be the same for this "
            "Wavefunction.");
    }
    for (int h = 0; h < nirrep_; h++) {
        if ((soccpi[h] + doccpi_[h]) > nmopi_[h]) {
            throw PSIEXCEPTION(
                "Wavefunction::force_soccpi: Number of doubly and singly occupied orbitals in an irrep cannot exceed "
                "the total number of molecular orbitals.");
        }
        soccpi_[h] = soccpi[h];
        nalphapi_[h] = doccpi_[h] + soccpi_[h];
    }
    nalpha_ = doccpi_.sum() + soccpi_.sum();
}

void Wavefunction::set_frzvpi(const Dimension &frzvpi) {
    for (int h = 0; h < nirrep_; h++) {
        frzvpi_[h] = frzvpi[h];
    }
}

SharedMatrix Wavefunction::Ca() const {
    if (!Ca_) {
        if (!reference_wavefunction_)
            throw PSIEXCEPTION("Wavefunction::Ca: Unable to obtain MO coefficients.");
        else
            return reference_wavefunction_->Ca();
    }

    return Ca_;
}

SharedMatrix Wavefunction::Cb() const {
    if (!Cb_) {
        if (!reference_wavefunction_)
            throw PSIEXCEPTION("Wavefunction::Cb: Unable to obtain MO coefficients.");
        else
            return reference_wavefunction_->Cb();
    }

    return Cb_;
}

std::vector<std::vector<int>> Wavefunction::subset_occupation(const Dimension &noccpi,
                                                              const std::string &subset) const {
    if (!(subset == "FROZEN_OCC" || subset == "FROZEN_VIR" || subset == "ACTIVE_OCC" || subset == "ACTIVE_VIR" ||
          subset == "FROZEN" || subset == "ACTIVE" || subset == "OCC" || subset == "VIR" || subset == "ALL"))
        throw PSIEXCEPTION(
            "Orbital subset is not defined, should be ALL, OCC, VIR, FROZEN, ACTIVE, or look at doxygen");

    // Vector of relevent positions by irrep
    std::vector<std::vector<int>> positions;
    positions.resize(nirrep_);

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h]; i++) {
            if (subset == "ALL" || subset == "FROZEN_OCC" || subset == "FROZEN" || subset == "OCC")
                positions[h].push_back(i);
        }
        for (int i = frzcpi_[h]; i < noccpi[h]; i++) {
            if (subset == "ALL" || subset == "ACTIVE_OCC" || subset == "ACTIVE" || subset == "OCC")
                positions[h].push_back(i);
        }
        for (int i = noccpi[h]; i < nmopi_[h] - frzvpi_[h]; i++) {
            if (subset == "ALL" || subset == "ACTIVE_VIR" || subset == "ACTIVE" || subset == "VIR")
                positions[h].push_back(i);
        }
        for (int i = nmopi_[h] - frzvpi_[h]; i < nmopi_[h]; i++) {
            if (subset == "ALL" || subset == "FROZEN_VIR" || subset == "FROZEN" || subset == "VIR")
                positions[h].push_back(i);
        }
    }
    return positions;
}

SharedMatrix Wavefunction::C_subset_helper(SharedMatrix C, const Dimension &noccpi, SharedVector epsilon,
                                           const std::string &basis, const std::string &subset) const {
    std::vector<std::vector<int>> positions = subset_occupation(noccpi, subset);

    Dimension nmopi(nirrep_);
    for (int h = 0; h < (int)positions.size(); h++) {
        nmopi[h] = positions[h].size();
    }
    auto C2 = std::make_shared<Matrix>("C " + basis + " " + subset, nsopi_, nmopi);
    for (int h = 0; h < (int)positions.size(); h++) {
        for (int i = 0; i < (int)positions[h].size(); i++) {
            C_DCOPY(nsopi_[h], &C->pointer(h)[0][positions[h][i]], nmopi_[h], &C2->pointer(h)[0][i], nmopi[h]);
        }
    }

    if (basis == "AO") {
        auto C3 = std::make_shared<Matrix>("C " + basis + " " + subset, nso_, nmopi.sum());
        std::swap(C2, C3);

        std::vector<std::tuple<double, int, int>> order;
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < (int)positions[h].size(); i++) {
                order.push_back(std::tuple<double, int, int>(epsilon->get(h, positions[h][i]), i, h));
            }
        }

        std::sort(order.begin(), order.end(), std::less<std::tuple<double, int, int>>());

        for (int index = 0; index < (int)order.size(); index++) {
            int i = std::get<1>(order[index]);
            int h = std::get<2>(order[index]);

            int nao = nso_;
            int nso = nsopi_[h];

            if (!nso) continue;

            C_DGEMV('N', nao, nso, 1.0, AO2SO_->pointer(h)[0], nso, &C3->pointer(h)[0][i], nmopi[h], 0.0,
                    &C2->pointer()[0][index], nmopi.sum());
        }

    } else if (basis == "SO" || basis == "MO") {
        // Already done
    } else {
        throw PSIEXCEPTION("Invalid basis requested, use AO, SO, or MO");
    }

    return C2;
}

SharedVector Wavefunction::epsilon_subset_helper(SharedVector epsilon, const Dimension &noccpi,
                                                 const std::string &basis, const std::string &subset) const {
    std::vector<std::vector<int>> positions = subset_occupation(noccpi, subset);

    Dimension nmopi(nirrep_);
    for (int h = 0; h < (int)positions.size(); h++) {
        nmopi[h] = positions[h].size();
    }

    SharedVector C2;

    if (basis == "AO") {
        C2 = std::make_shared<Vector>("Epsilon " + basis + " " + subset, nmopi.sum());

        std::vector<std::tuple<double, int, int>> order;
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < (int)positions[h].size(); i++) {
                order.push_back(std::tuple<double, int, int>(epsilon->get(h, positions[h][i]), i, h));
            }
        }

        std::sort(order.begin(), order.end(), std::less<std::tuple<double, int, int>>());

        for (int index = 0; index < (int)order.size(); index++) {
            C2->set(0, index, std::get<0>(order[index]));
        }

    } else if (basis == "SO" || basis == "MO") {
        C2 = std::make_shared<Vector>("Epsilon " + basis + " " + subset, nmopi);
        for (int h = 0; h < (int)positions.size(); h++) {
            for (int i = 0; i < (int)positions[h].size(); i++) {
                C2->set(h, i, epsilon->get(h, positions[h][i]));
            }
        }

    } else {
        throw PSIEXCEPTION("Invalid basis requested, use AO, SO, or MO");
    }

    return C2;
}

SharedMatrix Wavefunction::F_subset_helper(SharedMatrix F, SharedMatrix C, const std::string &basis) const {
    if (basis == "AO") {
        double *temp = new double[AO2SO_->max_ncol() * AO2SO_->max_nrow()];
        auto F2 = std::make_shared<Matrix>("Fock (AO basis)", basisset_->nbf(), basisset_->nbf());
        int symm = F->symmetry();
        for (int h = 0; h < AO2SO_->nirrep(); ++h) {
            int nao = AO2SO_->rowspi()[0];
            int nsol = AO2SO_->colspi()[h];
            int nsor = AO2SO_->colspi()[h ^ symm];
            if (!nsol || !nsor) continue;
            double **Ulp = AO2SO_->pointer(h);
            double **Urp = AO2SO_->pointer(h ^ symm);
            double **FSOp = F->pointer(h);
            double **FAOp = F2->pointer();
            C_DGEMM('N', 'T', nsol, nao, nsor, 1.0, FSOp[0], nsor, Urp[0], nsor, 0.0, temp, nao);
            C_DGEMM('N', 'N', nao, nao, nsol, 1.0, Ulp[0], nsol, temp, nao, 1.0, FAOp[0], nao);
        }
        delete[] temp;
        return F2;
    } else if (basis == "SO") {
        return SharedMatrix(F->clone());
    } else if (basis == "MO") {
        auto F2 = std::make_shared<Matrix>("Fock (MO Basis)", C->colspi(), C->colspi());

        int symm = F->symmetry();
        int nirrep = F->nirrep();

        double *SC = new double[C->max_ncol() * C->max_nrow()];
        double *temp = new double[C->max_ncol() * C->max_nrow()];
        for (int h = 0; h < nirrep; h++) {
            int nmol = C->colspi()[h];
            int nmor = C->colspi()[h ^ symm];
            int nsol = C->rowspi()[h];
            int nsor = C->rowspi()[h ^ symm];
            if (!nmol || !nmor || !nsol || !nsor) continue;
            double **Slp = S_->pointer(h);
            double **Srp = S_->pointer(h ^ symm);
            double **Clp = C->pointer(h);
            double **Crp = C->pointer(h ^ symm);
            double **Fmop = F2->pointer(h);
            double **Fsop = F->pointer(h);

            C_DGEMM('N', 'N', nsor, nmor, nsor, 1.0, Srp[0], nsor, Crp[0], nmor, 0.0, SC, nmor);
            C_DGEMM('N', 'N', nsol, nmor, nsor, 1.0, Fsop[0], nsor, SC, nmor, 0.0, temp, nmor);
            C_DGEMM('N', 'N', nsol, nmol, nsol, 1.0, Slp[0], nsol, Clp[0], nmol, 0.0, SC, nmol);
            C_DGEMM('T', 'N', nmol, nmor, nsol, 1.0, SC, nmol, temp, nmor, 0.0, Fmop[0], nmor);
        }
        delete[] temp;
        delete[] SC;
        return F2;
    } else {
        throw PSIEXCEPTION("Invalid basis requested, use AO, SO, or MO");
    }
}

SharedMatrix Wavefunction::D_subset_helper(SharedMatrix D, SharedMatrix C, const std::string &basis) const {
    if (basis == "AO") {
        double *temp = new double[AO2SO_->max_ncol() * AO2SO_->max_nrow()];
        auto D2 = std::make_shared<Matrix>("D (AO basis)", basisset_->nbf(), basisset_->nbf());
        int symm = D->symmetry();
        for (int h = 0; h < AO2SO_->nirrep(); ++h) {
            int nao = AO2SO_->rowspi()[0];
            int nsol = AO2SO_->colspi()[h];
            int nsor = AO2SO_->colspi()[h ^ symm];
            if (!nsol || !nsor) continue;
            double **Ulp = AO2SO_->pointer(h);
            double **Urp = AO2SO_->pointer(h ^ symm);
            double **DSOp = D->pointer(h);
            double **DAOp = D2->pointer();
            C_DGEMM('N', 'T', nsol, nao, nsor, 1.0, DSOp[0], nsor, Urp[0], nsor, 0.0, temp, nao);
            C_DGEMM('N', 'N', nao, nao, nsol, 1.0, Ulp[0], nsol, temp, nao, 1.0, DAOp[0], nao);
        }
        delete[] temp;
        return D2;
    } else if (basis == "CartAO") {
        /*
         * Added by ACS. Rob's definition of AO is simply a desymmetrized SO (i.e. using spherical basis
         * functions).  In cases like EFP and PCM where many OE integral evaluations are needed, we want
         * to avoid the spherical transformation, so we need to back transform the density matrix all the
         * way back to Cartesian AOs.
         */

        PetiteList petite(basisset_, integral_, true);
        SharedMatrix my_aotoso = petite.aotoso();
        double *temp = new double[my_aotoso->max_ncol() * my_aotoso->max_nrow()];
        auto D2 = std::make_shared<Matrix>("D (ao basis)", basisset_->nao(), basisset_->nao());
        int symm = D->symmetry();
        for (int h = 0; h < my_aotoso->nirrep(); ++h) {
            int nao = my_aotoso->rowspi()[0];
            int nsol = my_aotoso->colspi()[h];
            int nsor = my_aotoso->colspi()[h ^ symm];
            if (!nsol || !nsor) continue;
            double **Ulp = my_aotoso->pointer(h);
            double **Urp = my_aotoso->pointer(h ^ symm);
            double **DSOp = D->pointer(h);
            double **DAOp = D2->pointer();
            C_DGEMM('N', 'T', nsol, nao, nsor, 1.0, DSOp[0], nsor, Urp[0], nsor, 0.0, temp, nao);
            C_DGEMM('N', 'N', nao, nao, nsol, 1.0, Ulp[0], nsol, temp, nao, 1.0, DAOp[0], nao);
        }
        delete[] temp;
        return D2;
    } else if (basis == "SO") {
        return SharedMatrix(D->clone());
    } else if (basis == "MO") {
        auto D2 = std::make_shared<Matrix>("D (MO Basis)", C->colspi(), C->colspi());

        int symm = D->symmetry();
        int nirrep = D->nirrep();

        double *SC = new double[C->max_ncol() * C->max_nrow()];
        double *temp = new double[C->max_ncol() * C->max_nrow()];
        for (int h = 0; h < nirrep; h++) {
            int nmol = C->colspi()[h];
            int nmor = C->colspi()[h ^ symm];
            int nsol = C->rowspi()[h];
            int nsor = C->rowspi()[h ^ symm];
            if (!nmol || !nmor || !nsol || !nsor) continue;
            double **Slp = S_->pointer(h);
            double **Srp = S_->pointer(h ^ symm);
            double **Clp = C->pointer(h);
            double **Crp = C->pointer(h ^ symm);
            double **Dmop = D2->pointer(h);
            double **Dsop = D->pointer(h);

            C_DGEMM('N', 'N', nsor, nmor, nsor, 1.0, Srp[0], nsor, Crp[0], nmor, 0.0, SC, nmor);
            C_DGEMM('N', 'N', nsol, nmor, nsor, 1.0, Dsop[0], nsor, SC, nmor, 0.0, temp, nmor);
            C_DGEMM('N', 'N', nsol, nmol, nsol, 1.0, Slp[0], nsol, Clp[0], nmol, 0.0, SC, nmol);
            C_DGEMM('T', 'N', nmol, nmor, nsol, 1.0, SC, nmol, temp, nmor, 0.0, Dmop[0], nmor);
        }
        delete[] temp;
        delete[] SC;
        return D2;
    } else {
        throw PSIEXCEPTION("Invalid basis requested, use AO, CartAO, SO, or MO");
    }
}

SharedMatrix Wavefunction::matrix_subset_helper(SharedMatrix M, SharedMatrix C, const std::string &basis,
                                                const std::string matrix_basename) const {
    if (basis == "AO") {
        double *temp = new double[AO2SO_->max_ncol() * AO2SO_->max_nrow()];
        std::string m2_name = matrix_basename + " (AO basis)";
        SharedMatrix M2 = SharedMatrix(new Matrix(m2_name, basisset_->nbf(), basisset_->nbf()));
        int symm = M->symmetry();
        for (int h = 0; h < AO2SO_->nirrep(); ++h) {
            int nao = AO2SO_->rowspi()[0];
            int nsol = AO2SO_->colspi()[h];
            int nsor = AO2SO_->colspi()[h ^ symm];
            if (!nsol || !nsor) continue;
            double **Ulp = AO2SO_->pointer(h);
            double **Urp = AO2SO_->pointer(h ^ symm);
            double **MSOp = M->pointer(h);
            double **MAOp = M2->pointer();
            C_DGEMM('N', 'T', nsol, nao, nsor, 1.0, MSOp[0], nsor, Urp[0], nsor, 0.0, temp, nao);
            C_DGEMM('N', 'N', nao, nao, nsol, 1.0, Ulp[0], nsol, temp, nao, 1.0, MAOp[0], nao);
        }
        delete[] temp;
        return M2;
    } else if (basis == "CartAO") {
        /*
         * Added by ACS. Rob's definition of AO is simply a desymmetrized SO (i.e. using spherical basis
         * functions).  In cases like EFP and PCM where many OE integral evaluations are needed, we want
         * to avoid the spherical transformation, so we need to back transform the density matrix all the
         * way back to Cartesian AOs.
         */

        PetiteList petite(basisset_, integral_, true);
        SharedMatrix my_aotoso = petite.aotoso();
        double *temp = new double[my_aotoso->max_ncol() * my_aotoso->max_nrow()];
        std::string m2_name = matrix_basename + " (CartAO basis)";
        SharedMatrix M2 = SharedMatrix(new Matrix(m2_name, basisset_->nao(), basisset_->nao()));
        int symm = M->symmetry();
        for (int h = 0; h < my_aotoso->nirrep(); ++h) {
            int nao = my_aotoso->rowspi()[0];
            int nsol = my_aotoso->colspi()[h];
            int nsor = my_aotoso->colspi()[h ^ symm];
            if (!nsol || !nsor) continue;
            double **Ulp = my_aotoso->pointer(h);
            double **Urp = my_aotoso->pointer(h ^ symm);
            double **MSOp = M->pointer(h);
            double **MAOp = M2->pointer();
            C_DGEMM('N', 'T', nsol, nao, nsor, 1.0, MSOp[0], nsor, Urp[0], nsor, 0.0, temp, nao);
            C_DGEMM('N', 'N', nao, nao, nsol, 1.0, Ulp[0], nsol, temp, nao, 1.0, MAOp[0], nao);
        }
        delete[] temp;
        return M2;
    } else if (basis == "SO") {
        SharedMatrix M2 = M->clone();
        std::string m2_name = matrix_basename + " (SO basis)";
        M2->set_name(m2_name);
        return M2;
    } else if (basis == "MO") {
        std::string m2_name = matrix_basename + " (MO basis)";
        SharedMatrix M2(new Matrix(m2_name, C->colspi(), C->colspi()));

        int symm = M->symmetry();
        int nirrep = M->nirrep();

        double *SC = new double[C->max_ncol() * C->max_nrow()];
        double *temp = new double[C->max_ncol() * C->max_nrow()];
        for (int h = 0; h < nirrep; h++) {
            int nmol = C->colspi()[h];
            int nmor = C->colspi()[h ^ symm];
            int nsol = C->rowspi()[h];
            int nsor = C->rowspi()[h ^ symm];
            if (!nmol || !nmor || !nsol || !nsor) continue;
            double **Slp = S_->pointer(h);
            double **Srp = S_->pointer(h ^ symm);
            double **Clp = C->pointer(h);
            double **Crp = C->pointer(h ^ symm);
            double **Mmop = M2->pointer(h);
            double **Msop = M->pointer(h);

            C_DGEMM('N', 'N', nsor, nmor, nsor, 1.0, Srp[0], nsor, Crp[0], nmor, 0.0, SC, nmor);
            C_DGEMM('N', 'N', nsol, nmor, nsor, 1.0, Msop[0], nsor, SC, nmor, 0.0, temp, nmor);
            C_DGEMM('N', 'N', nsol, nmol, nsol, 1.0, Slp[0], nsol, Clp[0], nmol, 0.0, SC, nmol);
            C_DGEMM('T', 'N', nmol, nmor, nsol, 1.0, SC, nmol, temp, nmor, 0.0, Mmop[0], nmor);
        }
        delete[] temp;
        delete[] SC;
        return M2;
    } else {
        throw PSIEXCEPTION("Invalid basis requested, use AO, CartAO, SO, or MO");
    }
}

SharedMatrix Wavefunction::basis_projection(SharedMatrix C_A, Dimension noccpi, std::shared_ptr<BasisSet> old_basis,
                                            std::shared_ptr<BasisSet> new_basis) {
    // Based on Werner's method from Mol. Phys. 102, 21-22, 2311
    std::shared_ptr<IntegralFactory> newfactory =
        std::make_shared<IntegralFactory>(new_basis, new_basis, new_basis, new_basis);
    std::shared_ptr<IntegralFactory> hybfactory =
        std::make_shared<IntegralFactory>(old_basis, new_basis, old_basis, new_basis);
    std::shared_ptr<OneBodySOInt> intBB(newfactory->so_overlap());
    std::shared_ptr<OneBodySOInt> intAB(hybfactory->so_overlap());

    auto pet = std::make_shared<PetiteList>(new_basis, newfactory);
    SharedMatrix AO2USO(pet->aotoso());

    auto SAB = std::make_shared<Matrix>("S_AB", C_A->nirrep(), C_A->rowspi(), AO2USO->colspi());
    auto SBB = std::make_shared<Matrix>("S_BB", C_A->nirrep(), AO2USO->colspi(), AO2USO->colspi());

    intAB->compute(SAB);
    intBB->compute(SBB);

    // C_A->print();
    // noccpi.print();
    // SAB->print();
    // SBB->print();

    newfactory.reset();
    hybfactory.reset();
    intAB.reset();
    intBB.reset();
    pet.reset();

    // Constrained to the same symmetry at the moment, we can relax this soon
    auto C_B = std::make_shared<Matrix>("C_B", C_A->nirrep(), AO2USO->colspi(), noccpi);

    // Block over irreps (soon united irreps)
    for (int h = 0; h < C_A->nirrep(); h++) {
        int nocc = noccpi[h];
        int na = C_A->rowspi()[h];
        int nb = AO2USO->colspi()[h];
        int nc = C_A->colspi()[h];

        if (nocc == 0 || na == 0 || nb == 0) continue;

        auto Ca = C_A->pointer(h);
        auto Cb = C_B->pointer(h);
        auto Sab = SAB->pointer(h);
        auto Sbb = SBB->pointer(h);

        int CholError = C_DPOTRF('L', nb, Sbb[0], nb);
        if (CholError != 0) throw std::domain_error("S_BB Matrix Cholesky failed!");

        // Inversion (in place)
        int IError = C_DPOTRI('L', nb, Sbb[0], nb);
        if (IError != 0) throw std::domain_error("S_BB Inversion Failed!");

        // LAPACK is smart and all, only uses half of the thing
        for (int m = 0; m < nb; m++)
            for (int n = 0; n < m; n++) Sbb[m][n] = Sbb[n][m];

        // Form T
        double **Temp1 = block_matrix(nb, nocc);
        C_DGEMM('T', 'N', nb, nocc, na, 1.0, Sab[0], nb, Ca[0], nc, 0.0, Temp1[0], nocc);

        // outfile->Printf(" Temp1:\n");
        // print_mat(Temp1,nb,nocc,"outfile");

        double **Temp2 = block_matrix(nb, nocc);
        C_DGEMM('N', 'N', nb, nocc, nb, 1.0, Sbb[0], nb, Temp1[0], nocc, 0.0, Temp2[0], nocc);

        // outfile->Printf(" Temp2:\n");
        // print_mat(Temp2,nb,nocc,outfile);

        double **Temp3 = block_matrix(na, nocc);
        C_DGEMM('N', 'N', na, nocc, nb, 1.0, Sab[0], nb, Temp2[0], nocc, 0.0, Temp3[0], nocc);

        // outfile->Printf(" Temp3:\n");
        // print_mat(Temp3,na,nocc,outfile);

        double **T = block_matrix(nocc, nocc);
        C_DGEMM('T', 'N', nocc, nocc, na, 1.0, Ca[0], nc, Temp3[0], nocc, 0.0, T[0], nocc);

        // outfile->Printf(" T:\n");
        // print_mat(T,nocc,nocc,"outfile");

        // Find T^-1/2
        // First, diagonalize T
        // the C_DSYEV call replaces the original matrix T with its eigenvectors
        double *eigval = init_array(nocc);
        int lwork = nocc * 3;
        double *work = init_array(lwork);
        int stat = C_DSYEV('v', 'u', nocc, T[0], nocc, eigval, work, lwork);
        if (stat != 0) {
            outfile->Printf("C_DSYEV failed\n");
            exit(PSI_RETURN_FAILURE);
        }
        delete[] work;

        // Now T contains the eigenvectors of the original T
        // Copy T to T_copy
        double **T_mhalf = block_matrix(nocc, nocc);
        double **T_copy = block_matrix(nocc, nocc);
        C_DCOPY(nocc * nocc, T[0], 1, T_copy[0], 1);

        // Now form T^{-1/2} = U(T)*t^{-1/2}*U,
        // where t^{-1/2} is the diagonal matrix of the inverse square roots
        // of the eigenvalues, and U is the matrix of eigenvectors of T
        for (int i = 0; i < nocc; i++) {
            if (eigval[i] < 1.0E-10)
                eigval[i] = 0.0;
            else
                eigval[i] = 1.0 / sqrt(eigval[i]);

            // scale one set of eigenvectors by the diagonal elements t^{-1/2}
            C_DSCAL(nocc, eigval[i], T[i], 1);
        }
        free(eigval);

        // T_mhalf = T_copy(T) * T
        C_DGEMM('t', 'n', nocc, nocc, nocc, 1.0, T_copy[0], nocc, T[0], nocc, 0.0, T_mhalf[0], nocc);

        // Form CB
        C_DGEMM('N', 'N', nb, nocc, nocc, 1.0, Temp2[0], nocc, T_mhalf[0], nocc, 0.0, Cb[0], nocc);

        free_block(Temp1);
        free_block(Temp2);
        free_block(Temp3);
        free_block(T);
        free_block(T_copy);
        free_block(T_mhalf);
    }
    return C_B;
}

OrbitalSpace Wavefunction::alpha_orbital_space(const std::string &id, const std::string &basis,
                                               const std::string &subset) {
    return OrbitalSpace(id, subset, Ca_subset(basis, subset), epsilon_a_subset(basis, subset), basisset_, integral_);
}

OrbitalSpace Wavefunction::beta_orbital_space(const std::string &id, const std::string &basis,
                                              const std::string &subset) {
    return OrbitalSpace(id, subset, Cb_subset(basis, subset), epsilon_b_subset(basis, subset), basisset_, integral_);
}

SharedMatrix Wavefunction::Ca_subset(const std::string &basis, const std::string &subset) const {
    return C_subset_helper(Ca_, nalphapi_, epsilon_a_, basis, subset);
}

SharedMatrix Wavefunction::Cb_subset(const std::string &basis, const std::string &subset) const {
    return C_subset_helper(Cb_, nbetapi_, epsilon_b_, basis, subset);
}

SharedMatrix Wavefunction::Da_subset(const std::string &basis) const {
    return matrix_subset_helper(Da_, Ca_, basis, "D");
}

SharedMatrix Wavefunction::Db_subset(const std::string &basis) const {
    return matrix_subset_helper(Db_, Cb_, basis, "D");
}

SharedMatrix Wavefunction::Fa_subset(const std::string &basis) const {
    return matrix_subset_helper(Fa_, Ca_, basis, "Fock");
}

SharedMatrix Wavefunction::Fb_subset(const std::string &basis) const {
    return matrix_subset_helper(Fb_, Cb_, basis, "Fock");
}

SharedVector Wavefunction::epsilon_a_subset(const std::string &basis, const std::string &subset) const {
    return epsilon_subset_helper(epsilon_a_, nalphapi_, basis, subset);
}

SharedVector Wavefunction::epsilon_b_subset(const std::string &basis, const std::string &subset) const {
    return epsilon_subset_helper(epsilon_b_, nbetapi_, basis, subset);
}

SharedMatrix Wavefunction::Fa() const { return Fa_; }

SharedMatrix Wavefunction::Fb() const { return Fb_; }

SharedMatrix Wavefunction::Lagrangian() const { return Lagrangian_; }

SharedVector Wavefunction::epsilon_a() const { return epsilon_a_; }

SharedVector Wavefunction::epsilon_b() const { return epsilon_b_; }

const SharedMatrix Wavefunction::Da() const { return Da_; }

SharedMatrix Wavefunction::Db() const { return Db_; }

SharedMatrix Wavefunction::X() const { return Lagrangian_; }

void Wavefunction::set_energy(double ene) {
    set_scalar_variable("CURRENT ENERGY", ene);
}

SharedMatrix Wavefunction::gradient() const { return gradient_; }

void Wavefunction::set_gradient(SharedMatrix grad) {
    set_array_variable("CURRENT GRADIENT", grad);
}

SharedMatrix Wavefunction::hessian() const { return hessian_; }

void Wavefunction::set_hessian(SharedMatrix hess) {
    set_array_variable("CURRENT HESSIAN", hess);
}

SharedVector Wavefunction::frequencies() const { return frequencies_; }

void Wavefunction::set_frequencies(std::shared_ptr<Vector> freqs) { frequencies_ = freqs; }

void Wavefunction::save() const {}

std::shared_ptr<ExternalPotential> Wavefunction::external_pot() const { return external_pot_; }

std::shared_ptr<Vector> Wavefunction::get_esp_at_nuclei() const {
    std::shared_ptr<std::vector<double>> v = esp_at_nuclei();

    int n = molecule_->natom();
    std::shared_ptr<Vector> v_vector(new Vector(n));
    for (int i = 0; i < n; ++i) v_vector->set(i, (*v)[i]);
    return v_vector;
}

std::vector<SharedVector> Wavefunction::get_mo_extents() const {
    std::vector<SharedVector> m = mo_extents();

    int n = nmo_;
    std::vector<SharedVector> mo_vectors;
    mo_vectors.push_back(SharedVector(new Vector("<x^2>", basisset_->nbf())));
    mo_vectors.push_back(SharedVector(new Vector("<y^2>", basisset_->nbf())));
    mo_vectors.push_back(SharedVector(new Vector("<z^2>", basisset_->nbf())));
    mo_vectors.push_back(SharedVector(new Vector("<r^2>", basisset_->nbf())));
    for (int i = 0; i < n; i++) {
        mo_vectors[0]->set(0, i, m[0]->get(0, i));
        mo_vectors[1]->set(0, i, m[1]->get(0, i));
        mo_vectors[2]->set(0, i, m[2]->get(0, i));
        mo_vectors[3]->set(0, i, m[3]->get(0, i));
    }

    return mo_vectors;
}

std::shared_ptr<Vector> Wavefunction::get_atomic_point_charges() const {
    std::shared_ptr<std::vector<double>> q = atomic_point_charges();

    int n = molecule_->natom();
    auto q_vector = std::make_shared<Vector>(n);
    for (int i = 0; i < n; ++i) {
        q_vector->set(i, (*q)[i]);
    }
    return q_vector;
}

std::vector<std::vector<std::tuple<double, int, int>>> Wavefunction::get_no_occupations() const {
    std::vector<std::vector<std::tuple<double, int, int>>> nos = no_occupations();
    int nfsym = nos.size();
    std::vector<std::vector<std::tuple<double, int, int>>> no_occs;
    if (nfsym == 3) {
        no_occs.push_back(nos[0]);
        no_occs.push_back(nos[1]);
        no_occs.push_back(nos[2]);
    } else {
        no_occs.push_back(nos[0]);
    }

    return no_occs;
}

bool Wavefunction::has_scalar_variable(const std::string &key) { return variables_.count(to_upper_copy(key)); }

bool Wavefunction::has_array_variable(const std::string &key) { return arrays_.count(to_upper_copy(key)); }

double Wavefunction::scalar_variable(const std::string &key) {
    std::string uc_key = to_upper_copy(key);

    auto search = variables_.find(uc_key);
    if (search != variables_.end()) {
        return search->second;
    } else {
        throw PSIEXCEPTION("Wavefunction::scalar_variable: Requested variable " + uc_key + " was not set!\n");
    }
}

SharedMatrix Wavefunction::array_variable(const std::string &key) {
    std::string uc_key = to_upper_copy(key);

    auto search = arrays_.find(uc_key);
    if (search != arrays_.end()) {
        return search->second->clone();
    } else {
        throw PSIEXCEPTION("Wavefunction::array_variable: Requested variable " + uc_key + " was not set!\n");
    }
}

void Wavefunction::set_scalar_variable(const std::string &key, double val) {
    variables_[to_upper_copy(key)] = val;

    if (to_upper_copy(key) == "CURRENT ENERGY") energy_ = val;
}

void Wavefunction::set_array_variable(const std::string &key, SharedMatrix val) {
    arrays_[to_upper_copy(key)] = val->clone();

    if (to_upper_copy(key) == "CURRENT GRADIENT") gradient_ = val->clone();
    if (to_upper_copy(key) == "CURRENT HESSIAN") hessian_ = val->clone();
}

int Wavefunction::del_scalar_variable(const std::string &key) { return variables_.erase(to_upper_copy(key)); }

int Wavefunction::del_array_variable(const std::string &key) { return arrays_.erase(to_upper_copy(key)); }

std::map<std::string, double> Wavefunction::scalar_variables() { return variables_; }

std::map<std::string, SharedMatrix> Wavefunction::array_variables() { return arrays_; }

double Wavefunction::get_variable(const std::string &key) { return scalar_variable(key); }
SharedMatrix Wavefunction::get_array(const std::string &key) { return array_variable(key); }
void Wavefunction::set_variable(const std::string &key, double val) { return set_scalar_variable(key, val); }
void Wavefunction::set_array(const std::string &key, SharedMatrix val) { set_array_variable(key, val); }
std::map<std::string, double> Wavefunction::variables() { return scalar_variables(); }
std::map<std::string, SharedMatrix> Wavefunction::arrays() { return array_variables(); }

void Wavefunction::set_PCM(const std::shared_ptr<PCM> &pcm) {
    PCM_ = pcm;
    PCM_enabled_ = true;
}

std::shared_ptr<PCM> Wavefunction::get_PCM() const { return PCM_; }
