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

#include "sapt.h"
#include "psi4/lib3index/3index.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/potential.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libqt/qt.h"

#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {
namespace sapt {

SAPT::SAPT(SharedWavefunction Dimer, SharedWavefunction MonomerA, SharedWavefunction MonomerB, Options &options,
           std::shared_ptr<PSIO> psio)
    : Wavefunction(options) {
    shallow_copy(Dimer);

    if ((Dimer->nirrep() != 1) || (MonomerA->nirrep() != 1) || (MonomerA->nirrep() != 1)) {
        throw PSIEXCEPTION("SAPT must be run in C1 symmetry. Period.");
    }

    if ((Dimer->soccpi().sum() != 0) || (MonomerA->soccpi().sum() != 0) || (MonomerA->soccpi().sum() != 0)) {
        throw PSIEXCEPTION("This is a RHF SAPT constructor. Pair those electrons up cracker!");
    }

    psio_ = psio;

#ifdef USING_LAPACK_MKL
    mkl_set_dynamic(1);
#endif

#ifdef _OPENMP
    omp_set_max_active_levels(1);
#endif

    initialize(MonomerA, MonomerB);
    get_denom();
}

SAPT::~SAPT() {
    if (evalsA_ != nullptr) free(evalsA_);
    if (evalsB_ != nullptr) free(evalsB_);
    if (diagAA_ != nullptr) free(diagAA_);
    if (diagBB_ != nullptr) free(diagBB_);
    if (CA_ != nullptr) free_block(CA_);
    if (CB_ != nullptr) free_block(CB_);
    if (CHFA_ != nullptr) free_block(CHFA_);
    if (CHFB_ != nullptr) free_block(CHFB_);
    if (sAB_ != nullptr) free_block(sAB_);
    if (vABB_ != nullptr) free_block(vABB_);
    if (vBAA_ != nullptr) free_block(vBAA_);
    if (vAAB_ != nullptr) free_block(vAAB_);
    if (vBAB_ != nullptr) free_block(vBAB_);
    zero_.reset();
}

void SAPT::initialize(SharedWavefunction MonomerA, SharedWavefunction MonomerB) {
    evalsA_ = nullptr;
    evalsB_ = nullptr;
    diagAA_ = nullptr;
    diagBB_ = nullptr;
    CA_ = nullptr;
    CB_ = nullptr;
    CHFA_ = nullptr;
    CHFB_ = nullptr;
    sAB_ = nullptr;
    vABB_ = nullptr;
    vBAA_ = nullptr;
    vAAB_ = nullptr;
    vBAB_ = nullptr;

    // We inherit from the dimer basis
    ribasis_ = get_basisset("DF_BASIS_SAPT");
    elstbasis_ = get_basisset("DF_BASIS_ELST");

    // Compare pointers
    if (ribasis_ == elstbasis_) {
        elst_basis_ = false;
    } else {
        elst_basis_ = true;
    }

    zero_ = std::shared_ptr<BasisSet>(BasisSet::zero_ao_basis_set());

    if (options_.get_str("EXCH_SCALE_ALPHA") == "FALSE") {
        exch_scale_alpha_ = 0.0;
    } else if (options_.get_str("EXCH_SCALE_ALPHA") == "TRUE") {
        exch_scale_alpha_ = 1.0;  // Default value for alpha
    } else {
        exch_scale_alpha_ = std::atof(options_.get_str("EXCH_SCALE_ALPHA").c_str());
    }
    Process::environment.globals["SAPT ALPHA"] = exch_scale_alpha_;
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");
    schwarz_ = options_.get_double("INTS_TOLERANCE");
    mem_ = (long int)((double)memory_ * options_.get_double("SAPT_MEM_SAFETY"));
    mem_ /= 8L;

    std::vector<int> realsA;
    realsA.push_back(0);
    std::vector<int> ghostsA;
    ghostsA.push_back(1);
    auto monomerA = std::shared_ptr<Molecule>(molecule_->extract_subsets(realsA, ghostsA));
    foccA_ = MonomerA->basisset()->n_frozen_core(options_.get_str("FREEZE_CORE"), monomerA);

    std::vector<int> realsB;
    realsB.push_back(1);
    std::vector<int> ghostsB;
    ghostsB.push_back(0);
    auto monomerB = std::shared_ptr<Molecule>(molecule_->extract_subsets(realsB, ghostsB));
    foccB_ = MonomerB->basisset()->n_frozen_core(options_.get_str("FREEZE_CORE"), monomerB);

    natomsA_ = 0;
    natomsB_ = 0;

    for (int n = 0; n < monomerA->natom(); n++)
        if (monomerA->Z(n)) natomsA_++;
    for (int n = 0; n < monomerB->natom(); n++)
        if (monomerB->Z(n)) natomsB_++;

    ndf_ = ribasis_->nbf();

    double enucD, enucA, enucB;
    double eHFD, eHFA, eHFB;

    eHFD = energy_;
    enucD = molecule_->nuclear_repulsion_energy(dipole_field_strength_);

    // Monomer A info
    nsoA_ = MonomerA->nso();
    nmoA_ = MonomerA->nmo();
    noccA_ = MonomerA->doccpi().sum();
    nvirA_ = nmoA_ - noccA_;
    NA_ = 2 * noccA_;
    eHFA = MonomerA->energy();
    enucA = MonomerA->molecule()->nuclear_repulsion_energy(dipole_field_strength_);
    aoccA_ = noccA_ - foccA_;

    // Monomer B info
    nsoB_ = MonomerB->nso();
    nmoB_ = MonomerB->nmo();
    noccB_ = MonomerB->doccpi().sum();
    nvirB_ = nmoB_ - noccB_;
    NB_ = 2 * noccB_;
    eHFB = MonomerB->energy();
    enucB = MonomerB->molecule()->nuclear_repulsion_energy(dipole_field_strength_);
    aoccB_ = noccB_ - foccB_;

    enuc_ = enucD - enucA - enucB;
    eHF_ = eHFD - eHFA - eHFB;

    evalsA_ = init_array(nmoA_);
    std::memcpy(evalsA_, MonomerA->epsilon_a()->pointer(), sizeof(double) * nmoA_);

    evalsB_ = init_array(nmoB_);
    std::memcpy(evalsB_, MonomerB->epsilon_a()->pointer(), sizeof(double) * nmoB_);

    CA_ = block_matrix(nso_, nmoA_);
    double **tempA = block_matrix(nsoA_, nmoA_);
    std::memcpy(tempA[0], MonomerA->Ca()->pointer()[0], sizeof(double) * nmoA_ * nsoA_);

    if (nsoA_ != nso_) {
        for (int n = 0; n < nsoA_; n++) C_DCOPY(nmoA_, tempA[n], 1, CA_[n], 1);
    } else
        C_DCOPY(nso_ * nmoA_, tempA[0], 1, CA_[0], 1);
    free_block(tempA);

    CB_ = block_matrix(nso_, nmoB_);
    double **tempB = block_matrix(nsoB_, nmoB_);
    std::memcpy(tempB[0], MonomerB->Ca()->pointer()[0], sizeof(double) * nmoB_ * nsoB_);
    if (nsoB_ != nso_) {
        for (int n = 0; n < nsoB_; n++) C_DCOPY(nmoB_, tempB[n], 1, CB_[n + nsoA_], 1);
    } else
        C_DCOPY(nso_ * nmoB_, tempB[0], 1, CB_[0], 1);
    free_block(tempB);

    int nbf[8];
    nbf[0] = nso_;
    auto fact = std::make_shared<MatrixFactory>();
    fact->init_with(1, nbf, nbf);

    auto intfact = std::make_shared<IntegralFactory>(basisset_, basisset_, basisset_, basisset_);

    std::shared_ptr<OneBodyAOInt> Sint(intfact->ao_overlap());
    Smat_ = std::make_shared<Matrix>(fact->create_matrix("Overlap"));
    Sint->compute(Smat_);

    double **sIJ = Smat_->pointer();
    double **sAJ = block_matrix(nmoA_, nso_);
    sAB_ = block_matrix(nmoA_, nmoB_);

    C_DGEMM('T', 'N', nmoA_, nso_, nso_, 1.0, CA_[0], nmoA_, sIJ[0], nso_, 0.0, sAJ[0], nso_);
    C_DGEMM('N', 'N', nmoA_, nmoB_, nso_, 1.0, sAJ[0], nso_, CB_[0], nmoB_, 0.0, sAB_[0], nmoB_);

    free_block(sAJ);

    auto potA = std::shared_ptr<PotentialInt>(dynamic_cast<PotentialInt *>(intfact->ao_potential().release()));
    std::vector<std::pair<double, std::array<double, 3>>> ZxyzA;
    for (int n = 0, p = 0; n < monomerA->natom(); n++) {
        if (monomerA->Z(n)) {
            ZxyzA.push_back({(double)monomerA->Z(n), {{monomerA->x(n), monomerA->y(n), monomerA->z(n)}}});
        }
    }
    potA->set_charge_field(ZxyzA);
    VAmat_ = std::make_shared<Matrix>(fact->create_matrix("Nuclear Attraction (Monomer A)"));
    potA->compute(VAmat_);

    auto potB = std::shared_ptr<PotentialInt>(dynamic_cast<PotentialInt *>(intfact->ao_potential().release()));
    std::vector<std::pair<double, std::array<double, 3>>> ZxyzB;
    for (int n = 0, p = 0; n < monomerB->natom(); n++) {
        if (monomerB->Z(n)) {
            ZxyzB.push_back({(double)monomerB->Z(n), {{monomerB->x(n), monomerB->y(n), monomerB->z(n)}}});
        }
    }
    potB->set_charge_field(ZxyzB);
    VBmat_ = std::make_shared<Matrix>(fact->create_matrix("Nuclear Attraction (Monomer B)"));
    potB->compute(VBmat_);

    double **vIB = block_matrix(nso_, nmoB_);
    double **vAJ = block_matrix(nmoA_, nso_);
    vAAB_ = block_matrix(nmoA_, nmoB_);
    vABB_ = block_matrix(nmoB_, nmoB_);
    vBAA_ = block_matrix(nmoA_, nmoA_);
    vBAB_ = block_matrix(nmoA_, nmoB_);

    double **vIJ = VAmat_->pointer();

    C_DGEMM('N', 'N', nso_, nmoB_, nso_, 1.0, vIJ[0], nso_, CB_[0], nmoB_, 0.0, vIB[0], nmoB_);
    C_DGEMM('T', 'N', nmoA_, nmoB_, nso_, 1.0, CA_[0], nmoA_, vIB[0], nmoB_, 0.0, vAAB_[0], nmoB_);
    C_DGEMM('T', 'N', nmoB_, nmoB_, nso_, 1.0, CB_[0], nmoB_, vIB[0], nmoB_, 0.0, vABB_[0], nmoB_);

    vIJ = VBmat_->pointer();

    C_DGEMM('T', 'N', nmoA_, nso_, nso_, 1.0, CA_[0], nmoA_, vIJ[0], nso_, 0.0, vAJ[0], nso_);
    C_DGEMM('N', 'N', nmoA_, nmoA_, nso_, 1.0, vAJ[0], nso_, CA_[0], nmoA_, 0.0, vBAA_[0], nmoA_);
    C_DGEMM('N', 'N', nmoA_, nmoB_, nso_, 1.0, vAJ[0], nso_, CB_[0], nmoB_, 0.0, vBAB_[0], nmoB_);

    free_block(vIB);
    free_block(vAJ);

//Konrad: extra stuff for nonapproximate E(30)ex-ind in AOs
    CoccA_ = MonomerA->Ca_subset("AO", "OCC");
    CoccB_ = MonomerB->Ca_subset("AO", "OCC");
    CvirA_ = MonomerA->Ca_subset("AO", "VIR");
    CvirB_ = MonomerB->Ca_subset("AO", "VIR");
}

SharedMatrix SAPT::get_metric(std::shared_ptr<BasisSet> basis) const {
    FittingMetric metric(basis);
    metric.form_eig_inverse(options_.get_double("DF_FITTING_CONDITION"));
    return metric.get_metric();
}

void SAPT::get_denom() {
    auto evals_aoccA = std::make_shared<Vector>(aoccA_);
    auto evals_virA = std::make_shared<Vector>(nvirA_);
    auto evals_aoccB = std::make_shared<Vector>(aoccB_);
    auto evals_virB = std::make_shared<Vector>(nvirB_);

    for (int a = 0; a < aoccA_; a++) evals_aoccA->set(0, a, evalsA_[a + foccA_]);
    for (int r = 0; r < nvirA_; r++) evals_virA->set(0, r, evalsA_[r + noccA_]);
    for (int b = 0; b < aoccB_; b++) evals_aoccB->set(0, b, evalsB_[b + foccB_]);
    for (int s = 0; s < nvirB_; s++) evals_virB->set(0, s, evalsB_[s + noccB_]);

    denom_ =
        SAPTDenominator::buildDenominator(options_.get_str("DENOMINATOR_ALGORITHM"), evals_aoccA, evals_virA,
                                          evals_aoccB, evals_virB, options_.get_double("DENOMINATOR_DELTA"), debug_);

    if (debug_ > 1) denom_->debug();

    SharedMatrix tauAR = denom_->denominatorA();
    SharedMatrix tauBS = denom_->denominatorB();

    dAR_ = tauAR->pointer();
    dBS_ = tauBS->pointer();

    nvec_ = denom_->nvector();
}
}  // namespace sapt
}  // namespace psi
