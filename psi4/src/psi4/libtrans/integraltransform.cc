/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#include "integraltransform.h"
#include "mospace.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

using namespace psi;

IntegralTransform::IntegralTransform(std::shared_ptr<Wavefunction> wfn, SpaceVec spaces,
                                     TransformationType transformationType, OutputType outputType,
                                     MOOrdering moOrdering, FrozenOrbitals frozenOrbitals, bool init)
    : initialized_(false),
      psio_(_default_psio_lib_),
      sobasis_(wfn->sobasisset()),
      transformationType_(transformationType),
      uniqueSpaces_(spaces),
      moOrdering_(moOrdering),
      outputType_(outputType),
      frozenOrbitals_(frozenOrbitals),
      alreadyPresorted_(false),
      dpdIntFile_(PSIF_LIBTRANS_DPD),
      aHtIntFile_(PSIF_LIBTRANS_A_HT),
      bHtIntFile_(PSIF_LIBTRANS_B_HT),
      nTriSo_(0),
      nTriMo_(0),
      nfzc_(0),
      nfzv_(0),
      spaces_(nullptr),
      labels_(0),
      tolerance_(1.0E-16),
      moIntFileAA_(0),
      moIntFileAB_(0),
      moIntFileBB_(0),
      myDPDNum_(1),
      print_(1),
      zeros_(nullptr),
      sosym_(nullptr),
      mosym_(nullptr),
      aQT_(nullptr),
      bQT_(nullptr),
      cacheFiles_(nullptr),
      cacheList_(nullptr),
      Ca_(wfn->Ca()),
      Cb_(wfn->Cb()),
      H_(wfn->H()),
      keepIwlSoInts_(false),
      keepIwlMoTpdm_(true),
      keepDpdSoInts_(false),
      keepDpdMoTpdm_(true),
      keepHtInts_(true),
      keepHtTpdm_(true),
      tpdmAlreadyPresorted_(false),
      soIntTEIFile_(PSIF_SO_TEI) {
    // Implement set/get functions to customize any of this stuff.  Delayed initialization
    // is possible in case any of these variables need to be changed before setup.
    memory_ = Process::environment.get_memory();

    labels_ = wfn->molecule()->irrep_labels();
    nirreps_ = wfn->nirrep();
    nmo_ = wfn->nmo();
    nso_ = wfn->nso();
    sopi_ = wfn->nsopi();
    mopi_ = wfn->nmopi();
    clsdpi_ = wfn->doccpi();
    openpi_ = wfn->soccpi();
    frzcpi_ = wfn->frzcpi();
    frzvpi_ = wfn->frzvpi();
    nalphapi_ = wfn->nalphapi();
    nbetapi_ = wfn->nbetapi();
    frozen_core_energy_ = 0.0;

    common_initialize();

    if (init) initialize();
}

IntegralTransform::IntegralTransform(SharedMatrix H, SharedMatrix c, SharedMatrix i, SharedMatrix a, SharedMatrix v,
                                     SpaceVec spaces, TransformationType transformationType, OutputType outputType,
                                     MOOrdering moOrdering, FrozenOrbitals frozenOrbitals, bool init)
    : initialized_(false),
      psio_(_default_psio_lib_),
      transformationType_(transformationType),
      uniqueSpaces_(spaces),
      moOrdering_(moOrdering),
      outputType_(outputType),
      frozenOrbitals_(frozenOrbitals),
      alreadyPresorted_(false),
      dpdIntFile_(PSIF_LIBTRANS_DPD),
      aHtIntFile_(PSIF_LIBTRANS_A_HT),
      bHtIntFile_(PSIF_LIBTRANS_B_HT),
      nirreps_(0),
      nmo_(0),
      nso_(0),
      nTriSo_(0),
      nTriMo_(0),
      nfzc_(0),
      nfzv_(0),
      spaces_(nullptr),
      labels_(0),
      tolerance_(1.0E-16),
      memory_(250 * 1024 * 1024),
      moIntFileAA_(0),
      moIntFileAB_(0),
      moIntFileBB_(0),
      myDPDNum_(1),
      print_(1),
      zeros_(nullptr),
      sopi_(0),
      sosym_(nullptr),
      mosym_(nullptr),
      mopi_(0),
      clsdpi_(0),
      openpi_(0),
      frzcpi_(0),
      frzvpi_(0),
      aQT_(nullptr),
      bQT_(nullptr),
      cacheFiles_(nullptr),
      cacheList_(nullptr),
      aCorrToPitzer_(nullptr),
      bCorrToPitzer_(nullptr),
      keepIwlSoInts_(false),
      keepIwlMoTpdm_(true),
      keepDpdSoInts_(false),
      keepDpdMoTpdm_(true),
      keepHtInts_(true),
      keepHtTpdm_(true),
      tpdmAlreadyPresorted_(false),
      soIntTEIFile_(PSIF_SO_TEI) {
    memory_ = Process::environment.get_memory();

    nirreps_ = c->nirrep();
    nmo_ = c->ncol() + i->ncol() + a->ncol() + v->ncol();
    nso_ = i->nrow();
    sopi_ = i->rowspi();  // use i for this since there will always be occupied orbitals
    mopi_ = c->colspi() + i->colspi() + a->colspi() + v->colspi();
    clsdpi_ = i->colspi() + c->colspi();
    openpi_ = Dimension(nirreps_);  // This is the restricted constructor, there are no unpaired electrons
    frzcpi_ = c->colspi();
    frzvpi_ = v->colspi();
    nalphapi_ = clsdpi_ + frzcpi_;  // Restricted constructor
    nbetapi_ = clsdpi_ + frzcpi_;
    frozen_core_energy_ = 0.0;

    // Need to smash together the C's only for them to be ripped apart elsewhere.
    std::vector<SharedMatrix> Cs;
    Cs.push_back(c);
    Cs.push_back(i);
    Cs.push_back(a);
    Cs.push_back(v);
    Ca_ = linalg::horzcat(Cs);
    Cb_ = Ca_;
    H_ = H;

    common_initialize();

    if (init) initialize();
}

/**
 * Sets up the DPD buffers and performs semicanonicalization, if necessary.
 */
void IntegralTransform::initialize() {
    print_ = Process::environment.options.get_int("PRINT");
    printTei_ = print_ > 5;
    useIWL_ = outputType_ == OutputType::IWLAndDPD || outputType_ == OutputType::IWLOnly;
    useDPD_ = outputType_ == OutputType::IWLAndDPD || outputType_ == OutputType::DPDOnly;
    iwlAAIntFile_ = transformationType_ == TransformationType::Restricted ? PSIF_MO_TEI : PSIF_MO_AA_TEI;
    iwlABIntFile_ = transformationType_ == TransformationType::Restricted ? PSIF_MO_TEI : PSIF_MO_AB_TEI;
    iwlBBIntFile_ = transformationType_ == TransformationType::Restricted ? PSIF_MO_TEI : PSIF_MO_BB_TEI;

    tpdm_buffer_ = nullptr;

    aQT_ = init_int_array(nmo_);
    if (transformationType_ == TransformationType::Restricted) {
        reorder_qt(clsdpi_, openpi_, frzcpi_, frzvpi_, aQT_, mopi_, nirreps_);
        bQT_ = aQT_;
    } else {
        bQT_ = init_int_array(nmo_);
        reorder_qt_uhf(clsdpi_, openpi_, frzcpi_, frzvpi_, aQT_, bQT_, mopi_, nirreps_);
    }
    // Set up the correlated to Pitzer arrays.  These have to include the occupied core terms, because
    // the reference contributions are already folded into the TPDM.  However, they don't include frozen
    // virtuals
    aCorrToPitzer_ = init_int_array(nmo_);
    if (transformationType_ != TransformationType::Restricted) {
        bCorrToPitzer_ = init_int_array(nmo_);
    } else {
        bCorrToPitzer_ = aCorrToPitzer_;
    }

    int nFzvFound = 0;
    int pitzerCount = 0;
    for (int h = 0; h < nirreps_; ++h) {
        for (int p = 0; p < mopi_[h]; p++) {
            if (p < mopi_[h] - frzvpi_[h]) {
                // This is active, count it
                int q = aQT_[pitzerCount];
                aCorrToPitzer_[q] = pitzerCount - nFzvFound;
                if (transformationType_ != TransformationType::Restricted) {
                    int q = bQT_[pitzerCount];
                    bCorrToPitzer_[q] = pitzerCount - nFzvFound;
                }
            } else {
                nFzvFound++;
            }
            pitzerCount++;
        }
    }

    if (print_ > 4) {
        outfile->Printf("\tThe Alpha Pitzer to QT mapping array:\n\t\t");
        for (int p = 0; p < nmo_; ++p) {
            outfile->Printf("%d ", aQT_[p]);
        }
        outfile->Printf("\n");
        outfile->Printf("\tThe Beta Pitzer to QT mapping array:\n\t\t");
        for (int p = 0; p < nmo_; ++p) {
            outfile->Printf("%d ", bQT_[p]);
        }
        outfile->Printf("\n");
        outfile->Printf("\tThe Alpha Correlated to Pitzer mapping array:\n\t\t");
        for (int p = 0; p < nmo_; ++p) {
            outfile->Printf("%d ", aCorrToPitzer_[p]);
        }
        outfile->Printf("\n");
        outfile->Printf("\tThe Beta Correlated to Pitzer mapping array:\n\t\t");
        for (int p = 0; p < nmo_; ++p) {
            outfile->Printf("%d ", bCorrToPitzer_[p]);
        }
        outfile->Printf("\n");
    }

    process_spaces();

    // Set up the DPD library
    // TODO implement caching of files
    int numSpaces = spacesUsed_.size();
    int numIndexArrays = numSpaces * (numSpaces - 1) + 5 * numSpaces;
    cacheFiles_ = init_int_array(PSIO_MAXUNIT);
    cacheList_ = init_int_matrix(numIndexArrays, numIndexArrays);
    int currentActiveDPD = psi::dpd_default;
    dpd_init(myDPDNum_, nirreps_, memory_, 0, cacheFiles_, cacheList_, nullptr, numSpaces, spaceArray_);

    // We have to redefine the MO coefficients for a UHF-like treatment
    if (transformationType_ == TransformationType::SemiCanonical) {
        throw PSIEXCEPTION(
            "Semicanonical is deprecated in Libtrans. Please pre-semicanonicalize before passing to libtrans.");
    }
    process_eigenvectors();

    // Return DPD control to the user
    dpd_set_default(currentActiveDPD);

    initialized_ = true;
}

std::shared_ptr<PSIO> IntegralTransform::get_psio() const { return psio_; }

void IntegralTransform::set_psio(std::shared_ptr<PSIO> psio) { psio_ = psio; }

IntegralTransform::~IntegralTransform() {
    if (initialized_) {
        dpd_close(myDPDNum_);
        free_int_matrix(cacheList_);
        free(cacheFiles_);
        free(zeros_);
        free(aQT_);
        free(aCorrToPitzer_);
        if (transformationType_ != TransformationType::Restricted) {
            free(bQT_);
            free(bCorrToPitzer_);
        }
    }
    if (tpdm_buffer_) delete[] tpdm_buffer_;
}

void IntegralTransform::check_initialized() {
    if (initialized_ == false)
        throw PSIEXCEPTION("IntegralTransform::check_initialized: This instance is not initialized.");
}
