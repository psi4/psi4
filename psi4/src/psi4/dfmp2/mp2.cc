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

#include "mp2.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/psi4-dec.h"
#include "psi4/physconst.h"
#include "psi4/psifiles.h"

#include "psi4/lib3index/3index.h"
#include "psi4/libfock/apps.h"
#include "psi4/libfock/jk.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/extern.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/oeprop.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"

#include "corr_grad.h"

namespace psi {
namespace dfmp2 {

void DFMP2::compute_opdm_and_nos(const SharedMatrix Dnosym, SharedMatrix Dso, SharedMatrix Cno, SharedVector occ) {
    // The density matrix
    auto c1MO_c1NO = std::make_shared<Matrix>("NOs", nmo_, nmo_);
    auto occ_c1 = std::make_shared<Vector>("NO Occupations", nmo_);
    Dnosym->diagonalize(c1MO_c1NO, occ_c1, descending);
    // Rotate the canonical MOs to NOs
    auto AO_c1MO = reference_wavefunction_->Ca_subset("AO");
    auto AO_c1NO = AO_c1MO->clone();
    AO_c1NO->gemm(false, false, 1.0, AO_c1MO, c1MO_c1NO, 0.0);
    // Reapply the symmetry to the AO dimension
    auto AO_SO = reference_wavefunction_->aotoso();
    auto SO_c1NO = std::make_shared<Matrix>(nirrep_, (const int*)nsopi_, nmo_);
    SO_c1NO->set_name("SO to C1 NO");
    for (int h = 0; h < nirrep_; ++h) {
        int so_h = nsopi_[h];
        if (so_h) {
            double** pAOSO = AO_SO->pointer(h);
            double** pAONO = AO_c1NO->pointer(0);
            double** pSONO = SO_c1NO->pointer(h);
            C_DGEMM('T', 'N', so_h, nmo_, nso_, 1.0, pAOSO[0], so_h, pAONO[0], nmo_, 0.0, pSONO[0], nmo_);
        }
    }
    // Now, copy over the full matrix, whenever nonzero columns are
    for (int h = 0; h < nirrep_; ++h) {
        if (nsopi_[h] == 0) continue;
        std::vector<double> CStemp(nsopi_[h]);
        auto pC1 = SO_c1NO->pointer(h);
        auto Smat = S_->pointer(h);
        int symcol = 0;
        for (int col = 0; col < nmo_; ++col) {
            // Compute orthonormalized self-overlap, to see if it's nonzero.
            // If it is, grab this orbital and store it in the symmetry NO matrix.
            C_DGEMV('n', nsopi_[h], nsopi_[h], 1.0, Smat[0], nsopi_[h], &(pC1[0][col]), nmo_, 0.0, CStemp.data(), 1);
            double overlap = C_DDOT(nsopi_[h], CStemp.data(), 1, &(pC1[0][col]), nmo_);
            if (overlap > 0.8) {
                for (int row = 0; row < nsopi_[h]; ++row) {
                    Cno->set(h, row, symcol, pC1[row][col]);
                }
                occ->set(h, symcol, occ_c1->get(col));
                symcol++;
            }
        }
        if (symcol != nmopi_[h]) {
            outfile->Printf(
                "Problem determining natural orbital and density matrix symmetries.\n"
                "Future calls to oeprop using this density will not work. Try disabling symmetry.\n\n");
            occ->zero();
            Cno->zero();
            Dso->zero();
            return;
        }
    }

    // Backtransform Density matrix to the SO basis
    //
    // D(SO) = C[SO->MO(c1)] D[MO] SO->MO[C1]t
    //       = C[SO->MO(c1)] U[MO(c1)->NO(c1)] diag(NO occ) U[MO(c1)->NO(c1)]t C[SO->MO(c1)]t
    //       = C[SO->NO(c1)] diag(NO occ) C[SO->NO(c1)]t
    //
    // Now we're finished with the MO(c1)->NO(c1) matrix, use it as scratch for diag(NO occ)
    c1MO_c1NO->set_diagonal(occ_c1);
    auto temp = std::make_shared<Matrix>(nirrep_, (const int*)nsopi_, nmo_);
    for (int h = 0; h < nirrep_; ++h) {
        int so_h = nsopi_[h];
        if (so_h) {
            double** ptemp = temp->pointer(h);
            double** pDso = Dso->pointer(h);
            double** pC = SO_c1NO->pointer(h);
            double** pOcc = c1MO_c1NO->pointer(0);
            // temp = C[SO->MO(c1)] X diag(NO occ)
            C_DGEMM('N', 'N', so_h, nmo_, nmo_, 1.0, pC[0], nmo_, pOcc[0], nmo_, 0.0, ptemp[0], nmo_);
            // Da = C[SO->MO(c1)] X diag(NO occ)
            C_DGEMM('N', 'T', so_h, so_h, nmo_, 1.0, ptemp[0], nmo_, pC[0], nmo_, 0.0, pDso[0], so_h);
        }
    }
}

void DFMP2::block_status(std::vector<int> inds, const char* file, int line) {
    bool gimp = false;
    if (inds.size() > 2) {
        gimp = ((inds[inds.size() - 1] - inds[inds.size() - 2]) != (inds[1] - inds[0]));
    }
    printf("%s:%d %zu %s %d %d\n", file, line, inds.size(), (gimp ? "GIMP" : "NOT GIMP"), inds[1] - inds[0],
           inds[inds.size() - 1] - inds[inds.size() - 2]);
}
void DFMP2::block_status(std::vector<size_t> inds, const char* file, int line) {
    bool gimp = false;
    if (inds.size() > 2) {
        gimp = ((inds[inds.size() - 1] - inds[inds.size() - 2]) != (inds[1] - inds[0]));
    }
    printf("%s:%d %zu %s %zu %zu\n", file, line, inds.size(), (gimp ? "GIMP" : "NOT GIMP"), inds[1] - inds[0],
           inds[inds.size() - 1] - inds[inds.size() - 2]);
}

DFMP2::DFMP2(SharedWavefunction ref_wfn, Options& options, std::shared_ptr<PSIO> psio) : Wavefunction(options) {
    shallow_copy(ref_wfn);
    reference_wavefunction_ = ref_wfn;
    psio_ = psio;

    common_init();
}
DFMP2::~DFMP2() {}
void DFMP2::common_init() {
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");

    // if (options_.get_str("REFERENCE") == "ROHF" || options_.get_str("REFERENCE") == "CUHF")
    //    throw PSIEXCEPTION("SemiCanonical transform does not work at the moment");
    // reference_wavefunction_->semicanonicalize();

    // copy(reference_wavefunction_);
    name_ = "DF-MP2";
    module_ = "dfmp2";

    variables_["MP2 SINGLES ENERGY"] = 0.0;
    variables_["MP2 DOUBLES ENERGY"] = 0.0;
    variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = 0.0;
    variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = 0.0;
    variables_["SCF TOTAL ENERGY"] = reference_wavefunction_->energy();

    sss_ = options_.get_double("MP2_SS_SCALE");
    oss_ = options_.get_double("MP2_OS_SCALE");

    ribasis_ = get_basisset("DF_BASIS_MP2");
}
double DFMP2::compute_energy() {
    print_header();
    auto num_alpha_excit = std::min(Ca_subset("AO", "ACTIVE_OCC")->colspi()[0], Ca_subset("AO", "ACTIVE_VIR")->colspi()[0]);
    auto num_beta_excit = std::min(Cb_subset("AO", "ACTIVE_OCC")->colspi()[0], Cb_subset("AO", "ACTIVE_VIR")->colspi()[0]);
    if (num_alpha_excit + num_beta_excit < 2) {
        outfile->Printf("  Warning: no MP2 calculation performed on system with alpha ACTIVE_OCC (%d), beta ACTIVE_OCC (%d), alpha ACTIVE_VIR (%d), beta ACTIVE_VIR (%d); returning 0.0 by definition.\n", Ca_subset("AO", "ACTIVE_OCC")->colspi()[0], Cb_subset("AO", "ACTIVE_OCC")->colspi()[0], Ca_subset("AO", "ACTIVE_VIR")->colspi()[0], Cb_subset("AO", "ACTIVE_VIR")->colspi()[0]);
        variables_["MP2 SINGLES ENERGY"] = 0.0;
        variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = 0.0;
        variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = 0.0;
        variables_["MP2 CORRELATION ENERGY"] = 0.0;
        print_energies();
        energy_ = variables_["MP2 TOTAL ENERGY"];
        return variables_["MP2 TOTAL ENERGY"];
    }
    timer_on("DFMP2 Singles");
    form_singles();
    timer_off("DFMP2 Singles");
    timer_on("DFMP2 Aia");
    form_Aia();
    timer_off("DFMP2 Aia");
    timer_on("DFMP2 Bia");
    form_Bia();
    timer_off("DFMP2 Bia");
    timer_on("DFMP2 Energy");
    form_energy();
    timer_off("DFMP2 Energy");
    print_energies();
    energy_ = variables_["MP2 TOTAL ENERGY"];

    return variables_["MP2 TOTAL ENERGY"];
}
SharedMatrix DFMP2::compute_gradient() {
    print_header();

    if (Ca_subset("AO", "ACTIVE_OCC")->colspi()[0] == 0) {
        if (Cb_subset("AO", "ACTIVE_OCC")->colspi()[0] == 0) {
            throw PSIEXCEPTION("There are no occupied orbitals with alpha or beta spin.");
        }
        throw PSIEXCEPTION("There are no occupied orbitals with alpha spin.");
    }
    if (Cb_subset("AO", "ACTIVE_OCC")->colspi()[0] == 0) {
        throw PSIEXCEPTION("There are no occupied orbitals with beta spin.");
    }

    if (Ca_subset("AO", "ACTIVE_VIR")->colspi()[0] == 0) {
        if (Cb_subset("AO", "ACTIVE_VIR")->colspi()[0] == 0) {
            throw PSIEXCEPTION("There are no virtual orbitals with alpha or beta spin.");
        }
        throw PSIEXCEPTION("There are no virtual orbitals with alpha spin.");
    }
    if (Cb_subset("AO", "ACTIVE_VIR")->colspi()[0] == 0) {
        throw PSIEXCEPTION("There are no virtual orbitals with beta spin.");
    }

    timer_on("DFMP2 Singles");
    form_singles();
    timer_off("DFMP2 Singles");

    timer_on("DFMP2 Aia");
    form_Aia();
    timer_off("DFMP2 Aia");

    timer_on("DFMP2 iaB");
    form_Bia_Cia();
    timer_off("DFMP2 iaB");

    timer_on("DFMP2 aiB");
    form_Bia_transpose();
    timer_off("DFMP2 aiB");

    timer_on("DFMP2 Tij");
    form_Pab();
    timer_off("DFMP2 Tij");

    timer_on("DFMP2 Tab");
    form_Pij();
    timer_off("DFMP2 Tab");

    timer_on("DFMP2 gamma");
    form_gamma();
    timer_off("DFMP2 gamma");

    timer_on("DFMP2 Gia");
    form_G_transpose();
    timer_off("DFMP2 Gia");

    timer_on("DFMP2 AB^x");
    form_AB_x_terms();
    timer_off("DFMP2 AB^x");

    timer_on("DFMP2 Amn^x");
    form_Amn_x_terms();
    timer_off("DFMP2 Amn^x");

    timer_on("DFMP2 L");
    form_L();
    timer_off("DFMP2 L");

    timer_on("DFMP2 P");
    form_P();
    timer_off("DFMP2 P");

    timer_on("DFMP2 W");
    form_W();
    timer_off("DFMP2 W");

    timer_on("DFMP2 Z");
    form_Z();
    timer_off("DFMP2 Z");

    if (options_.get_bool("ONEPDM")) {
        print_energies();
        energy_ = variables_["MP2 TOTAL ENERGY"];
        return std::make_shared<Matrix>("nullptr", 0, 0);
    }

    timer_on("DFMP2 grad");
    form_gradient();
    timer_off("DFMP2 grad");

    print_energies();
    energy_ = variables_["MP2 TOTAL ENERGY"];

    // Symmetrize the gradient to remove noise
    gradients_["Total"]->symmetrize_gradient(molecule_);

    print_gradients();

    return gradients_["Total"];
}
void DFMP2::form_singles() {
    double E_singles_a = 0.0;
    double E_singles_b = 0.0;

    auto Caocc_a = Ca_subset("SO", "ACTIVE_OCC");
    auto Cavir_a = Ca_subset("SO", "ACTIVE_VIR");
    auto Caocc_b = Cb_subset("SO", "ACTIVE_OCC");
    auto Cavir_b = Cb_subset("SO", "ACTIVE_VIR");

    auto eps_aocc_a = epsilon_a_subset("SO", "ACTIVE_OCC");
    auto eps_avir_a = epsilon_a_subset("SO", "ACTIVE_VIR");
    auto eps_aocc_b = epsilon_b_subset("SO", "ACTIVE_OCC");
    auto eps_avir_b = epsilon_b_subset("SO", "ACTIVE_VIR");

    auto Fia_a = linalg::triplet(Caocc_a, Fa_, Cavir_a, true, false, false);
    auto Fia_b = linalg::triplet(Caocc_b, Fb_, Cavir_b, true, false, false);

    // DEFINITION: E_singles = - f^i_a f^a_i / (ea - ei)
    // Alpha spin
    for (int h = 0; h < Caocc_a->nirrep(); h++) {
        int nso = Fa_->rowspi()[h];
        int naocc = Caocc_a->colspi()[h];
        int navir = Cavir_a->colspi()[h];

        if (!nso || !naocc || !navir) continue;

        auto Fmop = Fia_a->pointer(h);

        auto eps_i = eps_aocc_a->pointer(h);
        auto eps_a = eps_avir_a->pointer(h);

        for (int i = 0; i < naocc; i++) {
            for (int a = 0; a < navir; a++) {
                E_singles_a -= Fmop[i][a] * Fmop[i][a] / (eps_a[a] - eps_i[i]);
            }
        }
    }

    // Beta spin
    for (int h = 0; h < Caocc_b->nirrep(); h++) {
        int nso = Fb_->rowspi()[h];
        int naocc = Caocc_b->colspi()[h];
        int navir = Cavir_b->colspi()[h];

        if (!nso || !naocc || !navir) continue;

        auto Fmop = Fia_b->pointer(h);

        auto eps_i = eps_aocc_b->pointer(h);
        auto eps_a = eps_avir_b->pointer(h);

        for (int i = 0; i < naocc; i++) {
            for (int a = 0; a < navir; a++) {
                E_singles_b -= Fmop[i][a] * Fmop[i][a] / (eps_a[a] - eps_i[i]);
            }
        }
    }

    variables_["MP2 SINGLES ENERGY"] = E_singles_a + E_singles_b;

    if (debug_) {
        Caocc_a->print();
        Cavir_a->print();
        eps_aocc_a->print();
        eps_avir_a->print();
        Caocc_b->print();
        Cavir_b->print();
        eps_aocc_b->print();
        eps_avir_b->print();

        Fia_a->print();
        Fia_b->print();
        outfile->Printf("  Alpha singles energy = %24.16E\n", E_singles_a);
        outfile->Printf("  Beta  singles energy = %24.16E\n\n", E_singles_b);
    }
}
SharedMatrix DFMP2::form_inverse_metric() {
    timer_on("DFMP2 Metric");

    int naux = ribasis_->nbf();

    // Load inverse metric from the SCF three-index integral file if it exists
    if (options_.get_str("DF_INTS_IO") == "LOAD") {
        auto Jm12 = std::make_shared<Matrix>("SO Basis Fitting Inverse (Eig)", naux, naux);
        outfile->Printf("\t Will attempt to load fitting metric from file %d.\n\n", PSIF_DFSCF_BJ);
        psio_->open(PSIF_DFSCF_BJ, PSIO_OPEN_OLD);
        psio_->read_entry(PSIF_DFSCF_BJ, "DFMP2 Jm12", (char*)Jm12->pointer()[0], sizeof(double) * naux * naux);
        psio_->close(PSIF_DFSCF_BJ, 1);

        timer_off("DFMP2 Metric");

        return Jm12;

    } else {
        // Form the inverse metric manually
        auto metric = std::make_shared<FittingMetric>(ribasis_, true);
        metric->form_eig_inverse(options_.get_double("DF_FITTING_CONDITION"));
        auto Jm12 = metric->get_metric();

        // Save inverse metric to the SCF three-index integral file if it exists
        if (options_.get_str("DF_INTS_IO") == "SAVE") {
            outfile->Printf("\t Will save fitting metric to file %d.\n\n", PSIF_DFSCF_BJ);
            psio_->open(PSIF_DFSCF_BJ, PSIO_OPEN_OLD);
            psio_->write_entry(PSIF_DFSCF_BJ, "DFMP2 Jm12", (char*)Jm12->pointer()[0], sizeof(double) * naux * naux);
            psio_->close(PSIF_DFSCF_BJ, 1);
        }

        timer_off("DFMP2 Metric");

        return Jm12;
    }
}
void DFMP2::apply_fitting(SharedMatrix Jm12, size_t file, size_t naux, size_t nia) {
    // Memory constraints
    size_t Jmem = naux * naux;
    size_t doubles = (size_t)(options_.get_double("DFMP2_MEM_FACTOR") * (memory_ / 8L));
    if (doubles < 2L * Jmem) {
        throw PSIEXCEPTION("DFMP2: More memory required for tractable disk transpose");
    }
    size_t rem = (doubles - Jmem) / 2L;
    size_t max_nia = (rem / naux);
    max_nia = (max_nia > nia ? nia : max_nia);
    max_nia = (max_nia < 1L ? 1L : max_nia);

    // Block sizing
    std::vector<size_t> ia_starts;
    ia_starts.push_back(0);
    for (size_t ia = 0L; ia < nia; ia += max_nia) {
        if (ia + max_nia >= nia) {
            ia_starts.push_back(nia);
        } else {
            ia_starts.push_back(ia + max_nia);
        }
    }
    // block_status(ia_starts, __FILE__,__LINE__);

    // Tensor blocks
    auto Aia = std::make_shared<Matrix>("A(Q|ia)", naux, max_nia);
    auto Bia = std::make_shared<Matrix>("B(Q|ia)", max_nia, naux);
    auto Aiap = Aia->pointer();
    auto Biap = Bia->pointer();
    auto Jp = Jm12->pointer();

    // Loop through blocks
    psio_->open(file, PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;
    psio_address next_BIA = PSIO_ZERO;
    for (int block = 0; block < ia_starts.size() - 1; block++) {
        // Sizing
        size_t ia_start = ia_starts[block];
        size_t ia_stop = ia_starts[block + 1];
        size_t ncols = ia_stop - ia_start;

        // Read Aia
        timer_on("DFMP2 Aia Read");
        for (size_t Q = 0; Q < naux; Q++) {
            next_AIA = psio_get_address(PSIO_ZERO, sizeof(double) * (Q * nia + ia_start));
            psio_->read(file, "A(Q|ia)", (char*)Aiap[Q], sizeof(double) * ncols, next_AIA, &next_AIA);
        }
        timer_off("DFMP2 Aia Read");

        // DEFINITION: B(ia|Q) := A(ia|A)(A|Q)
        timer_on("DFMP2 (ia|A)(A|Q)");
        C_DGEMM('T', 'N', ncols, naux, naux, 1.0, Aiap[0], max_nia, Jp[0], naux, 0.0, Biap[0], naux);
        timer_off("DFMP2 (ia|A)(A|Q)");

        timer_on("DFMP2 Bia Write");
        psio_->write(file, "B(ia|Q)", (char*)Biap[0], sizeof(double) * ncols * naux, next_BIA, &next_BIA);
        timer_off("DFMP2 Bia Write");
    }
    psio_->close(file, 1);
}
void DFMP2::apply_fitting_grad(SharedMatrix Jm12, size_t file, size_t naux, size_t nia) {
    // Memory constraints
    size_t Jmem = naux * naux;
    size_t doubles = (size_t)(options_.get_double("DFMP2_MEM_FACTOR") * (memory_ / 8L));
    if (doubles < 2L * Jmem) {
        throw PSIEXCEPTION("DFMP2: More memory required for tractable disk transpose");
    }
    size_t rem = (doubles - Jmem) / 2L;
    size_t max_nia = (rem / naux);
    max_nia = (max_nia > nia ? nia : max_nia);
    max_nia = (max_nia < 1L ? 1L : max_nia);

    // Block sizing
    std::vector<size_t> ia_starts;
    ia_starts.push_back(0);
    for (size_t ia = 0L; ia < nia; ia += max_nia) {
        if (ia + max_nia >= nia) {
            ia_starts.push_back(nia);
        } else {
            ia_starts.push_back(ia + max_nia);
        }
    }
    // block_status(ia_starts, __FILE__,__LINE__);

    // Tensor blocks
    auto Bia = std::make_shared<Matrix>("B(ia|Q)", max_nia, naux);
    auto Cia = std::make_shared<Matrix>("C(ia|Q)", max_nia, naux);
    auto Biap = Bia->pointer();
    auto Ciap = Cia->pointer();
    auto Jp = Jm12->pointer();

    // Loop through blocks
    psio_->open(file, PSIO_OPEN_OLD);
    psio_address next_BIA = PSIO_ZERO;
    psio_address next_CIA = PSIO_ZERO;
    for (int block = 0; block < ia_starts.size() - 1; block++) {
        // Sizing
        size_t ia_start = ia_starts[block];
        size_t ia_stop = ia_starts[block + 1];
        size_t ncols = ia_stop - ia_start;

        timer_on("DFMP2 Bia Read");
        psio_->read(file, "B(ia|Q)", (char*)Biap[0], sizeof(double) * ncols * naux, next_BIA, &next_BIA);
        timer_off("DFMP2 Bia Read");

        // DEFINITION: C(ia|Q) := B(ia|B)(C|Q)
        timer_on("DFMP2 (ia|B)(B|Q)");
        C_DGEMM('N', 'N', ncols, naux, naux, 1.0, Biap[0], naux, Jp[0], naux, 0.0, Ciap[0], naux);
        timer_off("DFMP2 (ia|B)(B|Q)");

        timer_on("DFMP2 Cia Write");
        psio_->write(file, "C(ia|Q)", (char*)Ciap[0], sizeof(double) * ncols * naux, next_CIA, &next_CIA);
        timer_off("DFMP2 Cia Write");
    }
    psio_->close(file, 1);
}
void DFMP2::apply_gamma(size_t file, size_t naux, size_t nia) {
    size_t Jmem = naux * naux;
    size_t doubles = (size_t)(options_.get_double("DFMP2_MEM_FACTOR") * (memory_ / 8L));
    if (doubles < 1L * Jmem) {
        throw PSIEXCEPTION("DFMP2: More memory required for gamma");
    }
    size_t rem = (doubles - Jmem) / 2L;
    size_t max_nia = (rem / naux);
    max_nia = (max_nia > nia ? nia : max_nia);
    max_nia = (max_nia < 1L ? 1L : max_nia);

    // Block sizing
    std::vector<size_t> ia_starts;
    ia_starts.push_back(0);
    for (size_t ia = 0L; ia < nia; ia += max_nia) {
        if (ia + max_nia >= nia) {
            ia_starts.push_back(nia);
        } else {
            ia_starts.push_back(ia + max_nia);
        }
    }
    // block_status(ia_starts, __FILE__,__LINE__);

    // Tensor blocks
    auto Gia = std::make_shared<Matrix>("G(ia|Q)", max_nia, naux);
    auto Cia = std::make_shared<Matrix>("C(ia|Q)", max_nia, naux);
    auto G = std::make_shared<Matrix>("G_PQ", naux, naux);
    double** Giap = Gia->pointer();
    double** Ciap = Cia->pointer();
    double** Gp = G->pointer();

    // Loop through blocks
    psio_->open(file, PSIO_OPEN_OLD);
    psio_address next_GIA = PSIO_ZERO;
    psio_address next_CIA = PSIO_ZERO;
    for (int block = 0; block < ia_starts.size() - 1; block++) {
        // Sizing
        size_t ia_start = ia_starts[block];
        size_t ia_stop = ia_starts[block + 1];
        size_t ncols = ia_stop - ia_start;

        timer_on("DFMP2 Gia Read");
        psio_->read(file, "G(ia|Q)", (char*)Giap[0], sizeof(double) * ncols * naux, next_GIA, &next_GIA);
        timer_off("DFMP2 Gia Read");

        timer_on("DFMP2 Cia Read");
        psio_->read(file, "C(ia|Q)", (char*)Ciap[0], sizeof(double) * ncols * naux, next_CIA, &next_CIA);
        timer_off("DFMP2 Cia Read");

        // DEFINITION: G_PQ := G(ia|P)C(ia|Q)
        // Eq. 3, 23 of DiStasio
        // This is spin-summed becuase G(ia|P) is
        // This term is contracted against metric derivatives
        timer_on("DFMP2 g");
        C_DGEMM('T', 'N', naux, naux, ncols, 1.0, Giap[0], naux, Ciap[0], naux, 1.0, Gp[0], naux);
        timer_off("DFMP2 g");
    }

    G->save(psio_, file, Matrix::SaveType::SubBlocks);

    psio_->close(file, 1);
}
void DFMP2::apply_G_transpose(size_t file, size_t naux, size_t nia) {
    // Memory constraints
    size_t doubles = (size_t)(options_.get_double("DFMP2_MEM_FACTOR") * (memory_ / 8L));
    size_t max_nia = (doubles / naux);
    max_nia = (max_nia > nia ? nia : max_nia);
    max_nia = (max_nia < 1L ? 1L : max_nia);

    // Block sizing
    std::vector<size_t> ia_starts;
    ia_starts.push_back(0);
    for (size_t ia = 0L; ia < nia; ia += max_nia) {
        if (ia + max_nia >= nia) {
            ia_starts.push_back(nia);
        } else {
            ia_starts.push_back(ia + max_nia);
        }
    }
    // block_status(ia_starts, __FILE__,__LINE__);

    // Prestripe
    psio_->open(file, PSIO_OPEN_OLD);
    psio_address next_QIA = PSIO_ZERO;
    std::vector<double> temp(nia, 0);
    for (int Q = 0; Q < naux; Q++) {
        psio_->write(file, "G(Q|ia)", (char*)temp.data(), sizeof(double) * nia, next_QIA, &next_QIA);
    }
    std::vector<double>().swap(temp); // Dirty trick to clear temp memory immediately.
    next_QIA = PSIO_ZERO;
    psio_address next_IAQ = PSIO_ZERO;

    // Tensor blocks
    auto Qia = std::make_shared<Matrix>("G (Q|ia)", naux, max_nia);
    auto iaQ = std::make_shared<Matrix>("G (ia|Q)", max_nia, naux);
    auto Qiap = Qia->pointer();
    auto iaQp = iaQ->pointer();

    // Loop through blocks
    for (int block = 0; block < ia_starts.size() - 1; block++) {
        // Sizing
        size_t ia_start = ia_starts[block];
        size_t ia_stop = ia_starts[block + 1];
        size_t ncols = ia_stop - ia_start;

        // Read Gia
        timer_on("DFMP2 Gia Read");
        psio_->read(file, "G(ia|Q)", (char*)iaQp[0], sizeof(double) * ncols * naux, next_IAQ, &next_IAQ);
        timer_off("DFMP2 Gia Read");

        // Transpose
        for (int Q = 0; Q < naux; Q++) {
            C_DCOPY(ncols, &iaQp[0][Q], naux, Qiap[Q], 1);
        }

        // DEFINITION: G(Q|ia) := G(ia|Q)

        // Write Gia^\dagger
        timer_on("DFMP2 aiG Write");
        for (size_t Q = 0; Q < naux; Q++) {
            next_QIA = psio_get_address(PSIO_ZERO, sizeof(double) * (Q * nia + ia_start));
            psio_->write(file, "G(Q|ia)", (char*)Qiap[Q], sizeof(double) * ncols, next_QIA, &next_QIA);
        }
        timer_off("DFMP2 aiG Write");
    }
    psio_->close(file, 1);
}
void DFMP2::apply_B_transpose(size_t file, size_t naux, size_t naocc, size_t navir) {
    // Memory constraints
    size_t doubles = (size_t)(options_.get_double("DFMP2_MEM_FACTOR") * (memory_ / 8L));
    size_t rows = doubles / (1L * naocc * naux);
    int max_A = (rows <= 0L ? 1L : rows);
    max_A = (rows > navir ? navir : rows);

    // Block sizing
    std::vector<int> a_starts;
    a_starts.push_back(0);
    for (int a = 0; a < navir; a += max_A) {
        if (a + max_A >= navir) {
            a_starts.push_back(navir);
        } else {
            a_starts.push_back(a + max_A);
        }
    }
    // block_status(a_starts, __FILE__,__LINE__);

    // Buffers
    auto Bia = std::make_shared<Matrix>("B(ia|Q)", max_A * naocc, naux);
    double** Biap = Bia->pointer();

    // Loop through blocks
    psio_->open(file, PSIO_OPEN_OLD);
    psio_address next_BIA = PSIO_ZERO;
    psio_address next_BAI = PSIO_ZERO;
    for (int block = 0; block < a_starts.size() - 1; block++) {
        // Sizing
        int a_start = a_starts[block];
        int a_stop = a_starts[block + 1];
        int na = a_stop - a_start;

        for (int a = 0; a < na; a++) {
            for (int i = 0; i < naocc; i++) {
                next_BIA = psio_get_address(PSIO_ZERO, sizeof(double) * (i * navir * naux + (a + a_start) * naux));
                psio_->read(file, "B(ia|Q)", (char*)Biap[a * naocc + i], sizeof(double) * naux, next_BIA, &next_BIA);
            }
        }

        psio_->write(file, "B(ai|Q)", (char*)Biap[0], sizeof(double) * na * naocc * naux, next_BAI, &next_BAI);
    }
    psio_->close(file, 1);
}
void DFMP2::print_energies() {
    variables_["MP2 DOUBLES ENERGY"] = variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] +
                                       variables_["MP2 SAME-SPIN CORRELATION ENERGY"];
    variables_["MP2 CORRELATION ENERGY"] = variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] +
                                           variables_["MP2 SAME-SPIN CORRELATION ENERGY"] +
                                           variables_["MP2 SINGLES ENERGY"];
    variables_["MP2 TOTAL ENERGY"] = variables_["SCF TOTAL ENERGY"] + variables_["MP2 CORRELATION ENERGY"];

    variables_["SCS-MP2 CORRELATION ENERGY"] = 6.0/5.0 * variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] +
                                               1.0/3.0 * variables_["MP2 SAME-SPIN CORRELATION ENERGY"] +
                                                         variables_["MP2 SINGLES ENERGY"];
    variables_["SCS-MP2 TOTAL ENERGY"] = variables_["SCF TOTAL ENERGY"] + variables_["SCS-MP2 CORRELATION ENERGY"];
    variables_["CUSTOM SCS-MP2 CORRELATION ENERGY"] = oss_ * variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] +
                                                      sss_ * variables_["MP2 SAME-SPIN CORRELATION ENERGY"] +
                                                             variables_["MP2 SINGLES ENERGY"];
    variables_["CUSTOM SCS-MP2 TOTAL ENERGY"] = variables_["SCF TOTAL ENERGY"] + variables_["CUSTOM SCS-MP2 CORRELATION ENERGY"];

    outfile->Printf("\t-----------------------------------------------------------\n");
    outfile->Printf("\t ==================> DF-MP2 Energies <==================== \n");
    outfile->Printf("\t-----------------------------------------------------------\n");
    outfile->Printf("\t %-25s = %24.16f [Eh]\n", "Reference Energy", variables_["SCF TOTAL ENERGY"]);
    outfile->Printf("\t %-25s = %24.16f [Eh]\n", "Singles Energy", variables_["MP2 SINGLES ENERGY"]);
    outfile->Printf("\t %-25s = %24.16f [Eh]\n", "Same-Spin Energy", variables_["MP2 SAME-SPIN CORRELATION ENERGY"]);
    outfile->Printf("\t %-25s = %24.16f [Eh]\n", "Opposite-Spin Energy",
                    variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"]);
    outfile->Printf("\t %-25s = %24.16f [Eh]\n", "Correlation Energy", variables_["MP2 CORRELATION ENERGY"]);
    outfile->Printf("\t %-25s = %24.16f [Eh]\n", "Total Energy", variables_["MP2 TOTAL ENERGY"]);
    outfile->Printf("\t-----------------------------------------------------------\n");
    outfile->Printf("\t ================> DF-SCS-MP2 Energies <================== \n");
    outfile->Printf("\t-----------------------------------------------------------\n");
    outfile->Printf("\t %-25s = %24.16f [-]\n", "SCS Same-Spin Scale", sss_);
    outfile->Printf("\t %-25s = %24.16f [-]\n", "SCS Opposite-Spin Scale", oss_);
    outfile->Printf("\t %-25s = %24.16f [Eh]\n", "SCS Same-Spin Energy",
                    sss_ * variables_["MP2 SAME-SPIN CORRELATION ENERGY"]);
    outfile->Printf("\t %-25s = %24.16f [Eh]\n", "SCS Opposite-Spin Energy",
                    oss_ * variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"]);
    outfile->Printf("\t %-25s = %24.16f [Eh]\n", "SCS Correlation Energy", variables_["SCS-MP2 CORRELATION ENERGY"]);
    outfile->Printf("\t %-25s = %24.16f [Eh]\n", "SCS Total Energy", variables_["SCS-MP2 TOTAL ENERGY"]);
    outfile->Printf("\t-----------------------------------------------------------\n");
    outfile->Printf("\n");
}

void DFMP2::print_gradients() {
    std::vector<std::string> gradient_terms;
    gradient_terms.push_back("Nuclear");
    gradient_terms.push_back("Kinetic");
    gradient_terms.push_back("Potential");
    gradient_terms.push_back("Overlap");
    gradient_terms.push_back("Coulomb");
    gradient_terms.push_back("Exchange");
    gradient_terms.push_back("Correlation");
    gradient_terms.push_back("Total");

    if (print_ > 1) {
        for (int i = 0; i < gradient_terms.size(); i++) {
            if (gradients_.count(gradient_terms[i])) {
                gradients_[gradient_terms[i]]->print_atom_vector();
            }
        }
    } else {
        gradients_["Total"]->print_atom_vector();
    }
}

RDFMP2::RDFMP2(SharedWavefunction ref_wfn, Options& options, std::shared_ptr<PSIO> psio)
    : DFMP2(ref_wfn, options, psio) {
    common_init();
}
RDFMP2::~RDFMP2() {}
void RDFMP2::common_init() {
    Cfocc_ = Ca_subset("AO", "FROZEN_OCC");
    Caocc_ = Ca_subset("AO", "ACTIVE_OCC");
    Cavir_ = Ca_subset("AO", "ACTIVE_VIR");
    Cfvir_ = Ca_subset("AO", "FROZEN_VIR");

    eps_focc_ = epsilon_a_subset("AO", "FROZEN_OCC");
    eps_aocc_ = epsilon_a_subset("AO", "ACTIVE_OCC");
    eps_avir_ = epsilon_a_subset("AO", "ACTIVE_VIR");
    eps_fvir_ = epsilon_a_subset("AO", "FROZEN_VIR");
}
void RDFMP2::print_header() {
    int nthread = 1;
#ifdef _OPENMP
    nthread = Process::environment.get_n_threads();
#endif

    outfile->Printf("\t --------------------------------------------------------\n");
    outfile->Printf("\t                          DF-MP2                         \n");
    outfile->Printf("\t      2nd-Order Density-Fitted Moller-Plesset Theory     \n");
    outfile->Printf("\t              RMP2 Wavefunction, %3d Threads             \n", nthread);
    outfile->Printf("\t                                                         \n");
    outfile->Printf("\t        Rob Parrish, Justin Turney, Andy Simmonett,      \n");
    outfile->Printf("\t           Ed Hohenstein, and C. David Sherrill          \n");
    outfile->Printf("\t --------------------------------------------------------\n");
    outfile->Printf("\n");

    int focc = frzcpi_.sum();
    int fvir = frzvpi_.sum();
    int aocc = Caocc_->colspi()[0];
    int avir = Cavir_->colspi()[0];
    int occ = focc + aocc;
    int vir = fvir + avir;

    if (print_ >= 1) {
        outfile->Printf("   => Auxiliary Basis Set <=\n\n");
        ribasis_->print_by_level("outfile", print_);
    }

    outfile->Printf("\t --------------------------------------------------------\n");
    outfile->Printf("\t                 NBF = %5d, NAUX = %5d\n", basisset_->nbf(), ribasis_->nbf());
    outfile->Printf("\t --------------------------------------------------------\n");
    outfile->Printf("\t %7s %7s %7s %7s %7s %7s %7s\n", "CLASS", "FOCC", "OCC", "AOCC", "AVIR", "VIR", "FVIR");
    outfile->Printf("\t %7s %7d %7d %7d %7d %7d %7d\n", "PAIRS", focc, occ, aocc, avir, vir, fvir);
    outfile->Printf("\t --------------------------------------------------------\n\n");
}
void RDFMP2::form_Aia() {
    // ERI objects
    int nthread = 1;
#ifdef _OPENMP
    if (options_.get_int("DF_INTS_NUM_THREADS") == 0) {
        nthread = Process::environment.get_n_threads();
    } else {
        nthread = options_.get_int("DF_INTS_NUM_THREADS");
    }
#endif
    IntegralFactory factory(ribasis_, BasisSet::zero_ao_basis_set(), basisset_, basisset_);
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    std::vector<const double*> buffer;
    for (int thread = 0; thread < nthread; thread++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(factory.eri()));
        buffer.push_back(eri[thread]->buffer());
    }

    // Schwarz Sieve
    const auto& shell_pairs = eri[0]->shell_pairs();
    const size_t npairs = shell_pairs.size();

    // Sizing
    int nso = basisset_->nbf();
    int naux = ribasis_->nbf();
    int naocc = Caocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];
    int maxQ = ribasis_->max_function_per_shell();

    // Max block size in naux
    size_t Amn_cost_per_row = nso * (size_t)nso;
    size_t Ami_cost_per_row = nso * (size_t)naocc;
    size_t Aia_cost_per_row = naocc * (size_t)navir;
    size_t total_cost_per_row = Amn_cost_per_row + Ami_cost_per_row + Aia_cost_per_row;
    size_t doubles = ((size_t)(options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    size_t max_temp = doubles / (total_cost_per_row);
    int max_naux = (max_temp > (size_t)naux ? naux : max_temp);
    max_naux = (max_naux < maxQ ? maxQ : max_naux);

    // Block extents
    std::vector<int> block_Q_starts;
    int counter = 0;
    block_Q_starts.push_back(0);
    for (int Q = 0; Q < ribasis_->nshell(); Q++) {
        int nQ = ribasis_->shell(Q).nfunction();
        if (counter + nQ > max_naux) {
            counter = 0;
            block_Q_starts.push_back(Q);
        }
        counter += nQ;
    }
    block_Q_starts.push_back(ribasis_->nshell());
    // block_status(block_Q_starts, __FILE__,__LINE__);

    // Tensor blocks
    auto Amn = std::make_shared<Matrix>("(A|mn) Block", max_naux, nso * (size_t)nso);
    auto Ami = std::make_shared<Matrix>("(A|mi) Block", max_naux, nso * (size_t)naocc);
    auto Aia = std::make_shared<Matrix>("(A|ia) Block", max_naux, naocc * (size_t)navir);
    auto Amnp = Amn->pointer();
    auto Amip = Ami->pointer();
    auto Aiap = Aia->pointer();

    // C Matrices
    auto Caoccp = Caocc_->pointer();
    auto Cavirp = Cavir_->pointer();

    psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_NEW);
    psio_address next_AIA = PSIO_ZERO;

    // Loop over blocks of Qshell
    for (int block = 0; block < block_Q_starts.size() - 1; block++) {
        // Block sizing/offsets
        int Qstart = block_Q_starts[block];
        int Qstop = block_Q_starts[block + 1];
        int qoff = ribasis_->shell(Qstart).function_index();
        int nrows = (Qstop == ribasis_->nshell()
                         ? ribasis_->nbf() - ribasis_->shell(Qstart).function_index()
                         : ribasis_->shell(Qstop).function_index() - ribasis_->shell(Qstart).function_index());

        // Clear Amn for Schwarz sieve
        ::memset((void*)Amnp[0], '\0', sizeof(double) * nrows * nso * nso);

        // Compute TEI tensor block (A|mn)
        timer_on("DFMP2 (A|mn)");
#pragma omp parallel for schedule(dynamic) num_threads(nthread)
        for (long int QMN = 0L; QMN < (Qstop - Qstart) * (size_t)npairs; QMN++) {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            int Q = QMN / npairs + Qstart;
            int MN = QMN % npairs;

            auto pair = shell_pairs[MN];
            int M = pair.first;
            int N = pair.second;

            int nq = ribasis_->shell(Q).nfunction();
            int nm = basisset_->shell(M).nfunction();
            int nn = basisset_->shell(N).nfunction();

            int sq = ribasis_->shell(Q).function_index();
            int sm = basisset_->shell(M).function_index();
            int sn = basisset_->shell(N).function_index();

            eri[thread]->compute_shell(Q, 0, M, N);
            buffer[thread] = eri[thread]->buffer();

            for (int oq = 0; oq < nq; oq++) {
                for (int om = 0; om < nm; om++) {
                    for (int on = 0; on < nn; on++) {
                        Amnp[sq + oq - qoff][(om + sm) * nso + (on + sn)] =
                            Amnp[sq + oq - qoff][(on + sn) * nso + (om + sm)] =
                                buffer[thread][oq * nm * nn + om * nn + on];
                    }
                }
            }
        }
        timer_off("DFMP2 (A|mn)");

        // Compute (A|mi) tensor block (A|mn) C_ni
        timer_on("DFMP2 (A|mn)C_mi");
        C_DGEMM('N', 'N', nrows * (size_t)nso, naocc, nso, 1.0, Amnp[0], nso, Caoccp[0], naocc, 0.0, Amip[0], naocc);
        timer_off("DFMP2 (A|mn)C_mi");

        // Compute (A|ia) tensor block (A|ia) = (A|mi) C_ma
        timer_on("DFMP2 (A|mi)C_na");
#pragma omp parallel for
        for (int row = 0; row < nrows; row++) {
            C_DGEMM('T', 'N', naocc, navir, nso, 1.0, Amip[row], naocc, Cavirp[0], navir, 0.0, Aiap[row], navir);
        }
        timer_off("DFMP2 (A|mi)C_na");

        // Stripe (A|ia) out to disk
        timer_on("DFMP2 Aia Write");
        psio_->write(PSIF_DFMP2_AIA, "A(Q|ia)", (char*)Aiap[0], sizeof(double) * nrows * naocc * navir, next_AIA,
                     &next_AIA);
        timer_off("DFMP2 Aia Write");
    }

    psio_->close(PSIF_DFMP2_AIA, 1);
}
void RDFMP2::form_Bia() {
    SharedMatrix Jm12 = form_inverse_metric();
    apply_fitting(Jm12, PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_->colspi()[0] * (size_t)Cavir_->colspi()[0]);
}
void RDFMP2::form_Bia_Cia() {
    SharedMatrix Jm12 = form_inverse_metric();
    apply_fitting(Jm12, PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_->colspi()[0] * (size_t)Cavir_->colspi()[0]);
    apply_fitting_grad(Jm12, PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_->colspi()[0] * (size_t)Cavir_->colspi()[0]);
}
void RDFMP2::form_Bia_transpose() {
    apply_B_transpose(PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_->colspi()[0], (size_t)Cavir_->colspi()[0]);
}
void RDFMP2::form_energy() {
    // Energy registers
    double e_ss = 0.0;
    double e_os = 0.0;

    // Sizing
    int naux = ribasis_->nbf();
    int naocc = Caocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];

    // Thread considerations
    int nthread = 1;
#ifdef _OPENMP
    nthread = Process::environment.get_n_threads();
#endif

    // Memory
    size_t Iab_memory = navir * (size_t)navir;
    size_t Qa_memory = naux * (size_t)navir;
    size_t doubles = ((size_t)(options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    if (doubles < nthread * Iab_memory) {
        throw PSIEXCEPTION("DFMP2: Insufficient memory for Iab buffers. Reduce OMP Threads or increase memory.");
    }
    size_t remainder = doubles - nthread * Iab_memory;
    size_t max_i = remainder / (2L * Qa_memory);
    max_i = (max_i > naocc ? naocc : max_i);
    max_i = (max_i < 1L ? 1L : max_i);

    // Blocks
    std::vector<size_t> i_starts;
    i_starts.push_back(0L);
    for (size_t i = 0; i < naocc; i += max_i) {
        if (i + max_i >= naocc) {
            i_starts.push_back(naocc);
        } else {
            i_starts.push_back(i + max_i);
        }
    }
    // block_status(i_starts, __FILE__,__LINE__);

    // Tensor blocks
    auto Bia = std::make_shared<Matrix>("B(ia|Q)", max_i * (size_t)navir, naux);
    auto Bjb = std::make_shared<Matrix>("B(jb|Q)", max_i * (size_t)navir, naux);
    double** Biap = Bia->pointer();
    double** Bjbp = Bjb->pointer();

    std::vector<SharedMatrix> Iab;
    for (int i = 0; i < nthread; i++) {
        Iab.push_back(std::make_shared<Matrix>("Iab", navir, navir));
    }

    double* eps_aoccp = eps_aocc_->pointer();
    double* eps_avirp = eps_avir_->pointer();

    // Loop through pairs of blocks
    psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
    psio_address next_BIA = PSIO_ZERO;
    for (int block_i = 0; block_i < i_starts.size() - 1; block_i++) {
        // Sizing
        size_t istart = i_starts[block_i];
        size_t istop = i_starts[block_i + 1];
        size_t ni = istop - istart;

        // Read iaQ chunk
        timer_on("DFMP2 Bia Read");
        next_BIA = psio_get_address(PSIO_ZERO, sizeof(double) * (istart * navir * naux));
        psio_->read(PSIF_DFMP2_AIA, "B(ia|Q)", (char*)Biap[0], sizeof(double) * (ni * navir * naux), next_BIA,
                    &next_BIA);
        timer_off("DFMP2 Bia Read");

        for (int block_j = 0; block_j <= block_i; block_j++) {
            // Sizing
            size_t jstart = i_starts[block_j];
            size_t jstop = i_starts[block_j + 1];
            size_t nj = jstop - jstart;

            // Read iaQ chunk (if unique)
            timer_on("DFMP2 Bia Read");
            if (block_i == block_j) {
                ::memcpy((void*)Bjbp[0], (void*)Biap[0], sizeof(double) * (ni * navir * naux));
            } else {
                next_BIA = psio_get_address(PSIO_ZERO, sizeof(double) * (jstart * navir * naux));
                psio_->read(PSIF_DFMP2_AIA, "B(ia|Q)", (char*)Bjbp[0], sizeof(double) * (nj * navir * naux), next_BIA,
                            &next_BIA);
            }
            timer_off("DFMP2 Bia Read");

#pragma omp parallel for schedule(dynamic) num_threads(nthread) reduction(+ : e_ss, e_os)
            for (long int ij = 0L; ij < ni * nj; ij++) {
                // Sizing
                size_t i = ij / nj + istart;
                size_t j = ij % nj + jstart;
                if (j > i) continue;

                double perm_factor = (i == j ? 1.0 : 2.0);

                // Which thread is this?
                int thread = 0;
#ifdef _OPENMP
                thread = omp_get_thread_num();
#endif
                double** Iabp = Iab[thread]->pointer();

                // Form the integral block (ia|jb) = (ia|Q)(Q|jb)
                C_DGEMM('N', 'T', navir, navir, naux, 1.0, Biap[(i - istart) * navir], naux, Bjbp[(j - jstart) * navir],
                        naux, 0.0, Iabp[0], navir);

                // Add the MP2 energy contributions
                for (int a = 0; a < navir; a++) {
                    for (int b = 0; b < navir; b++) {
                        double iajb = Iabp[a][b];
                        double ibja = Iabp[b][a];
                        double denom = -perm_factor / (eps_avirp[a] + eps_avirp[b] - eps_aoccp[i] - eps_aoccp[j]);

                        e_ss += (iajb * iajb - iajb * ibja) * denom;
                        e_os += (iajb * iajb) * denom;
                    }
                }
            }
        }
    }
    psio_->close(PSIF_DFMP2_AIA, 0);

    variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = e_ss;
    variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = e_os;
}
void RDFMP2::form_Pab() {
    // Energy registers
    double e_ss = 0.0;
    double e_os = 0.0;

    // Sizing
    int naux = ribasis_->nbf();
    int naocc = Caocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];

    // Thread considerations
    int nthread = 1;
#ifdef _OPENMP
    nthread = Process::environment.get_n_threads();
#endif

    // Memory
    size_t doubles = static_cast<size_t>(options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L);
    doubles -= static_cast<double>(navir) * static_cast<double>(navir);
    double C = -(double)doubles;
    double B = 4.0 * navir * naux;
    double A = 2.0 * navir * (double)navir;

    int max_i = (int)((-B + sqrt(B * B - 4.0 * A * C)) / (2.0 * A));
    if (max_i <= 0) {
        throw PSIEXCEPTION("Not enough memory in DFMP2");
    }
    max_i = (max_i <= 0 ? 1 : max_i);
    max_i = (max_i > naocc ? naocc : max_i);

    // Blocks
    std::vector<size_t> i_starts;
    i_starts.push_back(0L);
    for (size_t i = 0; i < naocc; i += max_i) {
        if (i + max_i >= naocc) {
            i_starts.push_back(naocc);
        } else {
            i_starts.push_back(i + max_i);
        }
    }
    // block_status(i_starts, __FILE__,__LINE__);

    // 2-Index Tensor blocks
    auto Pab = std::make_shared<Matrix>("P_ab", navir, navir);
    auto Pabp = Pab->pointer();

    // 3-Index Tensor blocks
    auto Bia = std::make_shared<Matrix>("B(ia|Q)", max_i * (size_t)navir, naux);
    auto Bjb = std::make_shared<Matrix>("B(jb|Q)", max_i * (size_t)navir, naux);
    auto Gia = std::make_shared<Matrix>("Gia", max_i * (size_t)navir, naux);
    auto Cjb = std::make_shared<Matrix>("C(jb|Q)", max_i * (size_t)navir, naux);

    auto Biap = Bia->pointer();
    auto Bjbp = Bjb->pointer();
    auto Giap = Gia->pointer();
    auto Cjbp = Cjb->pointer();

    // 4-index Tensor blocks
    auto I = std::make_shared<Matrix>("I", max_i * (size_t)navir, max_i * (size_t)navir);
    auto T = std::make_shared<Matrix>("T", max_i * (size_t)navir, max_i * (size_t)navir);
    auto Ip = I->pointer();
    auto Tp = T->pointer();

    size_t nIv = max_i * (size_t)navir;

    auto eps_aoccp = eps_aocc_->pointer();
    auto eps_avirp = eps_avir_->pointer();

    // Loop through pairs of blocks
    psio_address next_QIA = PSIO_ZERO;
    psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
    for (int block_i = 0; block_i < i_starts.size() - 1; block_i++) {
        // Sizing
        size_t istart = i_starts[block_i];
        size_t istop = i_starts[block_i + 1];
        size_t ni = istop - istart;

        // Read iaQ chunk
        timer_on("DFMP2 Bia Read");
        next_QIA = psio_get_address(PSIO_ZERO, sizeof(double) * (istart * navir * naux));
        psio_->read(PSIF_DFMP2_AIA, "B(ia|Q)", (char*)Biap[0], sizeof(double) * (ni * navir * naux), next_QIA,
                    &next_QIA);
        timer_off("DFMP2 Bia Read");

        // Zero Gamma for current ia
        Gia->zero();

        for (int block_j = 0; block_j < i_starts.size() - 1; block_j++) {
            // Sizing
            size_t jstart = i_starts[block_j];
            size_t jstop = i_starts[block_j + 1];
            size_t nj = jstop - jstart;

            // Read iaQ chunk (if unique)
            timer_on("DFMP2 Bia Read");
            if (block_i == block_j) {
                ::memcpy((void*)Bjbp[0], (void*)Biap[0], sizeof(double) * (ni * navir * naux));
            } else {
                next_QIA = psio_get_address(PSIO_ZERO, sizeof(double) * (jstart * navir * naux));
                psio_->read(PSIF_DFMP2_AIA, "B(ia|Q)", (char*)Bjbp[0], sizeof(double) * (nj * navir * naux), next_QIA,
                            &next_QIA);
            }
            timer_off("DFMP2 Bia Read");

            // Read iaC chunk
            timer_on("DFMP2 Cia Read");
            next_QIA = psio_get_address(PSIO_ZERO, sizeof(double) * (jstart * navir * naux));
            psio_->read(PSIF_DFMP2_AIA, "C(ia|Q)", (char*)Cjbp[0], sizeof(double) * (nj * navir * naux), next_QIA,
                        &next_QIA);
            timer_off("DFMP2 Cia Read");

            // Form the integrals (ia|jb) = B_ia^Q B_jb^Q
            timer_on("DFMP2 I");
            C_DGEMM('N', 'T', ni * (size_t)navir, nj * (size_t)navir, naux, 1.0, Biap[0], naux, Bjbp[0], naux, 0.0,
                    Ip[0], nIv);
            timer_off("DFMP2 I");

            timer_on("DFMP2 T2");
// Form the T amplitudes t_ia^jb = [2(ia|jb) - (ib|ja)] / (e_a + e_b - e_i - e_j)
// Form the I amplitudes I_ia^jb = (ia|jb) / (e_a + e_b - e_i - e_j);
// Form the energy contributions
#pragma omp parallel for schedule(dynamic) num_threads(nthread) reduction(+ : e_ss, e_os)
            for (long int ij = 0L; ij < ni * nj; ij++) {
                // Sizing
                size_t i_local = ij / nj;
                size_t j_local = ij % nj;
                size_t i = i_local + istart;
                size_t j = j_local + jstart;

                // Add the MP2 energy contributions and form the T amplitudes in place
                for (int a = 0; a < navir; a++) {
                    for (int b = 0; b <= a; b++) {
                        double iajb = Ip[i_local * navir + a][j_local * navir + b];
                        double ibja = Ip[i_local * navir + b][j_local * navir + a];
                        double denom = -1.0 / (eps_avirp[a] + eps_avirp[b] - eps_aoccp[i] - eps_aoccp[j]);
                        Tp[i_local * navir + a][j_local * navir + b] = denom * (2.0 * iajb - ibja);
                        Tp[i_local * navir + b][j_local * navir + a] = denom * (2.0 * ibja - iajb);
                        Ip[i_local * navir + a][j_local * navir + b] = denom * (iajb);
                        Ip[i_local * navir + b][j_local * navir + a] = denom * (ibja);

                        e_ss += (iajb * iajb - iajb * ibja) * denom;
                        e_os += (iajb * iajb) * denom;

                        if (a != b) {
                            e_ss += (ibja * ibja - ibja * iajb) * denom;
                            e_os += (ibja * ibja) * denom;
                        }
                    }
                }
            }
            timer_off("DFMP2 T2");

            // DEFINITION: G(ia|Q) = t^ab_ij C(jb|Q)
            // Eq. 22 (spin-summed case of Eq. 2) of DiStasio.
            // The three-index intermediate that is contracted against derivatives of A integrals.
            // Also used to construct the two-index intermediate contracted against metric derivatives.
            timer_on("DFMP2 G");
            C_DGEMM('N', 'N', ni * (size_t)navir, naux, nj * (size_t)navir, 2.0, Tp[0], nIv, Cjbp[0], naux, 1.0,
                    Giap[0], naux);
            timer_off("DFMP2 G");

            // Sort the gimp column blocks, if gimp occurred. The idea is to get a contiguous iajb tensor
            if (nj != max_i) {
                size_t counter = 0L;
                for (long int ind = 0L; ind < ni * (size_t)navir; ind++) {
                    ::memmove((void*)&Tp[0][counter], (void*)Tp[ind], sizeof(double) * nj * (size_t)navir);
                    ::memmove((void*)&Ip[0][counter], (void*)Ip[ind], sizeof(double) * nj * (size_t)navir);
                    counter += nj * (size_t)navir;
                }
            }

            // DEFINITION: P_ab := + 1/2 t^ac_ij t^bc_ij
            // Eq. 8 of DiStasio summed over both spin cases. We use spin-adapted intermediates.
            // The *2 in DGEMM accounts for spin summation. The /2 in Eq. 8 is accounted for when defining T and I.
            // This quantity is more commonly known as the Virt/Virt pure correlation OPDM.
            timer_on("DFMP2 Pab");
            C_DGEMM('T', 'N', navir, navir, ni * (size_t)nj * navir, 2.0, Tp[0], navir, Ip[0], navir, 1.0, Pabp[0],
                    navir);
            timer_off("DFMP2 Pab");
        }

        // Write iaG chunk
        timer_on("DFMP2 Gia Write");
        next_QIA = psio_get_address(PSIO_ZERO, sizeof(double) * (istart * navir * naux));
        psio_->write(PSIF_DFMP2_AIA, "G(ia|Q)", (char*)Giap[0], sizeof(double) * (ni * navir * naux), next_QIA,
                     &next_QIA);
        timer_off("DFMP2 Gia Write");
    }

    psio_->write_entry(PSIF_DFMP2_AIA, "P_ab", (char*)Pabp[0], sizeof(double) * navir * navir);

    psio_->close(PSIF_DFMP2_AIA, 1);

    variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = e_ss;
    variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = e_os;
}
void RDFMP2::form_Pij() {
    // Energy registers
    double e_ss = 0.0;
    double e_os = 0.0;

    // Sizing
    int naux = ribasis_->nbf();
    int naocc = Caocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];

    // Thread considerations
    int nthread = 1;
#ifdef _OPENMP
    nthread = Process::environment.get_n_threads();
#endif

    // Memory
    size_t doubles = static_cast<size_t>(options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L);
    doubles -= naocc * naocc;
    double C = -(double)doubles;
    double B = 2.0 * naocc * naux;
    double A = 2.0 * naocc * (double)naocc;

    int max_a = (int)((-B + sqrt(B * B - 4.0 * A * C)) / (2.0 * A));
    if (max_a <= 0) {
        throw PSIEXCEPTION("Not enough memory in DFMP2");
    }

    max_a = (max_a > navir ? navir : max_a);

    // Blocks
    std::vector<size_t> a_starts;
    a_starts.push_back(0L);
    for (size_t a = 0; a < navir; a += max_a) {
        if (a + max_a >= navir) {
            a_starts.push_back(navir);
        } else {
            a_starts.push_back(a + max_a);
        }
    }
    // block_status(a_starts, __FILE__,__LINE__);

    // 2-Index Tensor blocks
    auto Pij = std::make_shared<Matrix>("P_ij", naocc, naocc);
    auto Pijp = Pij->pointer();

    // 3-Index Tensor blocks
    auto Bia = std::make_shared<Matrix>("B(ia|Q)", max_a * (size_t)naocc, naux);
    auto Bjb = std::make_shared<Matrix>("B(jb|Q)", max_a * (size_t)naocc, naux);

    auto Biap = Bia->pointer();
    auto Bjbp = Bjb->pointer();

    // 4-index Tensor blocks
    auto I = std::make_shared<Matrix>("I", max_a * (size_t)naocc, max_a * (size_t)naocc);
    auto T = std::make_shared<Matrix>("T", max_a * (size_t)naocc, max_a * (size_t)naocc);
    auto Ip = I->pointer();
    auto Tp = T->pointer();

    size_t nVi = max_a * (size_t)naocc;

    auto eps_aoccp = eps_aocc_->pointer();
    auto eps_avirp = eps_avir_->pointer();

    // Loop through pairs of blocks
    psio_address next_BAI = PSIO_ZERO;
    psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
    for (int block_a = 0; block_a < a_starts.size() - 1; block_a++) {
        // Sizing
        size_t astart = a_starts[block_a];
        size_t astop = a_starts[block_a + 1];
        size_t na = astop - astart;

        // Read iaQ chunk
        timer_on("DFMP2 Bai Read");
        next_BAI = psio_get_address(PSIO_ZERO, sizeof(double) * (astart * naocc * naux));
        psio_->read(PSIF_DFMP2_AIA, "B(ai|Q)", (char*)Biap[0], sizeof(double) * (na * naocc * naux), next_BAI,
                    &next_BAI);
        timer_off("DFMP2 Bai Read");

        for (int block_b = 0; block_b < a_starts.size() - 1; block_b++) {
            // Sizing
            size_t bstart = a_starts[block_b];
            size_t bstop = a_starts[block_b + 1];
            size_t nb = bstop - bstart;

            // Read iaQ chunk (if unique)
            timer_on("DFMP2 Qai Read");
            if (block_a == block_b) {
                ::memcpy((void*)Bjbp[0], (void*)Biap[0], sizeof(double) * (na * naocc * naux));
            } else {
                next_BAI = psio_get_address(PSIO_ZERO, sizeof(double) * (bstart * naocc * naux));
                psio_->read(PSIF_DFMP2_AIA, "B(ai|Q)", (char*)Bjbp[0], sizeof(double) * (nb * naocc * naux), next_BAI,
                            &next_BAI);
            }
            timer_off("DFMP2 Qai Read");

            // Form the integrals (ia|jb) = B_ia^Q B_jb^Q
            timer_on("DFMP2 I");
            C_DGEMM('N', 'T', na * (size_t)naocc, nb * (size_t)naocc, naux, 1.0, Biap[0], naux, Bjbp[0], naux, 0.0,
                    Ip[0], nVi);
            timer_off("DFMP2 I");

            timer_on("DFMP2 T2");
// Form the T amplitudes t_ia^jb = [2(ia|jb) - (ib|ja)] / (e_a + e_b - e_i - e_j)
// Form the I amplitudes I_ia^jb = (ia|jb) / (e_a + e_b - e_i - e_j);
// Form the energy contributions
#pragma omp parallel for schedule(dynamic) num_threads(nthread) reduction(+ : e_ss, e_os)
            for (long int ab = 0L; ab < na * nb; ab++) {
                // Sizing
                size_t a_local = ab / nb;
                size_t b_local = ab % nb;
                size_t a = a_local + astart;
                size_t b = b_local + bstart;

                // Add the MP2 energy contributions and form the T amplitudes in place
                for (int i = 0; i < naocc; i++) {
                    for (int j = 0; j <= i; j++) {
                        double iajb = Ip[a_local * naocc + i][b_local * naocc + j];
                        double ibja = Ip[a_local * naocc + j][b_local * naocc + i];
                        double denom = -1.0 / (eps_avirp[a] + eps_avirp[b] - eps_aoccp[i] - eps_aoccp[j]);
                        Tp[a_local * naocc + i][b_local * naocc + j] = denom * (2.0 * iajb - ibja);
                        Tp[a_local * naocc + j][b_local * naocc + i] = denom * (2.0 * ibja - iajb);
                        Ip[a_local * naocc + i][b_local * naocc + j] = denom * (iajb);
                        Ip[a_local * naocc + j][b_local * naocc + i] = denom * (ibja);

                        e_ss += (iajb * iajb - iajb * ibja) * denom;
                        e_os += (iajb * iajb) * denom;

                        if (i != j) {
                            e_ss += (ibja * ibja - ibja * iajb) * denom;
                            e_os += (ibja * ibja) * denom;
                        }
                    }
                }
            }
            timer_off("DFMP2 T2");

            // Sort the gimp column blocks, if gimp occurred. The idea is to get a contiguous iajb tensor
            if (nb != max_a) {
                size_t counter = 0L;
                for (long int ind = 0L; ind < na * (size_t)naocc; ind++) {
                    ::memmove((void*)&Tp[0][counter], (void*)Tp[ind], sizeof(double) * nb * (size_t)naocc);
                    ::memmove((void*)&Ip[0][counter], (void*)Ip[ind], sizeof(double) * nb * (size_t)naocc);
                    counter += nb * (size_t)naocc;
                }
            }

            // DEFINITION: P_ij := - 1/2 * t^ab_ik t^ab_jk
            // Eq. 7 of DiStasio summed over both spin cases. We use spin-adapted intermediates.
            // The *2 in DGEMM accounts for spin summation. The /2 in Eq. 7 is accounted for when defining T and I.
            // This quantity is more commonly known as the Occ/Occ pure correlation OPDM.
            timer_on("DFMP2 Pij");
            C_DGEMM('T', 'N', naocc, naocc, na * (size_t)nb * naocc, -2.0, Tp[0], naocc, Ip[0], naocc, 1.0, Pijp[0],
                    naocc);
            timer_off("DFMP2 Pij");
        }
    }

    psio_->write_entry(PSIF_DFMP2_AIA, "P_ij", (char*)Pijp[0], sizeof(double) * naocc * naocc);

    psio_->close(PSIF_DFMP2_AIA, 1);
}
void RDFMP2::form_gamma() {
    apply_gamma(PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_->colspi()[0] * (size_t)Cavir_->colspi()[0]);
}
void RDFMP2::form_G_transpose() {
    apply_G_transpose(PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_->colspi()[0] * (size_t)Cavir_->colspi()[0]);
}
void RDFMP2::form_AB_x_terms() {
    auto naux = ribasis_->nbf();

    auto V = std::make_shared<Matrix>("G_PQ", naux, naux);
    V->load(psio_, PSIF_DFMP2_AIA, Matrix::SaveType::SubBlocks);
    V->hermitivitize();
    V->scale(2);
    std::map<std::string, SharedMatrix> densities = {{"(A|B)^x", V}};

    auto results = mintshelper_->metric_grad(densities, "DF_BASIS_MP2");

    for (const auto& kv: results) {
        gradients_[kv.first] = kv.second;
    }
}
void RDFMP2::form_Amn_x_terms() {
    // => Sizing <= //

    int natom = basisset_->molecule()->natom();
    int nso = basisset_->nbf();
    int naocc = Caocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];
    int nia = Caocc_->colspi()[0] * Cavir_->colspi()[0];
    int naux = ribasis_->nbf();

    // => Thread Count <= //

    int num_threads = 1;
#ifdef _OPENMP
    num_threads = Process::environment.get_n_threads();
#endif

    // => Integrals <= //

    IntegralFactory rifactory(ribasis_, BasisSet::zero_ao_basis_set(), basisset_, basisset_);
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < num_threads; t++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory.eri(1)));
    }

    // => ERI Sieve <= //
    const auto& shell_pairs = eri[0]->shell_pairs();
    int npairs = shell_pairs.size();

    // => Gradient Contribution <= //

    gradients_["(A|mn)^x"] = std::make_shared<Matrix>("(A|mn)^x Gradient", natom, 3);

    // => Memory Constraints <= //

    size_t memory = ((size_t)(options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    int max_rows;
    int maxP = ribasis_->max_function_per_shell();
    size_t row_cost = 0L;
    row_cost += nso * (size_t)nso;
    row_cost += nso * (size_t)naocc;
    row_cost += naocc * (size_t)navir;
    size_t rows = memory / row_cost;
    rows = (rows > naux ? naux : rows);
    rows = (rows < maxP ? maxP : rows);
    max_rows = (int)rows;

    // => Block Sizing <= //

    std::vector<int> Pstarts;
    int counter = 0;
    Pstarts.push_back(0);
    for (int P = 0; P < ribasis_->nshell(); P++) {
        int nP = ribasis_->shell(P).nfunction();
        if (counter + nP > max_rows) {
            counter = 0;
            Pstarts.push_back(P);
        }
        counter += nP;
    }
    Pstarts.push_back(ribasis_->nshell());
    // block_status(Pstarts, __FILE__,__LINE__);

    // => Temporary Buffers <= //

    auto Gia = std::make_shared<Matrix>("Gia", max_rows, naocc * navir);
    auto Gmi = std::make_shared<Matrix>("Gmi", max_rows, nso * naocc);
    auto Gmn = std::make_shared<Matrix>("Gmn", max_rows, nso * (size_t)nso);

    double** Giap = Gia->pointer();
    double** Gmip = Gmi->pointer();
    double** Gmnp = Gmn->pointer();

    double** Caoccp = Caocc_->pointer();
    double** Cavirp = Cavir_->pointer();

    // => Temporary Gradients <= //

    std::vector<SharedMatrix> Ktemps;
    for (int t = 0; t < num_threads; t++) {
        Ktemps.push_back(std::make_shared<Matrix>("Ktemp", natom, 3));
    }

    // => PSIO <= //

    psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;

    // => Master Loop <= //

    for (int block = 0; block < Pstarts.size() - 1; block++) {
        // > Sizing < //

        int Pstart = Pstarts[block];
        int Pstop = Pstarts[block + 1];
        int NP = Pstop - Pstart;

        int pstart = ribasis_->shell(Pstart).function_index();
        int pstop = (Pstop == ribasis_->nshell() ? naux : ribasis_->shell(Pstop).function_index());
        int np = pstop - pstart;

        // > G_ia^P -> G_mn^P < //

        psio_->read(PSIF_DFMP2_AIA, "G(Q|ia)", (char*)Giap[0], sizeof(double) * np * nia, next_AIA, &next_AIA);

#pragma omp parallel for num_threads(num_threads)
        for (int p = 0; p < np; p++) {
            C_DGEMM('N', 'T', nso, naocc, navir, 1.0, Cavirp[0], navir, Giap[p], navir, 0.0, Gmip[p], naocc);
        }

        C_DGEMM('N', 'T', np * (size_t)nso, nso, naocc, 1.0, Gmip[0], naocc, Caoccp[0], naocc, 0.0, Gmnp[0], nso);

        // On Prefactors:
        // One factor of 2 is built into the definition of the term.
        // Combined, the permutational factor and 0.5 * (G(P|mn) + G(P|nm)) account for us summing over ordered shell pairs. 
        // Spin cases are accounted for because G(Q|ia) is spin-summed

// > Integrals < //
#pragma omp parallel for schedule(dynamic) num_threads(num_threads)
        for (long int PMN = 0L; PMN < static_cast<long int>(NP) * npairs; PMN++) {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            int P = PMN / npairs + Pstart;
            int MN = PMN % npairs;
            int M = shell_pairs[MN].first;
            int N = shell_pairs[MN].second;

            eri[thread]->compute_shell_deriv1(P, 0, M, N);

            const auto buffer = eri[thread]->buffer();
            const auto &buffers = eri[thread]->buffers();

            int nP = ribasis_->shell(P).nfunction();
            int cP = ribasis_->shell(P).ncartesian();
            int aP = ribasis_->shell(P).ncenter();
            int oP = ribasis_->shell(P).function_index() - pstart;

            int nM = basisset_->shell(M).nfunction();
            int cM = basisset_->shell(M).ncartesian();
            int aM = basisset_->shell(M).ncenter();
            int oM = basisset_->shell(M).function_index();

            int nN = basisset_->shell(N).nfunction();
            int cN = basisset_->shell(N).ncartesian();
            int aN = basisset_->shell(N).ncenter();
            int oN = basisset_->shell(N).function_index();

            const double* Px = buffers[0];
            const double* Py = buffers[1];
            const double* Pz = buffers[2];
            const double* Mx = buffers[3];
            const double* My = buffers[4];
            const double* Mz = buffers[5];
            const double* Nx = buffers[6];
            const double* Ny = buffers[7];
            const double* Nz = buffers[8];

            double perm = (M == N ? 1.0 : 2.0);

            double** grad_Kp = Ktemps[thread]->pointer();

            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {
                        double Jval =
                            2.0 * perm *
                            (0.5 * (Gmnp[p + oP][(m + oM) * nso + (n + oN)] + Gmnp[p + oP][(n + oN) * nso + (m + oM)]));
                        grad_Kp[aP][0] += Jval * (*Px);
                        grad_Kp[aP][1] += Jval * (*Py);
                        grad_Kp[aP][2] += Jval * (*Pz);
                        grad_Kp[aM][0] += Jval * (*Mx);
                        grad_Kp[aM][1] += Jval * (*My);
                        grad_Kp[aM][2] += Jval * (*Mz);
                        grad_Kp[aN][0] += Jval * (*Nx);
                        grad_Kp[aN][1] += Jval * (*Ny);
                        grad_Kp[aN][2] += Jval * (*Nz);

                        Px++;
                        Py++;
                        Pz++;
                        Mx++;
                        My++;
                        Mz++;
                        Nx++;
                        Ny++;
                        Nz++;
                    }
                }
            }
        }
    }

    // => Temporary Gradient Reduction <= //

    for (int t = 0; t < num_threads; t++) {
        gradients_["(A|mn)^x"]->add(Ktemps[t]);
    }

    psio_->close(PSIF_DFMP2_AIA, 1);
}
void RDFMP2::form_L() {
    // => Sizing <= //

    int nso = basisset_->nbf();
    int naocc = Caocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];
    int nia = Caocc_->colspi()[0] * Cavir_->colspi()[0];
    int naux = ribasis_->nbf();

    // => Thread Count <= //

    int num_threads = 1;
#ifdef _OPENMP
    num_threads = Process::environment.get_n_threads();
#endif

    // => Integrals <= //

    IntegralFactory rifactory(ribasis_, BasisSet::zero_ao_basis_set(), basisset_, basisset_);
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < num_threads; t++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory.eri()));
    }

    // => ERI Sieve <= //

    const auto& shell_pairs = eri[0]->shell_pairs();
    int npairs = shell_pairs.size();

    // => Memory Constraints <= //

    size_t memory = static_cast<size_t>((options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    memory -= static_cast<size_t>(naocc) * static_cast<size_t>(nso);
    memory -= static_cast<size_t>(navir) * static_cast<size_t>(nso);
    memory -= static_cast<size_t>(naocc) * static_cast<size_t>(navir);
    int max_rows;
    int maxP = ribasis_->max_function_per_shell();
    size_t row_cost = 0L;
    row_cost += static_cast<size_t>(nso)   * static_cast<size_t>(nso);
    row_cost += static_cast<size_t>(nso)   * static_cast<size_t>(naocc);
    row_cost += static_cast<size_t>(nso)   * static_cast<size_t>(navir);
    row_cost += static_cast<size_t>(naocc) * static_cast<size_t>(navir);
    size_t rows = memory / row_cost;
    rows = (rows > naux ? naux : rows);
    rows = (rows < maxP ? maxP : rows);
    max_rows = static_cast<int>(rows);

    // => Block Sizing <= //

    std::vector<int> Pstarts;
    int counter = 0;
    Pstarts.push_back(0);
    for (int P = 0; P < ribasis_->nshell(); P++) {
        int nP = ribasis_->shell(P).nfunction();
        if (counter + nP > max_rows) {
            counter = 0;
            Pstarts.push_back(P);
        }
        counter += nP;
    }
    Pstarts.push_back(ribasis_->nshell());
    // block_status(Pstarts, __FILE__,__LINE__);

    // => Temporary Buffers <= //

    auto Gia = std::make_shared<Matrix>("Gia", max_rows, naocc * navir);
    auto Gim = std::make_shared<Matrix>("Pim", max_rows, nso * naocc);
    auto Gam = std::make_shared<Matrix>("Pam", max_rows, nso * navir);
    auto Gmn = std::make_shared<Matrix>("Pmn", max_rows, nso * (size_t)nso);

    auto Giap = Gia->pointer();
    auto Gimp = Gim->pointer();
    auto Gamp = Gam->pointer();
    auto Gmnp = Gmn->pointer();

    auto Caoccp = Caocc_->pointer();
    auto Cavirp = Cavir_->pointer();

    std::vector<double> temp(naocc * navir);

    // => Targets <= //

    auto Lmi = std::make_shared<Matrix>("L_ma", nso, naocc);
    auto Lma = std::make_shared<Matrix>("L_mi", nso, navir);
    auto Lmip = Lmi->pointer();
    auto Lmap = Lma->pointer();

    // => PSIO <= //

    psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;

    // => Master Loop <= //

    for (int block = 0; block < Pstarts.size() - 1; block++) {
        // > Sizing < //

        int Pstart = Pstarts[block];
        int Pstop = Pstarts[block + 1];
        int NP = Pstop - Pstart;

        int pstart = ribasis_->shell(Pstart).function_index();
        int pstop = (Pstop == ribasis_->nshell() ? naux : ribasis_->shell(Pstop).function_index());
        int np = pstop - pstart;

        // > G_ia^P Read < //

        psio_->read(PSIF_DFMP2_AIA, "G(Q|ia)", (char*)Giap[0], sizeof(double) * np * nia, next_AIA, &next_AIA);

        // > Integrals < //
        Gmn->zero();
#pragma omp parallel for schedule(dynamic) num_threads(num_threads)
        for (long int PMN = 0L; PMN < static_cast<long int>(NP) * npairs; PMN++) {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            int P = PMN / npairs + Pstart;
            int MN = PMN % npairs;
            int M = shell_pairs[MN].first;
            int N = shell_pairs[MN].second;

            eri[thread]->compute_shell(P, 0, M, N);

            const double* buffer = eri[thread]->buffer();

            int nP = ribasis_->shell(P).nfunction();
            int oP = ribasis_->shell(P).function_index() - pstart;

            int nM = basisset_->shell(M).nfunction();
            int oM = basisset_->shell(M).function_index();

            int nN = basisset_->shell(N).nfunction();
            int oN = basisset_->shell(N).function_index();

            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {
                        Gmnp[p + oP][(m + oM) * nso + (n + oN)] = Gmnp[p + oP][(n + oN) * nso + (m + oM)] = (*buffer++);
                    }
                }
            }
        }

        // DEFINITION: L_ma := B(mi|Q)G(Q|ia)
        // Eq. 20, 29 of DiStasio
        // Used to construct Z-vector terms
        // This is spin-summed because G(Q|ia) is
        // N.B. Compared to DiStasio, this equation has an extra factor of -1. Adjust equations using this intermediate accordingly.

#pragma omp parallel for
        for (int p = 0; p < np; p++) {
            C_DGEMM('T', 'N', naocc, nso, nso, 1.0, Caoccp[0], naocc, Gmnp[p], nso, 0.0, Gimp[p], nso);
        }

        C_DGEMM('T', 'N', nso, navir, naocc * (size_t)np, 1.0, Gimp[0], nso, Giap[0], navir, 1.0, Lmap[0], navir);

        // Sort G_P^ia to G_P^ai
        for (int p = 0; p < np; p++) {
            ::memcpy((void*)temp.data(), (void*)Giap[p], sizeof(double) * naocc * navir);
            for (int i = 0; i < naocc; i++) {
                C_DCOPY(navir, &temp[i * navir], 1, &Giap[p][i], naocc);
            }
        }

        // DEFINITION: L_mi := B(ma|Q)G(Q|ia)
        // Eq. 19, 28 of DiStasio
        // Used to construct Z-vector terms
        // This is spin-summed because G(Q|ia) is

#pragma omp parallel for
        for (int p = 0; p < np; p++) {
            C_DGEMM('T', 'N', navir, nso, nso, 1.0, Cavirp[0], navir, Gmnp[p], nso, 0.0, Gamp[p], nso);
        }

        C_DGEMM('T', 'N', nso, naocc, navir * (size_t)np, 1.0, Gamp[0], nso, Giap[0], naocc, 1.0, Lmip[0], naocc);
    }

    psio_->write_entry(PSIF_DFMP2_AIA, "L_mi", (char*)Lmip[0], sizeof(double) * nso * naocc);
    psio_->write_entry(PSIF_DFMP2_AIA, "L_ma", (char*)Lmap[0], sizeof(double) * nso * navir);

    psio_->close(PSIF_DFMP2_AIA, 1);
}
void RDFMP2::form_P() {
    // => Sizing <= //

    int nso = basisset_->nbf();
    int nfocc = Cfocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];
    int naocc = Caocc_->colspi()[0];
    int nfvir = Cfvir_->colspi()[0];
    int nmo = nfocc + naocc + navir + nfvir;

    // => Tensors <= //

    auto Pij = std::make_shared<Matrix>("Pij", naocc, naocc);
    auto Pab = std::make_shared<Matrix>("Pab", navir, navir);
    auto PIj = std::make_shared<Matrix>("PIj", nfocc, naocc);
    auto PAb = std::make_shared<Matrix>("PAb", nfvir, navir);
    auto Ppq = std::make_shared<Matrix>("Ppq", nmo, nmo);

    auto Pijp = Pij->pointer();
    auto Pabp = Pab->pointer();
    auto PIjp = PIj->pointer();
    auto PAbp = PAb->pointer();
    auto Ppqp = Ppq->pointer();

    auto Lmi = std::make_shared<Matrix>("L_mi", nso, naocc);
    auto Lma = std::make_shared<Matrix>("L_ma", nso, navir);

    double** Lmip = Lmi->pointer();
    double** Lmap = Lma->pointer();

    // => Read-in <= //

    psio_->open(PSIF_DFMP2_AIA, 1);
    psio_->read_entry(PSIF_DFMP2_AIA, "P_ij", (char*)Pijp[0], sizeof(double) * naocc * naocc);
    psio_->read_entry(PSIF_DFMP2_AIA, "P_ab", (char*)Pabp[0], sizeof(double) * navir * navir);
    psio_->read_entry(PSIF_DFMP2_AIA, "L_mi", (char*)Lmip[0], sizeof(double) * nso * naocc);
    psio_->read_entry(PSIF_DFMP2_AIA, "L_ma", (char*)Lmap[0], sizeof(double) * nso * navir);

    // => Occ-Occ <= //

    for (int i = 0; i < naocc; i++) {
        ::memcpy((void*)&Ppqp[nfocc + i][nfocc], (void*)Pijp[i], sizeof(double) * naocc);
    }

    // => Virt-Virt <= //

    for (int a = 0; a < navir; a++) {
        ::memcpy((void*)&Ppqp[nfocc + naocc + a][nfocc + naocc], (void*)Pabp[a], sizeof(double) * navir);
    }

    // => Frozen-Core/Occ <= //

    if (nfocc) {
        // P_iJ := C_mJ L_mi / (ei - eJ)
        // This is spin-summed because L_mi is
        // This is a block of the OPDM that is added directly to the entire OPDM
        double** Cfoccp = Cfocc_->pointer();
        double* eps_foccp = eps_focc_->pointer();
        double* eps_aoccp = eps_aocc_->pointer();

        C_DGEMM('T', 'N', nfocc, naocc, nso, 1.0, Cfoccp[0], nfocc, Lmip[0], naocc, 0.0, PIjp[0], naocc);
        for (int i = 0; i < naocc; i++) {
            for (int J = 0; J < nfocc; J++) {
                PIjp[J][i] /= (eps_aoccp[i] - eps_foccp[J]);
            }
        }

        for (int J = 0; J < nfocc; J++) {
            C_DCOPY(naocc, PIjp[J], 1, &Ppqp[J][nfocc], 1);
            C_DCOPY(naocc, PIjp[J], 1, &Ppqp[nfocc][J], nmo);
        }
    }

    // => Frozen-Virt/Virt <= //

    if (nfvir) {
        // P_Ab := C_mA L_mb / (eb - eA)
        // This is spin-summed because L_mb is
        // This is a block of the OPDM that is added directly to the entire OPDM
        double** Cfvirp = Cfvir_->pointer();
        double* eps_fvirp = eps_fvir_->pointer();
        double* eps_avirp = eps_avir_->pointer();

        C_DGEMM('T', 'N', nfvir, navir, nso, 1.0, Cfvirp[0], nfvir, Lmap[0], navir, 0.0, PAbp[0], navir);
        for (int b = 0; b < navir; b++) {
            for (int A = 0; A < nfvir; A++) {
                PAbp[A][b] /= (eps_avirp[b] - eps_fvirp[A]);
            }
        }

        for (int B = 0; B < nfvir; B++) {
            C_DCOPY(navir, PAbp[B], 1, &Ppqp[nfocc + naocc + navir + B][nfocc + naocc], 1);
            C_DCOPY(navir, PAbp[B], 1, &Ppqp[nfocc + naocc][nfocc + naocc + navir + B], nmo);
        }
    }

    // DEFINITION: P_pq
    // This is the spin-summed density matrix, eq 6-10 of DiStasio, plus (later) the reference density.
    // N.B. This matrix is incomplete! The occupied/virtual block, eq. 10, cannot be added until we solve for the Z-vector.
    // Other formulations, differ in how they assign terms to the OPDM vs the Z-vector. The distinction isn't well-defined.
    // Further, we will add the reference to this.
    psio_->write_entry(PSIF_DFMP2_AIA, "P_pq", (char*)Ppqp[0], sizeof(double) * nmo * nmo);
    psio_->close(PSIF_DFMP2_AIA, 1);
}
void RDFMP2::form_W() {
    // => Sizing <= //

    int nso = basisset_->nbf();
    int nfocc = Cfocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];
    int naocc = Caocc_->colspi()[0];
    int nfvir = Cfvir_->colspi()[0];
    int nmo = nfocc + naocc + navir + nfvir;

    // => Tensors <= //

    auto Wpq1 = std::make_shared<Matrix>("Wpq1", nmo, nmo);
    auto Wpq1p = Wpq1->pointer();

    auto Ppq = std::make_shared<Matrix>("Ppq", nmo, nmo);
    auto Ppqp = Ppq->pointer();

    auto Lmi = std::make_shared<Matrix>("L_mi", nso, naocc);
    auto Lma = std::make_shared<Matrix>("L_ma", nso, navir);
    auto Lia = std::make_shared<Matrix>("L_ia", naocc + nfocc, navir + nfvir);

    auto Lmip = Lmi->pointer();
    auto Lmap = Lma->pointer();
    auto Liap = Lia->pointer();

    auto Cfoccp = Cfocc_->pointer();
    auto Caoccp = Caocc_->pointer();
    auto Cavirp = Cavir_->pointer();
    auto Cfvirp = Cfvir_->pointer();

    // => Read-in <= //

    psio_->open(PSIF_DFMP2_AIA, 1);
    psio_->read_entry(PSIF_DFMP2_AIA, "P_pq", (char*)Ppqp[0], sizeof(double) * nmo * nmo);
    psio_->read_entry(PSIF_DFMP2_AIA, "L_mi", (char*)Lmip[0], sizeof(double) * nso * naocc);
    psio_->read_entry(PSIF_DFMP2_AIA, "L_ma", (char*)Lmap[0], sizeof(double) * nso * navir);

    // => Term 1 <= //
    // All formulas for Term 1 are given in the block of five unnumbered equations on pg. 842
    // Note that the Lagrangians are spin-summed, so the final equations are as well.

    // We divide by 1/2 now and correct for it later. 
    // In some cases, we correct when we hermtivitize.
    // In other cases, we correct with a "subtle trick" later in the code. (Search those words to see where.)
    // In yet other cases, the DiStasio paper is just wrong, and no correction is needed.

    // => Occ/Occ <= //
    // W_ij := -1/2 * C_mi L_mj (DiStasio 11); Corrected when we hermitvitize.
    C_DGEMM('T', 'N', naocc, naocc, nso, -0.5, Caoccp[0], naocc, Lmip[0], naocc, 0.0, &Wpq1p[nfocc][nfocc], nmo);
    // => Frozen-Core/Occ <= //
    if (nfocc) {
        // W_Ij := -1/2 * C_mI L_mj (DiStasio 11)
        // No correction needed. Equation 11 needs an extra factor of 1/2.
        // You can compare DiStasio's formula to Evangelista's (correct) formula to see this yourself.
        C_DGEMM('T', 'N', nfocc, naocc, nso, -0.5, Cfoccp[0], nfocc, Lmip[0], naocc, 0.0, &Wpq1p[0][nfocc], nmo);
    }

    // => Virt/Virt <= //
    // W_ab := -1/2 * C_ma L_mb (DiStasio 12); Corrected when we hermitivitize
    C_DGEMM('T', 'N', navir, navir, nso, -0.5, Cavirp[0], navir, Lmap[0], navir, 0.0,
            &Wpq1p[nfocc + naocc][nfocc + naocc], nmo);
    // => Frozen-Virt/Virt <= //
    if (nfvir) {
        // W_Ab := -1/2 * C_mA L_mb (DiStasio 12)
        // No correction needed. Equation 12 needs an extra factor of 1/2.
        // You can compare DiStasio's formula to Evangelista's (correct) formula to see this yourself.
        C_DGEMM('T', 'N', nfvir, navir, nso, -0.5, Cfvirp[0], nfvir, Lmap[0], navir, 0.0,
                &Wpq1p[nfocc + naocc + navir][nfocc + naocc], nmo);
    }


    // > Occ-Virt <= //
    // WARNING! DiStasio 13 is in error. The last expression should be W_ia, and j must in occ.
    // W_ia := -1/2 C_mi L_ma (DiStasio 13); Corrected with the trick.
    C_DGEMM('T', 'N', naocc, navir, nso, -0.5, Caoccp[0], naocc, Lmap[0], navir, 0.0, &Wpq1p[nfocc][nfocc + naocc],
            nmo);
    if (nfocc) {
        // W_Ia := -1/2 C_mI L_ma (DiStasio 13); Corrected with the trick
        C_DGEMM('T', 'N', nfocc, navir, nso, -0.5, Cfoccp[0], nfocc, Lmap[0], navir, 0.0, &Wpq1p[0][nfocc + naocc],
                nmo);
    }

    // > Vir-Occ < //
    // These terms don't actually belong in the W terms, but we need them in the L terms.
    // We'll keep it here for now and get rid of it later with the trick.
    C_DGEMM('T', 'N', navir, naocc, nso, -0.5, Cavirp[0], navir, Lmip[0], naocc, 0.0, &Wpq1p[nfocc + naocc][nfocc],
            nmo);
    if (nfvir) {
        C_DGEMM('T', 'N', nfvir, naocc, nso, -0.5, Cfvirp[0], nfvir, Lmip[0], naocc, 0.0,
                &Wpq1p[nfocc + naocc + navir][nfocc], nmo);
    }

    // => Lia (L contributions) <= //
    // First two terms of DiStasio 18... With some modfications

    // Multiply by 2 to remove factor of 0.5 applied above
    for (int i = 0; i < (nfocc + naocc); i++) {
        for (int a = 0; a < (nfvir + navir); a++) {
            Liap[i][a] = 2.0 * (Wpq1p[i][a + naocc + nfocc] - Wpq1p[a + naocc + nfocc][i]);
        }
    }

    // W_pq = W_pq + W_qp
    Wpq1->hermitivitize();
    Wpq1->scale(2.0);

    // => Write-out <= //

    psio_->write_entry(PSIF_DFMP2_AIA, "L_ia", (char*)Liap[0], sizeof(double) * (naocc + nfocc) * (navir + nfvir));
    psio_->write_entry(PSIF_DFMP2_AIA, "W", (char*)Wpq1p[0], sizeof(double) * nmo * nmo);
    psio_->close(PSIF_DFMP2_AIA, 1);
}
void RDFMP2::form_Z() {
    // => Sizing <= //

    int nso = basisset_->nbf();
    int nfocc = Cfocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];
    int naocc = Caocc_->colspi()[0];
    int nfvir = Cfvir_->colspi()[0];
    int nmo = nfocc + naocc + navir + nfvir;
    int nocc = nfocc + naocc;
    int nvir = nfvir + navir;

    // => Tensors <= //

    auto Wpq1 = std::make_shared<Matrix>("Wpq1", nmo, nmo);
    auto Wpq1p = Wpq1->pointer();
    auto Wpq2 = std::make_shared<Matrix>("Wpq2", nmo, nmo);
    auto Wpq2p = Wpq2->pointer();
    auto Wpq3 = std::make_shared<Matrix>("Wpq3", nmo, nmo);
    auto Wpq3p = Wpq3->pointer();

    auto Ppq = std::make_shared<Matrix>("Ppq", nmo, nmo);
    auto Ppqp = Ppq->pointer();

    auto Lia = std::make_shared<Matrix>("L_ia", naocc + nfocc, navir + nfvir);
    auto Liap = Lia->pointer();

    auto Cocc = Ca_subset("AO", "OCC");
    auto Cvir = Ca_subset("AO", "VIR");
    auto C = Ca_subset("AO", "ALL");

    auto Cfoccp = Cfocc_->pointer();
    auto Caoccp = Caocc_->pointer();
    auto Cavirp = Cavir_->pointer();
    auto Cfvirp = Cfvir_->pointer();

    auto Coccp = Cocc->pointer();
    auto Cvirp = Cvir->pointer();
    auto Cp = C->pointer();

    auto eps_occ = epsilon_a_subset("AO", "OCC");
    auto eps_vir = epsilon_a_subset("AO", "VIR");
    auto eps = epsilon_a_subset("AO", "ALL");
    auto epsp = eps->pointer();

    // => CPHF/JK Object <= //

    auto cphf = std::make_shared<RCPHF>(reference_wavefunction_, options_);
    cphf->set_C(C);
    cphf->set_Caocc(Cocc);
    cphf->set_Cavir(Cvir);
    cphf->set_eps_aocc(eps_occ);
    cphf->set_eps_avir(eps_vir);
    cphf->preiterations();

    auto jk = cphf->jk();
    auto& Cl = jk->C_left();
    auto& Cr = jk->C_right();
    const auto& J = jk->J();
    const auto& K = jk->K();

    // => Read-in <= //

    psio_->open(PSIF_DFMP2_AIA, 1);
    psio_->read_entry(PSIF_DFMP2_AIA, "P_pq", (char*)Ppqp[0], sizeof(double) * nmo * nmo);

    auto T = std::make_shared<Matrix>("T", nocc, nso);
    auto Tp = T->pointer();
    auto dPpq = std::make_shared<Matrix>("dP", nmo, nmo);
    auto dPpqp = dPpq->pointer();
    auto AP = std::make_shared<Matrix>("A_mn^ls P_ls^(2)", nso, nso);
    auto APp = AP->pointer();
    SharedMatrix Dtemp;

    if (options_.get_bool("OPDM_RELAX")) {
        psio_->read_entry(PSIF_DFMP2_AIA, "W", (char*)Wpq1p[0], sizeof(double) * nmo * nmo);
        psio_->read_entry(PSIF_DFMP2_AIA, "L_ia", (char*)Liap[0], sizeof(double) * (naocc + nfocc) * (navir + nfvir));

        // => Lia += 1/2 A_pqia P_pq (unrelaxed) <= //

        // We're using a Cholesky decomposition so we can use JK. See README.
        // > Factor the unrelaxed P^(2) (hopefully low rank) < //
        std::pair<SharedMatrix, SharedMatrix> factor1 =
            Ppq->partial_square_root(options_.get_double("DFMP2_P2_TOLERANCE"));
        auto P1 = factor1.first;
        auto N1 = factor1.second;

        // > Back-transform the transition orbitals < //
        auto P1AO = linalg::doublet(C, P1);
        auto N1AO = linalg::doublet(C, N1);

        // > Form the J/K-like matrices (P,N contributions are separable) < //
        Cl.clear();
        Cr.clear();
        Cl.push_back(P1AO);
        Cl.push_back(N1AO);

        jk->compute();

        SharedMatrix J1 = J[0];
        SharedMatrix K1 = K[0];
        J1->subtract(J[1]);
        K1->subtract(K[1]);
        auto J1p = J1->pointer();
        auto K1p = K1->pointer();

        J1->scale(2.0);
        K1->scale(1.0);
        AP->add(J1);
        AP->subtract(K1);

        // > Form the contribution to Lia from the J/K-like matrices < //

        // L_ia += -1.0 (spin) C_mi { [ 4(mn|pq) - (mq|pn) - (mq|nq)] P_pq } C_na
        C_DGEMM('T', 'N', nocc, nso, nso, 1.0, Coccp[0], nocc, APp[0], nso, 0.0, Tp[0], nso);
        C_DGEMM('N', 'N', nocc, nvir, nso, 1.0, Tp[0], nso, Cvirp[0], nvir, 1.0, Liap[0], nvir);

        // Lia->print();

        // => (\delta_ij \delta_ab (\epsilon_a - \epsilon_i) + A_ia,jb) Z_jb = L_ia <= //

        auto& b = cphf->b();
        auto& x = cphf->x();

        // The occ/vir (frozen or not) blocks of the density are given by b.
        // See Aikens eq. 67
        b["L_ia"] = Lia;
        cphf->compute_energy();
        auto Zia = x["L_ia"];
        Zia->scale(-1.0);

        // > Add Pia and Pai into the OPDM < //
        auto Ziap = Zia->pointer();
        for (int i = 0; i < nocc; i++) {
            for (int a = 0; a < nvir; a++) {
                dPpqp[i][a + nocc] = dPpqp[a + nocc][i] = Ziap[i][a];
            }
        }

        // UPDATE: P_pq: This quantity is now the relaxed density matrix. The occupied-virtual block is nonzero.
        Ppq->add(dPpq);

        Dtemp = Ppq->clone();

        // Scale the OPDM so it's no longer spin-summed.
        Dtemp->scale(0.5);

        // Add in the reference density
        for (int i = 0; i < nocc; ++i) Dtemp->add(i, i, 1.0);
        Ca_ = std::make_shared<Matrix>("DF-MP2 Natural Orbitals", nsopi_, nmopi_);
        epsilon_a_ = std::make_shared<Vector>("DF-MP2 NO Occupations", nmopi_);
        Da_ = std::make_shared<Matrix>("DF-MP2 relaxed density", nsopi_, nsopi_);

    } else {
        // Don't relax the OPDM
        Dtemp = Ppq->clone();

        // Scale the OPDM so it's no longer spin-summed.
        Dtemp->scale(0.5);

        // Add in the reference density
        for (int i = 0; i < nocc; ++i) Dtemp->add(i, i, 1.0);
        Ca_ = std::make_shared<Matrix>("DF-MP2 (unrelaxed) Natural Orbitals", nsopi_, nmopi_);
        epsilon_a_ = std::make_shared<Vector>("DF-MP2 (unrelaxed) NO Occupations", nmopi_);
        Da_ = std::make_shared<Matrix>("DF-MP2 unrelaxed density", nsopi_, nsopi_);
    }

    compute_opdm_and_nos(Dtemp, Da_, Ca_, epsilon_a_);
    Cb_ = Ca_;
    epsilon_b_ = epsilon_a_;
    Db_ = Da_;

    if (options_.get_bool("ONEPDM")) {
        // Shut everything down; only the OPDM was requested
        psio_->write_entry(PSIF_DFMP2_AIA, "P_pq", (char*)Ppqp[0], sizeof(double) * nmo * nmo);
        psio_->close(PSIF_DFMP2_AIA, 1);
        cphf->postiterations();

        return;
    }

    // => Wik -= 1/2 A_pqik P_pq (relaxed) <= //

    // We're using a Cholesky decomposition so we can use JK. See README.
    // > Factor the relaxation contribution dP^(2) (hopefully low rank) < //
    std::pair<SharedMatrix, SharedMatrix> factor2 =
        dPpq->partial_square_root(options_.get_double("DFMP2_P2_TOLERANCE"));
    auto P2 = factor2.first;
    auto N2 = factor2.second;

    // > Back-transform the transition orbitals < //
    auto P2AO = linalg::doublet(C, P2);
    auto N2AO = linalg::doublet(C, N2);

    // > Form the J/K-like matrices (P,N contributions are separable) < //
    Cl.clear();
    Cr.clear();
    Cl.push_back(P2AO);
    Cl.push_back(N2AO);

    jk->compute();

    auto J2 = J[0];
    auto K2 = K[0];
    J2->subtract(J[1]);
    K2->subtract(K[1]);
    auto J2p = J2->pointer();
    auto K2p = K2->pointer();

    J2->scale(2.0);
    K2->scale(1.0);
    AP->add(J2);
    AP->subtract(K2);

    // > Form the contribution to Lia from the J/K-like matrices < //
    // Form the spin-summed version of Distasio 17. It takes a spin-summed density and
    // returns a spin-summed W matrix.
    // Two important notes:
    // 1. Because P is symmetric, (mq|pn) P_pq = (mp|nq) P_pq, so we combine the last two terms.
    // 2. We absorb the 1/2 into the sum of integral, so 4J - 2K => 2J - K.

    // W_ik -= 1/2 (spin) C_mi { [ 4(mn|pq) - (mq|pn) - (mp|nq)] P_pq } C_nk
    C_DGEMM('T', 'N', nocc, nso, nso, 1.0, Coccp[0], nocc, APp[0], nso, 0.0, Tp[0], nso);

    // occ-occ term
    C_DGEMM('N', 'N', nocc, nocc, nso, -1.0, Tp[0], nso, Coccp[0], nocc, 0.0, &Wpq3p[0][0], nmo);

    // This next term is an extremely subtle trick. It's helpful to compare Wang to DiStasio.
    // All orbital energies correct the orbital energy term below. All other terms correct for
    // something back in form_W.
    // WANG TO DISTASIO CONVERSIONS
    // w^ij_ab => t^ij_ab
    // w^ab_ij => - t^ij_ab
    // z^P_Q => 2 P_PQ
    // v^pq_rs => g^pq_rs (Antisymmetrized integral, always valid)
    // v^qr_ps => 1/2 A_pq,rs (when contracted against a density; special case of above)
    // EXPLANATION OF THE CASES
    // Case 1: I frozen, A frozen
    // z^I_A (e_I - e_A) = v_qI,pA z_pq
    // 2 P^I_A (eI - eA) = A_pq,IA P_pq
    // P^I_A (eI - eA) = 1/2 A_pq,IA P_pq
    // We have 1/2 A_pq,IA P_pq, so end with -1/2 P^I_A(eI-eA).
    // Case 2: I frozen, A active
    // z^I_a (e_I - e_a) + v^Ib_ij t^ij_ab = v^qI_pa z_pq
    // 2 P^I_a (e_I - e_a) + g^Ib_ij t^ij_ab = A_pq,Ia P_pq
    // P^I_a (e_I - e_a) + 1/2 g^Ib_ij t^ij_ab = 1/2 A_pq,Ia P_pq
    // We have 1/2 A_pq,Ia P_pq, so end with -1/2 P^I_a (e_I - e_a) - 1/4 g^Ib_ij t^ij_ab.
    // Case 3: i active, A active
    // z^i_a (e_i - e_a) + v^ib_kl t^kl_ab - v^cd_aj t^ij_cd = v^qi_pa z^p_q
    // 2P^i_a (e_i - e_a) + g^ib_kl t^kl_ab - g^cd_aj t^ij_cd = A_pq,ia P_pq
    // P^i_a (e_i - e_a) + 1/2 g^ib_kl t^kl_ab - 1/2 g^cd_aj t^ij_cd = 1/2 A_pq,ia P_pq
    // We have 1/2 A_pq,ia P_pq, so end with -1/2 P^i_a (e_i - e_a) -1/4 g^ib_kl t^kl_ab + 1/4 g^cd_aj t^ij_cd
    // Case 4: i active, A frozen
    // z^i_A (e_i - e_a) - t^ij_ab v^ab_Aj = v^qi_pA z_pq
    // 2 P^i_A (e_i - e_a) - t^ij_ab g^ab_Aj = A_pq,iA P_pq
    // P^i_A (e_i - e_a) - 1/2 t^ij_ab g^ab_Aj = 1/2 A_pq,iA P_pq
    // We have 1/2 A_pq,iA P_pq, so end with -1/2 P^i_A (e_i - e_a) + 1/4 t^ij_ab g^ab_Aj
    C_DGEMM('N', 'N', nocc, nvir, nso, -0.5, Tp[0], nso, Cvirp[0], nvir, 0.0, &Wpq3p[0][nocc], nmo);
    C_DGEMM('T', 'T', nvir, nocc, nso, -0.5, Cvirp[0], nvir, Tp[0], nso, 0.0, &Wpq3p[nocc][0], nmo);

    // => W Term 2 <= //
    // Equations 14-16 of DiStasio
    // Although the summation range includes some terms that don't appear in Term 2, they all have zero Ppq.
    // The bigger concern is that for the occupied-virtual block, the correct formula should have a prefactor
    // of -1 and exclude the virtual orbital. This is corrected in the above terms.
    // P_pq is spin-summed, handling spin-integration for us.

    for (int p = 0; p < nmo; p++) {
        for (int q = 0; q < nmo; q++) {
            Wpq2p[p][q] = -0.5 * (epsp[p] + epsp[q]) * Ppqp[p][q];
        }
    }

    // => Final W <= //

    Wpq1->add(Wpq2);
    Wpq1->add(Wpq3);
    Wpq1->set_name("Wpq");

    psio_->write_entry(PSIF_DFMP2_AIA, "W", (char*)Wpq1p[0], sizeof(double) * nmo * nmo);

    // => Final P <= //

    psio_->write_entry(PSIF_DFMP2_AIA, "P_pq", (char*)Ppqp[0], sizeof(double) * nmo * nmo);

    // => Finalize <= //

    psio_->close(PSIF_DFMP2_AIA, 1);
    cphf->postiterations();
}
void RDFMP2::form_gradient() {
    // => Sizing <= //

    int nso = basisset_->nbf();
    int nfocc = Cfocc_->colspi()[0];
    int naocc = Caocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];
    int nfvir = Cfvir_->colspi()[0];
    int nmo = nfocc + naocc + navir + nfvir;
    int nocc = nfocc + naocc;
    int nvir = nfvir + navir;

    // => Tensors <= //

    auto W = std::make_shared<Matrix>("W", nmo, nmo);
    auto Wp = W->pointer();

    auto P2 = std::make_shared<Matrix>("P_pq", nmo, nmo);
    auto P2p = P2->pointer();

    auto Cocc = reference_wavefunction_->Ca_subset("AO", "OCC");
    auto C = reference_wavefunction_->Ca_subset("AO", "ALL");

    auto Coccp = Cocc->pointer();
    auto Cp = C->pointer();

    auto eps = reference_wavefunction_->epsilon_a_subset("AO", "ALL");
    auto epsp = eps->pointer();

    // => Read-in <= //

    psio_->open(PSIF_DFMP2_AIA, 1);
    psio_->read_entry(PSIF_DFMP2_AIA, "P_pq", (char*)P2p[0], sizeof(double) * nmo * nmo);
    psio_->read_entry(PSIF_DFMP2_AIA, "W", (char*)Wp[0], sizeof(double) * nmo * nmo);

    // => Dress for SCF <= //

    SharedMatrix P2F(P2->clone());
    double** P2Fp = P2F->pointer();
    P2F->scale(2.0);

    // UPDATE: W : The HF EWDM is added in.
    // UPDATE: P_pq : The spin-summed RHF density matrix is added to P_pq.
    W->scale(-1.0);
    for (int i = 0; i < nocc; i++) {
        Wp[i][i] += 2.0 * epsp[i];
        P2p[i][i] += 2.0;
        P2Fp[i][i] += 2.0;
    }

    psio_->write_entry(PSIF_DFMP2_AIA, "P_pq", (char*)P2p[0], sizeof(double) * nmo * nmo);
    psio_->write_entry(PSIF_DFMP2_AIA, "W", (char*)Wp[0], sizeof(double) * nmo * nmo);

    // => Factorize the P matrix <= //

    P2F->scale(0.5);
    auto factor = P2F->partial_square_root(options_.get_double("DFMP2_P_TOLERANCE"));
    P2F->scale(2.0);

    auto P1 = factor.first;
    auto N1 = factor.second;
    auto P1p = P1->pointer();
    auto N1p = N1->pointer();

    // => Back-transform <= //

    auto T1 = std::make_shared<Matrix>("T", nmo, nso);
    auto PAO = std::make_shared<Matrix>("P AO", nso, nso);
    auto PFAO = std::make_shared<Matrix>("PF AO", nso, nso);
    auto WAO = std::make_shared<Matrix>("W AO", nso, nso);
    auto P1AO = std::make_shared<Matrix>("P1 AO", nso, P1->colspi()[0]);
    auto N1AO = std::make_shared<Matrix>("N1 AO", nso, N1->colspi()[0]);

    auto T1p = T1->pointer();
    auto PAOp = PAO->pointer();
    auto PFAOp = PFAO->pointer();
    auto WAOp = WAO->pointer();
    auto P1AOp = P1AO->pointer();
    auto N1AOp = N1AO->pointer();

    C_DGEMM('N', 'T', nmo, nso, nmo, 1.0, P2p[0], nmo, Cp[0], nmo, 0.0, T1p[0], nso);
    C_DGEMM('N', 'N', nso, nso, nmo, 1.0, Cp[0], nmo, T1p[0], nso, 0.0, PAOp[0], nso);

    C_DGEMM('N', 'T', nmo, nso, nmo, 1.0, P2Fp[0], nmo, Cp[0], nmo, 0.0, T1p[0], nso);
    C_DGEMM('N', 'N', nso, nso, nmo, 1.0, Cp[0], nmo, T1p[0], nso, 0.0, PFAOp[0], nso);

    C_DGEMM('N', 'T', nmo, nso, nmo, 1.0, Wp[0], nmo, Cp[0], nmo, 0.0, T1p[0], nso);
    C_DGEMM('N', 'N', nso, nso, nmo, 1.0, Cp[0], nmo, T1p[0], nso, 0.0, WAOp[0], nso);

    if (P1->colspi()[0]) {
        C_DGEMM('N', 'N', nso, P1->colspi()[0], nmo, 1.0, Cp[0], nmo, P1p[0], P1->colspi()[0], 0.0, P1AOp[0],
                P1->colspi()[0]);
    }

    if (N1->colspi()[0]) {
        C_DGEMM('N', 'N', nso, N1->colspi()[0], nmo, 1.0, Cp[0], nmo, N1p[0], N1->colspi()[0], 0.0, N1AOp[0],
                N1->colspi()[0]);
    }

    // => Random extra drivers for the JK gradients <= //

    SharedMatrix D(PAO->clone());
    auto Dp = D->pointer();
    C_DGEMM('N', 'T', nso, nso, nocc, 2.0, Coccp[0], nocc, Coccp[0], nocc, 0.0, Dp[0], nso);

    SharedMatrix Ds(D->clone());
    Ds->scale(0.5);

    SharedMatrix PFAOs(PFAO->clone());
    PFAOs->scale(0.5);

    auto mints = std::make_shared<MintsHelper>(basisset_, options_);

    // => Gogo Gradients <= //

    std::vector<std::string> gradient_terms;
    gradient_terms.push_back("Nuclear");
    gradient_terms.push_back("Core");
    gradient_terms.push_back("Overlap");
    gradient_terms.push_back("Coulomb");
    gradient_terms.push_back("Exchange");
    gradient_terms.push_back("Correlation");
    gradient_terms.push_back("Total");

    // => Sizings <= //
    auto natom = molecule_->natom();

    // => Nuclear Gradient <= //
    gradients_["Nuclear"] = SharedMatrix(molecule_->nuclear_repulsion_energy_deriv1(dipole_field_strength_).clone());
    gradients_["Nuclear"]->set_name("Nuclear Gradient");

    // => Kinetic Gradient <= //
    timer_on("Grad: V T Perturb");
    gradients_["Core"] = mints->core_hamiltonian_grad(PAO);
    timer_off("Grad: V T Perturb");

    // If an external field exists, add it to the one-electron Hamiltonian
    if (external_pot_) {
        gradient_terms.push_back("External Potential");
        timer_on("Grad: External");
        gradients_["External Potential"] = external_pot_->computePotentialGradients(basisset_, PAO);
        timer_off("Grad: External");
    }  // end external

    // => Overlap Gradient <= //
    timer_on("Grad: S");
    gradients_["Overlap"] = mints->overlap_grad(WAO);
    gradients_["Overlap"]->scale(-1.0);
    timer_off("Grad: S");

    // => Two-Electron Gradient <= //

    timer_on("Grad: JK");

    auto jk = CorrGrad::build_CorrGrad(mintshelper_);
    jk->set_memory((size_t)(options_.get_double("SCF_MEM_SAFETY_FACTOR") * memory_ / 8L));

    jk->set_Ca(Cocc);
    jk->set_Cb(Cocc);
    jk->set_La(P1AO);
    jk->set_Lb(P1AO);
    jk->set_Ra(N1AO);
    jk->set_Rb(N1AO);
    jk->set_Da(Ds);
    jk->set_Db(Ds);
    jk->set_Dt(D);
    jk->set_Pa(PFAOs);
    jk->set_Pb(PFAOs);
    jk->set_Pt(PFAO);

    jk->print_header();
    jk->compute_gradient();

    auto& jk_gradients = jk->gradients();
    gradients_["Coulomb"] = jk_gradients["Coulomb"];
    gradients_["Exchange"] = jk_gradients["Exchange"];
    gradients_["Exchange"]->scale(-1.0);

    timer_off("Grad: JK");

    // => Correlation Gradient (Previously computed) <= //

    auto correlation = SharedMatrix(gradients_["Nuclear"]->clone());
    correlation->zero();
    correlation->add(gradients_["(A|mn)^x"]);
    correlation->add(gradients_["(A|B)^x"]);
    gradients_["Correlation"] = correlation;
    gradients_["Correlation"]->set_name("Correlation Gradient");

    // => Total Gradient <= //
    auto total = SharedMatrix(gradients_["Nuclear"]->clone());
    total->zero();

    for (int i = 0; i < gradient_terms.size(); i++) {
        if (gradients_.count(gradient_terms[i])) {
            total->add(gradients_[gradient_terms[i]]);
        }
    }

    gradients_["Total"] = total;
    gradients_["Total"]->set_name("Total Gradient");

    // => Finalize <= //

    psio_->close(PSIF_DFMP2_AIA, 1);
}

UDFMP2::UDFMP2(SharedWavefunction ref_wfn, Options& options, std::shared_ptr<PSIO> psio)
    : DFMP2(ref_wfn, options, psio) {
    common_init();
}
UDFMP2::~UDFMP2() {}
void UDFMP2::common_init() {
    Caocc_a_ = Ca_subset("AO", "ACTIVE_OCC");
    Cavir_a_ = Ca_subset("AO", "ACTIVE_VIR");
    Caocc_b_ = Cb_subset("AO", "ACTIVE_OCC");
    Cavir_b_ = Cb_subset("AO", "ACTIVE_VIR");

    eps_aocc_a_ = epsilon_a_subset("AO", "ACTIVE_OCC");
    eps_avir_a_ = epsilon_a_subset("AO", "ACTIVE_VIR");
    eps_aocc_b_ = epsilon_b_subset("AO", "ACTIVE_OCC");
    eps_avir_b_ = epsilon_b_subset("AO", "ACTIVE_VIR");
}
void UDFMP2::print_header() {
    int nthread = 1;
#ifdef _OPENMP
    nthread = Process::environment.get_n_threads();
#endif

    outfile->Printf("\t --------------------------------------------------------\n");
    outfile->Printf("\t                          DF-MP2                         \n");
    outfile->Printf("\t      2nd-Order Density-Fitted Moller-Plesset Theory     \n");
    outfile->Printf("\t              UMP2 Wavefunction, %3d Threads             \n", nthread);
    outfile->Printf("\t                                                         \n");
    outfile->Printf("\t        Rob Parrish, Justin Turney, Andy Simmonett,      \n");
    outfile->Printf("\t           Ed Hohenstein, and C. David Sherrill          \n");
    outfile->Printf("\t --------------------------------------------------------\n");
    outfile->Printf("\n");

    int focc_a = frzcpi_.sum();
    int fvir_a = frzvpi_.sum();
    int aocc_a = Caocc_a_->colspi()[0];
    int avir_a = Cavir_a_->colspi()[0];
    int occ_a = focc_a + aocc_a;
    int vir_a = fvir_a + avir_a;

    int focc_b = frzcpi_.sum();
    int fvir_b = frzvpi_.sum();
    int aocc_b = Caocc_b_->colspi()[0];
    int avir_b = Cavir_b_->colspi()[0];
    int occ_b = focc_b + aocc_b;
    int vir_b = fvir_b + avir_b;

    outfile->Printf("\t --------------------------------------------------------\n");
    outfile->Printf("\t                 NBF = %5d, NAUX = %5d\n", basisset_->nbf(), ribasis_->nbf());
    outfile->Printf("\t --------------------------------------------------------\n");
    outfile->Printf("\t %7s %7s %7s %7s %7s %7s %7s\n", "CLASS", "FOCC", "OCC", "AOCC", "AVIR", "VIR", "FVIR");
    outfile->Printf("\t %7s %7d %7d %7d %7d %7d %7d\n", "ALPHA", focc_a, occ_a, aocc_a, avir_a, vir_a, fvir_a);
    outfile->Printf("\t %7s %7d %7d %7d %7d %7d %7d\n", "BETA", focc_b, occ_b, aocc_b, avir_b, vir_b, fvir_b);
    outfile->Printf("\t --------------------------------------------------------\n\n");
}
void UDFMP2::form_Aia() {
    // ERI objects
    int nthread = 1;
#ifdef _OPENMP
    if (options_.get_int("DF_INTS_NUM_THREADS") == 0) {
        nthread = Process::environment.get_n_threads();
    } else {
        nthread = options_.get_int("DF_INTS_NUM_THREADS");
    }
#endif

    IntegralFactory factory(ribasis_, BasisSet::zero_ao_basis_set(), basisset_, basisset_);
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    std::vector<const double*> buffer;
    for (int thread = 0; thread < nthread; thread++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(factory.eri()));
        buffer.push_back(eri[thread]->buffer());
    }

    // Schwarz Sieve
    const auto& shell_pairs = eri[0]->shell_pairs();
    const size_t npairs = shell_pairs.size();

    // Sizing
    int nso = basisset_->nbf();
    int naux = ribasis_->nbf();
    int naocc_a = Caocc_a_->colspi()[0];
    int navir_a = Cavir_a_->colspi()[0];
    int naocc_b = Caocc_b_->colspi()[0];
    int navir_b = Cavir_b_->colspi()[0];
    int naocc = (naocc_a > naocc_b ? naocc_a : naocc_b);
    int navir = (navir_a > navir_b ? navir_a : navir_b);
    int maxQ = ribasis_->max_function_per_shell();

    // Max block size in naux
    size_t Amn_cost_per_row = nso * (size_t)nso;
    size_t Ami_cost_per_row = nso * (size_t)naocc;
    size_t Aia_cost_per_row = naocc * (size_t)navir;
    size_t total_cost_per_row = Amn_cost_per_row + Ami_cost_per_row + Aia_cost_per_row;
    size_t doubles = ((size_t)(options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    size_t max_temp = doubles / (total_cost_per_row);
    int max_naux = (max_temp > (size_t)naux ? naux : max_temp);
    max_naux = (max_naux < maxQ ? maxQ : max_naux);

    // Block extents
    std::vector<int> block_Q_starts;
    int counter = 0;
    block_Q_starts.push_back(0);
    for (int Q = 0; Q < ribasis_->nshell(); Q++) {
        int nQ = ribasis_->shell(Q).nfunction();
        if (counter + nQ > max_naux) {
            counter = 0;
            block_Q_starts.push_back(Q);
        }
        counter += nQ;
    }
    block_Q_starts.push_back(ribasis_->nshell());
    // block_status(block_Q_starts, __FILE__,__LINE__);

    // Tensor blocks
    auto Amn = std::make_shared<Matrix>("(A|mn) Block", max_naux, nso * (size_t)nso);
    auto Ami = std::make_shared<Matrix>("(A|mi) Block", max_naux, nso * (size_t)naocc);
    auto Aia = std::make_shared<Matrix>("(A|ia) Block", max_naux, naocc * (size_t)navir);
    auto Amnp = Amn->pointer();
    auto Amip = Ami->pointer();
    auto Aiap = Aia->pointer();

    // C Matrices
    auto Caoccap = Caocc_a_->pointer();
    auto Cavirap = Cavir_a_->pointer();
    auto Caoccbp = Caocc_b_->pointer();
    auto Cavirbp = Cavir_b_->pointer();

    psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_NEW);
    psio_address next_AIA = PSIO_ZERO;
    psio_->open(PSIF_DFMP2_QIA, PSIO_OPEN_NEW);
    psio_address next_QIA = PSIO_ZERO;

    // Loop over blocks of Qshell
    for (int block = 0; block < block_Q_starts.size() - 1; block++) {
        // Block sizing/offsets
        int Qstart = block_Q_starts[block];
        int Qstop = block_Q_starts[block + 1];
        int qoff = ribasis_->shell(Qstart).function_index();
        int nrows = (Qstop == ribasis_->nshell()
                         ? ribasis_->nbf() - ribasis_->shell(Qstart).function_index()
                         : ribasis_->shell(Qstop).function_index() - ribasis_->shell(Qstart).function_index());

        // Clear Amn for Schwarz sieve
        ::memset((void*)Amnp[0], '\0', sizeof(double) * nrows * nso * nso);

        // Compute TEI tensor block (A|mn)
        timer_on("DFMP2 (A|mn)");
#pragma omp parallel for schedule(dynamic) num_threads(nthread)
        for (long int QMN = 0L; QMN < (Qstop - Qstart) * (size_t)npairs; QMN++) {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            int Q = QMN / npairs + Qstart;
            int MN = QMN % npairs;

            auto pair = shell_pairs[MN];
            int M = pair.first;
            int N = pair.second;

            int nq = ribasis_->shell(Q).nfunction();
            int nm = basisset_->shell(M).nfunction();
            int nn = basisset_->shell(N).nfunction();

            int sq = ribasis_->shell(Q).function_index();
            int sm = basisset_->shell(M).function_index();
            int sn = basisset_->shell(N).function_index();

            eri[thread]->compute_shell(Q, 0, M, N);
            buffer[thread] = eri[thread]->buffer();

            for (int oq = 0; oq < nq; oq++) {
                for (int om = 0; om < nm; om++) {
                    for (int on = 0; on < nn; on++) {
                        Amnp[sq + oq - qoff][(om + sm) * nso + (on + sn)] =
                            Amnp[sq + oq - qoff][(on + sn) * nso + (om + sm)] =
                                buffer[thread][oq * nm * nn + om * nn + on];
                    }
                }
            }
        }
        timer_off("DFMP2 (A|mn)");

        // => Alpha Case <= //

        // Compute (A|mi) tensor block (A|mn) C_ni
        timer_on("DFMP2 (A|mn)C_mi");
        C_DGEMM('N', 'N', nrows * (size_t)nso, naocc_a, nso, 1.0, Amnp[0], nso, Caoccap[0], naocc_a, 0.0, Amip[0],
                naocc);
        timer_off("DFMP2 (A|mn)C_mi");

        // Compute (A|ia) tensor block (A|ia) = (A|mi) C_ma
        timer_on("DFMP2 (A|mi)C_na");
#pragma omp parallel for
        for (int row = 0; row < nrows; row++) {
            C_DGEMM('T', 'N', naocc_a, navir_a, nso, 1.0, Amip[row], naocc, Cavirap[0], navir_a, 0.0,
                    &Aiap[0][row * (size_t)naocc_a * navir_a], navir_a);
        }
        timer_off("DFMP2 (A|mi)C_na");

        // Stripe (A|ia) out to disk
        timer_on("DFMP2 Aia Write");
        psio_->write(PSIF_DFMP2_AIA, "A(Q|ia)", (char*)Aiap[0], sizeof(double) * nrows * naocc_a * navir_a, next_AIA,
                     &next_AIA);
        timer_off("DFMP2 Aia Write");

        // => Beta Case <= //

        // Compute (A|mi) tensor block (A|mn) C_ni
        timer_on("DFMP2 (A|mn)C_mi");
        C_DGEMM('N', 'N', nrows * (size_t)nso, naocc_b, nso, 1.0, Amnp[0], nso, Caoccbp[0], naocc_b, 0.0, Amip[0],
                naocc);
        timer_off("DFMP2 (A|mn)C_mi");

        // Compute (A|ia) tensor block (A|ia) = (A|mi) C_ma
        timer_on("DFMP2 (A|mi)C_na");
#pragma omp parallel for
        for (int row = 0; row < nrows; row++) {
            C_DGEMM('T', 'N', naocc_b, navir_b, nso, 1.0, Amip[row], naocc, Cavirbp[0], navir_b, 0.0,
                    &Aiap[0][row * (size_t)naocc_b * navir_b], navir_b);
        }
        timer_off("DFMP2 (A|mi)C_na");

        // Stripe (A|ia) out to disk
        timer_on("DFMP2 Aia Write");
        psio_->write(PSIF_DFMP2_QIA, "A(Q|ia)", (char*)Aiap[0], sizeof(double) * nrows * naocc_b * navir_b, next_QIA,
                     &next_QIA);
        timer_off("DFMP2 Aia Write");
    }

    psio_->close(PSIF_DFMP2_AIA, 1);
    psio_->close(PSIF_DFMP2_QIA, 1);
}
void UDFMP2::form_Bia() {
    SharedMatrix Jm12 = form_inverse_metric();
    apply_fitting(Jm12, PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_a_->colspi()[0] * (size_t)Cavir_a_->colspi()[0]);
    apply_fitting(Jm12, PSIF_DFMP2_QIA, ribasis_->nbf(), Caocc_b_->colspi()[0] * (size_t)Cavir_b_->colspi()[0]);
}
void UDFMP2::form_Bia_Cia() {
    SharedMatrix Jm12 = form_inverse_metric();
    apply_fitting(Jm12, PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_a_->colspi()[0] * (size_t)Cavir_a_->colspi()[0]);
    apply_fitting(Jm12, PSIF_DFMP2_QIA, ribasis_->nbf(), Caocc_b_->colspi()[0] * (size_t)Cavir_b_->colspi()[0]);
    apply_fitting_grad(Jm12, PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_a_->colspi()[0] * (size_t)Cavir_a_->colspi()[0]);
    apply_fitting_grad(Jm12, PSIF_DFMP2_QIA, ribasis_->nbf(), Caocc_b_->colspi()[0] * (size_t)Cavir_b_->colspi()[0]);
}
void UDFMP2::form_Bia_transpose() { throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented"); }
void UDFMP2::form_energy() {
    // Energy registers
    double e_ss = 0.0;
    double e_os = 0.0;

    /* => AA Terms <= */ {
        // Sizing
        int naux = ribasis_->nbf();
        int naocc = Caocc_a_->colspi()[0];
        int navir = Cavir_a_->colspi()[0];

        // Thread considerations
        int nthread = 1;
#ifdef _OPENMP
        nthread = Process::environment.get_n_threads();
#endif

        // Memory
        size_t Iab_memory = navir * (size_t)navir;
        size_t Qa_memory = naux * (size_t)navir;
        size_t doubles = ((size_t)(options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
        if (doubles < nthread * Iab_memory) {
            throw PSIEXCEPTION("DFMP2: Insufficient memory for Iab buffers. Reduce OMP Threads or increase memory.");
        }
        size_t remainder = doubles - nthread * Iab_memory;
        size_t max_i = remainder / (2L * Qa_memory);
        max_i = (max_i > naocc ? naocc : max_i);
        max_i = (max_i < 1L ? 1L : max_i);

        // Blocks
        std::vector<size_t> i_starts;
        i_starts.push_back(0L);
        for (size_t i = 0; i < naocc; i += max_i) {
            if (i + max_i >= naocc) {
                i_starts.push_back(naocc);
            } else {
                i_starts.push_back(i + max_i);
            }
        }
        // block_status(i_starts, __FILE__,__LINE__);

        // Tensor blocks
        auto Bia = std::make_shared<Matrix>("B(ia|Q)", max_i * (size_t)navir, naux);
        auto Bjb = std::make_shared<Matrix>("B(jb|Q)", max_i * (size_t)navir, naux);
        double** Biap = Bia->pointer();
        double** Bjbp = Bjb->pointer();

        std::vector<SharedMatrix> Iab;
        for (int i = 0; i < nthread; i++) {
            Iab.push_back(std::make_shared<Matrix>("Iab", navir, navir));
        }

        double* eps_aoccp = eps_aocc_a_->pointer();
        double* eps_avirp = eps_avir_a_->pointer();

        // Loop through pairs of blocks
        psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
        psio_address next_BIA = PSIO_ZERO;
        for (int block_i = 0; block_i < i_starts.size() - 1; block_i++) {
            // Sizing
            size_t istart = i_starts[block_i];
            size_t istop = i_starts[block_i + 1];
            size_t ni = istop - istart;

            // Read iaQ chunk
            timer_on("DFMP2 Bia Read");
            next_BIA = psio_get_address(PSIO_ZERO, sizeof(double) * (istart * navir * naux));
            psio_->read(PSIF_DFMP2_AIA, "B(ia|Q)", (char*)Biap[0], sizeof(double) * (ni * navir * naux), next_BIA,
                        &next_BIA);
            timer_off("DFMP2 Bia Read");

            for (int block_j = 0; block_j <= block_i; block_j++) {
                // Sizing
                size_t jstart = i_starts[block_j];
                size_t jstop = i_starts[block_j + 1];
                size_t nj = jstop - jstart;

                // Read iaQ chunk (if unique)
                timer_on("DFMP2 Bia Read");
                if (block_i == block_j) {
                    ::memcpy((void*)Bjbp[0], (void*)Biap[0], sizeof(double) * (ni * navir * naux));
                } else {
                    next_BIA = psio_get_address(PSIO_ZERO, sizeof(double) * (jstart * navir * naux));
                    psio_->read(PSIF_DFMP2_AIA, "B(ia|Q)", (char*)Bjbp[0], sizeof(double) * (nj * navir * naux),
                                next_BIA, &next_BIA);
                }
                timer_off("DFMP2 Bia Read");

#pragma omp parallel for schedule(dynamic) num_threads(nthread) reduction(+ : e_ss)
                for (long int ij = 0L; ij < ni * nj; ij++) {
                    // Sizing
                    size_t i = ij / nj + istart;
                    size_t j = ij % nj + jstart;
                    if (j > i) continue;

                    double perm_factor = (i == j ? 1.0 : 2.0);

                    // Which thread is this?
                    int thread = 0;
#ifdef _OPENMP
                    thread = omp_get_thread_num();
#endif
                    double** Iabp = Iab[thread]->pointer();

                    // Form the integral block (ia|jb) = (ia|Q)(Q|jb)
                    C_DGEMM('N', 'T', navir, navir, naux, 1.0, Biap[(i - istart) * navir], naux,
                            Bjbp[(j - jstart) * navir], naux, 0.0, Iabp[0], navir);

                    // Add the MP2 energy contributions
                    for (int a = 0; a < navir; a++) {
                        for (int b = 0; b < navir; b++) {
                            double iajb = Iabp[a][b];
                            double ibja = Iabp[b][a];
                            double denom = -perm_factor / (eps_avirp[a] + eps_avirp[b] - eps_aoccp[i] - eps_aoccp[j]);

                            e_ss += 0.5 * (iajb * iajb - iajb * ibja) * denom;
                        }
                    }
                }
            }
        }
        psio_->close(PSIF_DFMP2_AIA, 1);

    /* End AA Terms */ }

    /* => BB Terms <= */ {
        // Sizing
        int naux = ribasis_->nbf();
        int naocc = Caocc_b_->colspi()[0];
        int navir = Cavir_b_->colspi()[0];

        // Thread considerations
        int nthread = 1;
#ifdef _OPENMP
        nthread = Process::environment.get_n_threads();
#endif

        // Memory
        size_t Iab_memory = navir * (size_t)navir;
        size_t Qa_memory = naux * (size_t)navir;
        size_t doubles = ((size_t)(options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
        if (doubles < nthread * Iab_memory) {
            throw PSIEXCEPTION("DFMP2: Insufficient memory for Iab buffers. Reduce OMP Threads or increase memory.");
        }
        size_t remainder = doubles - nthread * Iab_memory;
        size_t max_i = remainder / (2L * Qa_memory);
        max_i = (max_i > naocc ? naocc : max_i);
        max_i = (max_i < 1L ? 1L : max_i);

        // Blocks
        std::vector<size_t> i_starts;
        i_starts.push_back(0L);
        for (size_t i = 0; i < naocc; i += max_i) {
            if (i + max_i >= naocc) {
                i_starts.push_back(naocc);
            } else {
                i_starts.push_back(i + max_i);
            }
        }
        // block_status(i_starts, __FILE__,__LINE__);

        // Tensor blocks
        auto Bia = std::make_shared<Matrix>("B(ia|Q)", max_i * (size_t)navir, naux);
        auto Bjb = std::make_shared<Matrix>("Q(jb|Q)", max_i * (size_t)navir, naux);
        double** Biap = Bia->pointer();
        double** Bjbp = Bjb->pointer();

        std::vector<SharedMatrix> Iab;
        for (int i = 0; i < nthread; i++) {
            Iab.push_back(std::make_shared<Matrix>("Iab", navir, navir));
        }

        double* eps_aoccp = eps_aocc_b_->pointer();
        double* eps_avirp = eps_avir_b_->pointer();

        // Loop through pairs of blocks
        psio_->open(PSIF_DFMP2_QIA, PSIO_OPEN_OLD);
        psio_address next_BIA = PSIO_ZERO;
        for (int block_i = 0; block_i < i_starts.size() - 1; block_i++) {
            // Sizing
            size_t istart = i_starts[block_i];
            size_t istop = i_starts[block_i + 1];
            size_t ni = istop - istart;

            // Read iaQ chunk
            timer_on("DFMP2 Bia Read");
            next_BIA = psio_get_address(PSIO_ZERO, sizeof(double) * (istart * navir * naux));
            psio_->read(PSIF_DFMP2_QIA, "B(ia|Q)", (char*)Biap[0], sizeof(double) * (ni * navir * naux), next_BIA,
                        &next_BIA);
            timer_off("DFMP2 Bia Read");

            for (int block_j = 0; block_j <= block_i; block_j++) {
                // Sizing
                size_t jstart = i_starts[block_j];
                size_t jstop = i_starts[block_j + 1];
                size_t nj = jstop - jstart;

                // Read iaQ chunk (if unique)
                timer_on("DFMP2 Bia Read");
                if (block_i == block_j) {
                    ::memcpy((void*)Bjbp[0], (void*)Biap[0], sizeof(double) * (ni * navir * naux));
                } else {
                    next_BIA = psio_get_address(PSIO_ZERO, sizeof(double) * (jstart * navir * naux));
                    psio_->read(PSIF_DFMP2_QIA, "B(ia|Q)", (char*)Bjbp[0], sizeof(double) * (nj * navir * naux),
                                next_BIA, &next_BIA);
                }
                timer_off("DFMP2 Bia Read");

#pragma omp parallel for schedule(dynamic) num_threads(nthread) reduction(+ : e_ss)
                for (long int ij = 0L; ij < ni * nj; ij++) {
                    // Sizing
                    size_t i = ij / nj + istart;
                    size_t j = ij % nj + jstart;
                    if (j > i) continue;

                    double perm_factor = (i == j ? 1.0 : 2.0);

                    // Which thread is this?
                    int thread = 0;
#ifdef _OPENMP
                    thread = omp_get_thread_num();
#endif
                    double** Iabp = Iab[thread]->pointer();

                    // Form the integral block (ia|jb) = (ia|Q)(Q|jb)
                    C_DGEMM('N', 'T', navir, navir, naux, 1.0, Biap[(i - istart) * navir], naux,
                            Bjbp[(j - jstart) * navir], naux, 0.0, Iabp[0], navir);

                    // Add the MP2 energy contributions
                    for (int a = 0; a < navir; a++) {
                        for (int b = 0; b < navir; b++) {
                            double iajb = Iabp[a][b];
                            double ibja = Iabp[b][a];
                            double denom = -perm_factor / (eps_avirp[a] + eps_avirp[b] - eps_aoccp[i] - eps_aoccp[j]);

                            e_ss += 0.5 * (iajb * iajb - iajb * ibja) * denom;
                        }
                    }
                }
            }
        }
        psio_->close(PSIF_DFMP2_QIA, 1);

    /* End BB Terms */ }

    /* => AB Terms <= */ {
        // Sizing
        int naux = ribasis_->nbf();
        int naocc_a = Caocc_a_->colspi()[0];
        int navir_a = Cavir_a_->colspi()[0];
        int naocc_b = Caocc_b_->colspi()[0];
        int navir_b = Cavir_b_->colspi()[0];
        int naocc = (naocc_a > naocc_b ? naocc_a : naocc_b);
        int navir = (navir_a > navir_b ? navir_a : navir_b);

        // Thread considerations
        int nthread = 1;
#ifdef _OPENMP
        nthread = Process::environment.get_n_threads();
#endif

        // Memory
        size_t Iab_memory = navir_a * (size_t)navir_b;
        size_t Qa_memory = naux * (size_t)navir_a;
        size_t Qb_memory = naux * (size_t)navir_b;
        size_t doubles = ((size_t)(options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
        if (doubles < nthread * Iab_memory) {
            throw PSIEXCEPTION("DFMP2: Insufficient memory for Iab buffers. Reduce OMP Threads or increase memory.");
        }
        size_t remainder = doubles - nthread * Iab_memory;
        size_t max_i = remainder / (Qa_memory + Qb_memory);
        max_i = (max_i > naocc ? naocc : max_i);
        max_i = (max_i < 1L ? 1L : max_i);

        // Blocks
        std::vector<size_t> i_starts_a;
        i_starts_a.push_back(0L);
        for (size_t i = 0; i < naocc_a; i += max_i) {
            if (i + max_i >= naocc_a) {
                i_starts_a.push_back(naocc_a);
            } else {
                i_starts_a.push_back(i + max_i);
            }
        }
        std::vector<size_t> i_starts_b;
        i_starts_b.push_back(0L);
        for (size_t i = 0; i < naocc_b; i += max_i) {
            if (i + max_i >= naocc_b) {
                i_starts_b.push_back(naocc_b);
            } else {
                i_starts_b.push_back(i + max_i);
            }
        }

        // Tensor blocks
        auto Bia = std::make_shared<Matrix>("B(ia|Q)", max_i * (size_t)navir_a, naux);
        auto Bjb = std::make_shared<Matrix>("B(jb|Q)", max_i * (size_t)navir_b, naux);
        double** Biap = Bia->pointer();
        double** Bjbp = Bjb->pointer();

        std::vector<SharedMatrix> Iab;
        for (int i = 0; i < nthread; i++) {
            Iab.push_back(std::make_shared<Matrix>("Iab", navir_a, navir_b));
        }

        double* eps_aoccap = eps_aocc_a_->pointer();
        double* eps_avirap = eps_avir_a_->pointer();
        double* eps_aoccbp = eps_aocc_b_->pointer();
        double* eps_avirbp = eps_avir_b_->pointer();

        // Loop through pairs of blocks
        psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
        psio_->open(PSIF_DFMP2_QIA, PSIO_OPEN_OLD);
        psio_address next_BIAa = PSIO_ZERO;
        psio_address next_BIAb = PSIO_ZERO;
        for (int block_i = 0; block_i < i_starts_a.size() - 1; block_i++) {
            // Sizing
            size_t istart = i_starts_a[block_i];
            size_t istop = i_starts_a[block_i + 1];
            size_t ni = istop - istart;

            // Read iaQ chunk
            timer_on("DFMP2 Bia Read");
            next_BIAa = psio_get_address(PSIO_ZERO, sizeof(double) * (istart * navir_a * naux));
            psio_->read(PSIF_DFMP2_AIA, "B(ia|Q)", (char*)Biap[0], sizeof(double) * (ni * navir_a * naux), next_BIAa,
                        &next_BIAa);
            timer_off("DFMP2 Bia Read");

            for (int block_j = 0; block_j < i_starts_b.size() - 1; block_j++) {
                // Sizing
                size_t jstart = i_starts_b[block_j];
                size_t jstop = i_starts_b[block_j + 1];
                size_t nj = jstop - jstart;

                // Read iaQ chunk
                timer_on("DFMP2 Bia Read");
                next_BIAb = psio_get_address(PSIO_ZERO, sizeof(double) * (jstart * navir_b * naux));
                psio_->read(PSIF_DFMP2_QIA, "B(ia|Q)", (char*)Bjbp[0], sizeof(double) * (nj * navir_b * naux), next_BIAb,
                            &next_BIAb);
                timer_off("DFMP2 Bia Read");

#pragma omp parallel for schedule(dynamic) num_threads(nthread) reduction(+ : e_os)
                for (long int ij = 0L; ij < ni * nj; ij++) {
                    // Sizing
                    size_t i = ij / nj + istart;
                    size_t j = ij % nj + jstart;

                    // Which thread is this?
                    int thread = 0;
#ifdef _OPENMP
                    thread = omp_get_thread_num();
#endif
                    double** Iabp = Iab[thread]->pointer();

                    // Form the integral block (ia|jb) = (ia|Q)(Q|jb)
                    C_DGEMM('N', 'T', navir_a, navir_b, naux, 1.0, Biap[(i - istart) * navir_a], naux,
                            Bjbp[(j - jstart) * navir_b], naux, 0.0, Iabp[0], navir_b);

                    // Add the MP2 energy contributions
                    for (int a = 0; a < navir_a; a++) {
                        for (int b = 0; b < navir_b; b++) {
                            double iajb = Iabp[a][b];
                            double denom = -1.0 / (eps_avirap[a] + eps_avirbp[b] - eps_aoccap[i] - eps_aoccbp[j]);
                            e_os += (iajb * iajb) * denom;
                        }
                    }
                }
            }
        }
        psio_->close(PSIF_DFMP2_AIA, 0);
        psio_->close(PSIF_DFMP2_QIA, 0);

    /* End AB Terms */ }

    variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = e_ss;
    variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = e_os;
}
void UDFMP2::form_Pab() { throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented"); }
void UDFMP2::form_Pij() { throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented"); }
void UDFMP2::form_gamma() { throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented"); }
void UDFMP2::form_G_transpose() { throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented"); }
void UDFMP2::form_AB_x_terms() { throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented"); }
void UDFMP2::form_Amn_x_terms() { throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented"); }
void UDFMP2::form_L() { throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented"); }
void UDFMP2::form_P() { throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented"); }
void UDFMP2::form_W() { throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented"); }
void UDFMP2::form_Z() { throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented"); }
void UDFMP2::form_gradient() { throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented"); }

RODFMP2::RODFMP2(SharedWavefunction ref_wfn, Options& options, std::shared_ptr<PSIO> psio)
    : UDFMP2(ref_wfn, options, psio) {
    common_init();
}
RODFMP2::~RODFMP2() {}
void RODFMP2::common_init() {}
void RODFMP2::print_header() {
    int nthread = 1;
#ifdef _OPENMP
    nthread = Process::environment.get_n_threads();
#endif

    outfile->Printf("\t --------------------------------------------------------\n");
    outfile->Printf("\t                          DF-MP2                         \n");
    outfile->Printf("\t      2nd-Order Density-Fitted Moller-Plesset Theory     \n");
    outfile->Printf("\t          ROHF-MBPT(2) Wavefunction, %3d Threads         \n", nthread);
    outfile->Printf("\t                                                         \n");
    outfile->Printf("\t        Rob Parrish, Justin Turney, Andy Simmonett,      \n");
    outfile->Printf("\t           Ed Hohenstein, and C. David Sherrill          \n");
    outfile->Printf("\t --------------------------------------------------------\n");
    outfile->Printf("\n");

    int focc_a = frzcpi_.sum();
    int fvir_a = frzvpi_.sum();
    int aocc_a = Caocc_a_->colspi()[0];
    int avir_a = Cavir_a_->colspi()[0];
    int occ_a = focc_a + aocc_a;
    int vir_a = fvir_a + avir_a;

    int focc_b = frzcpi_.sum();
    int fvir_b = frzvpi_.sum();
    int aocc_b = Caocc_b_->colspi()[0];
    int avir_b = Cavir_b_->colspi()[0];
    int occ_b = focc_b + aocc_b;
    int vir_b = fvir_b + avir_b;

    outfile->Printf("\t --------------------------------------------------------\n");
    outfile->Printf("\t                 NBF = %5d, NAUX = %5d\n", basisset_->nbf(), ribasis_->nbf());
    outfile->Printf("\t --------------------------------------------------------\n");
    outfile->Printf("\t %7s %7s %7s %7s %7s %7s %7s\n", "CLASS", "FOCC", "OCC", "AOCC", "AVIR", "VIR", "FVIR");
    outfile->Printf("\t %7s %7d %7d %7d %7d %7d %7d\n", "ALPHA", focc_a, occ_a, aocc_a, avir_a, vir_a, fvir_a);
    outfile->Printf("\t %7s %7d %7d %7d %7d %7d %7d\n", "BETA", focc_b, occ_b, aocc_b, avir_b, vir_b, fvir_b);
    outfile->Printf("\t --------------------------------------------------------\n\n");
}
}  // namespace dfmp2
}  // namespace psi
