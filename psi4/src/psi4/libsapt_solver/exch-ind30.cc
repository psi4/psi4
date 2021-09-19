/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

#include "sapt2p3.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/lib3index/dfhelper.h"
#include "psi4/libfock/jk.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/process.h"

namespace psi {
namespace sapt {

void SAPT2p3::exch_ind30() {
    double **tAR = block_matrix(noccA_, nvirA_);
    double **vAR = block_matrix(noccA_, nvirA_);

    psio_->read_entry(PSIF_SAPT_AMPS, "Ind30 uAR Amplitudes", (char *)tAR[0], sizeof(double) * noccA_ * nvirA_);
    psio_->read_entry(PSIF_SAPT_AMPS, "AR Exch-Ind Integrals", (char *)vAR[0], sizeof(double) * noccA_ * nvirA_);

    double ex_1 = -2.0 * C_DDOT(noccA_ * nvirA_, tAR[0], 1, vAR[0], 1);

    free_block(tAR);
    free_block(vAR);

    double **tBS = block_matrix(noccB_, nvirB_);
    double **vBS = block_matrix(noccB_, nvirB_);

    psio_->read_entry(PSIF_SAPT_AMPS, "Ind30 uBS Amplitudes", (char *)tBS[0], sizeof(double) * noccB_ * nvirB_);
    psio_->read_entry(PSIF_SAPT_AMPS, "BS Exch-Ind Integrals", (char *)vBS[0], sizeof(double) * noccB_ * nvirB_);

    double ex_2 = -2.0 * C_DDOT(noccB_ * nvirB_, tBS[0], 1, vBS[0], 1);

    free_block(tBS);
    free_block(vBS);

    double **sAR = block_matrix(noccA_, nvirA_);
    double **sBS = block_matrix(noccB_, nvirB_);

    for (int a = 0; a < noccA_; a++) {
        for (int r = 0; r < nvirA_; r++) {
            sAR[a][r] = wBAR_[a][r] / (evalsA_[a] - evalsA_[r + noccA_]);
        }
    }

    for (int b = 0; b < noccB_; b++) {
        for (int s = 0; s < nvirB_; s++) {
            sBS[b][s] = wABS_[b][s] / (evalsB_[b] - evalsB_[s + noccB_]);
        }
    }

    double ex_3 = exch_ind30_1(sAR, sBS);
    double ex_4 = exch_ind30_2(sAR);
    double ex_5 = exch_ind30_3(sBS);

    free_block(sAR);
    free_block(sBS);

    e_exch_ind30_ = ex_1 + ex_2 + ex_3 + ex_4 + ex_5;

    if (debug_) {
        outfile->Printf("\n    Exch-Ind_1          = %18.12lf [Eh]\n", ex_1);
        outfile->Printf("    Exch-Ind_2          = %18.12lf [Eh]\n", ex_2);
        outfile->Printf("    Exch-Ind_3          = %18.12lf [Eh]\n", ex_3);
        outfile->Printf("    Exch-Ind_4          = %18.12lf [Eh]\n", ex_4);
        outfile->Printf("    Exch-Ind_5          = %18.12lf [Eh]\n", ex_5);
    }
    if (print_) {
        outfile->Printf("    Exch-Ind30          = %18.12lf [Eh]\n", e_exch_ind30_);
    }
}

double SAPT2p3::exch_ind30_1(double **sAR, double **sBS) {
    double energy = 0.0;

    double **vARBS = block_matrix(noccA_ * nvirA_, noccB_ * nvirB_);
    psio_->read_entry(PSIF_SAPT_AMPS, "Exch-Disp V_ARBS", (char *)vARBS[0],
                      sizeof(double) * noccA_ * nvirA_ * noccB_ * nvirB_);

    for (int a = 0, ar = 0; a < noccA_; a++) {
        for (int r = 0; r < nvirA_; r++, ar++) {
            energy -= 2.0 * sAR[a][r] * C_DDOT(noccB_ * nvirB_, vARBS[ar], 1, sBS[0], 1);
        }
    }

    free_block(vARBS);

    double **xAR = block_matrix(noccA_, nvirA_);
    double **xBS = block_matrix(noccB_, nvirB_);

    C_DGEMM('N', 'T', noccA_, nvirA_, noccB_, 1.0, sAB_[0], nmoB_, sAB_[noccA_], nmoB_, 0.0, xAR[0], nvirA_);

    C_DGEMM('T', 'N', noccB_, nvirB_, noccA_, 1.0, sAB_[0], nmoB_, &(sAB_[0][noccB_]), nmoB_, 0.0, xBS[0], nvirB_);

    double **yAR = block_matrix(noccA_, nvirA_);
    double **yBS = block_matrix(noccB_, nvirB_);

    double **B_p_AR = get_AR_ints(1);
    double **B_p_BS = get_BS_ints(1);

    C_DGEMV('n', noccA_ * nvirA_, ndf_ + 3, 1.0, &(B_p_AR[0][0]), ndf_ + 3, diagBB_, 1, 0.0, &(yAR[0][0]), 1);

    C_DGEMV('n', noccB_ * nvirB_, ndf_ + 3, 1.0, &(B_p_BS[0][0]), ndf_ + 3, diagAA_, 1, 0.0, &(yBS[0][0]), 1);

    energy += 8.0 * C_DDOT(noccA_ * nvirA_, xAR[0], 1, sAR[0], 1) * C_DDOT(noccB_ * nvirB_, yBS[0], 1, sBS[0], 1);
    energy += 8.0 * C_DDOT(noccA_ * nvirA_, yAR[0], 1, sAR[0], 1) * C_DDOT(noccB_ * nvirB_, xBS[0], 1, sBS[0], 1);

    free_block(B_p_AR);
    free_block(B_p_BS);
    free_block(xAR);
    free_block(xBS);
    free_block(yAR);
    free_block(yBS);

    return (energy);
}

double SAPT2p3::exch_ind30_2(double **sAR) {
    double energy = 0.0;

    double **ssRB = block_matrix(nvirA_, noccB_);

    C_DGEMM('T', 'N', nvirA_, noccB_, noccA_, 1.0, sAR[0], nvirA_, sAB_[0], nmoB_, 0.0, ssRB[0], noccB_);

    double **A_p_AR = get_AR_ints(1);
    double **B_p_BB = get_BB_ints(1);
    double **B_p_RB = get_RB_ints(1);

    double **C_p_AB = block_matrix(noccA_ * noccB_, ndf_ + 3);
    double **D_p_AB = block_matrix(noccA_ * noccB_, ndf_ + 3);

    for (int a = 0; a < noccA_; a++) {
        C_DGEMM('T', 'N', noccB_, ndf_ + 3, nvirA_, 1.0, ssRB[0], noccB_, A_p_AR[a * nvirA_], ndf_ + 3, 0.0,
                C_p_AB[a * noccB_], ndf_ + 3);
    }

    C_DGEMM('N', 'N', noccA_, noccB_ * (ndf_ + 3), nvirA_, 1.0, sAR[0], nvirA_, B_p_RB[0], noccB_ * (ndf_ + 3), 0.0,
            D_p_AB[0], noccB_ * (ndf_ + 3));

    energy += 2.0 * C_DDOT(noccA_ * noccB_ * (ndf_ + 3), C_p_AB[0], 1, D_p_AB[0], 1);

    free_block(C_p_AB);
    free_block(D_p_AB);

    double *X = init_array(ndf_ + 3);
    double *Y = init_array(ndf_ + 3);

    C_DGEMV('t', noccA_ * nvirA_, ndf_ + 3, 1.0, A_p_AR[0], ndf_ + 3, sAR[0], 1, 0.0, X, 1);
    C_DGEMV('t', nvirA_ * noccB_, ndf_ + 3, 1.0, B_p_RB[0], ndf_ + 3, ssRB[0], 1, 0.0, Y, 1);

    energy -= 4.0 * C_DDOT(ndf_ + 3, X, 1, Y, 1);

    double **xAB = block_matrix(noccA_, noccB_);
    double **xAR = block_matrix(noccA_, nvirA_);
    double **yAR = block_matrix(noccA_, nvirA_);

    C_DGEMM('N', 'N', noccA_, noccB_, nvirA_, 1.0, sAR[0], nvirA_, sAB_[noccA_], nmoB_, 0.0, xAB[0], noccB_);

    C_DGEMM('N', 'T', noccA_, nvirA_, noccB_, 1.0, xAB[0], noccB_, ssRB[0], noccB_, 0.0, xAR[0], nvirA_);

    C_DGEMV('n', noccA_ * nvirA_, ndf_ + 3, 1.0, A_p_AR[0], ndf_ + 3, diagBB_, 1, 0.0, yAR[0], 1);

    energy += 4.0 * C_DDOT(noccA_ * nvirA_, xAR[0], 1, yAR[0], 1);

    free_block(xAR);
    free_block(yAR);

    double **E_p_AB = block_matrix(noccA_ * noccB_, ndf_ + 3);
    double **C_p_BB = block_matrix(noccB_ * noccB_, ndf_ + 3);

    for (int a = 0; a < noccA_; a++) {
        C_DGEMM('T', 'N', noccB_, ndf_ + 3, nvirA_, 1.0, ssRB[0], noccB_, A_p_AR[a * nvirA_], ndf_ + 3, 0.0,
                E_p_AB[a * noccB_], ndf_ + 3);
    }

    C_DGEMM('T', 'N', noccB_, noccB_ * (ndf_ + 3), noccA_, 1.0, xAB[0], noccB_, E_p_AB[0], noccB_ * (ndf_ + 3), 0.0,
            C_p_BB[0], noccB_ * (ndf_ + 3));

    energy -= 2.0 * C_DDOT(noccB_ * noccB_ * (ndf_ + 3), C_p_BB[0], 1, B_p_BB[0], 1);

    free_block(xAB);
    free_block(E_p_AB);
    free_block(C_p_BB);

    double **xBB = block_matrix(noccB_, noccB_);

    C_DGEMM('T', 'N', noccB_, noccB_, nvirA_, 1.0, ssRB[0], noccB_, sAB_[noccA_], nmoB_, 0.0, xBB[0], noccB_);

    C_DGEMV('t', noccB_ * noccB_, ndf_ + 3, 1.0, B_p_BB[0], ndf_ + 3, xBB[0], 1, 0.0, Y, 1);

    energy += 4.0 * C_DDOT(ndf_ + 3, X, 1, Y, 1);

    free_block(xBB);
    free_block(ssRB);
    free(X);
    free(Y);

    free_block(A_p_AR);
    free_block(B_p_RB);
    free_block(B_p_BB);

    return (energy);
}

double SAPT2p3::exch_ind30_3(double **sBS) {
    double energy = 0.0;

    double **ssAS = block_matrix(noccA_, nvirB_);

    C_DGEMM('N', 'N', noccA_, nvirB_, noccB_, 1.0, sAB_[0], nmoB_, sBS[0], nvirB_, 0.0, ssAS[0], nvirB_);

    double **A_p_AA = get_AA_ints(1);
    double **A_p_AS = get_AS_ints(1);
    double **B_p_BS = get_BS_ints(1);

    double **C_p_AB = block_matrix(noccA_ * noccB_, ndf_ + 3);
    double **D_p_AB = block_matrix(noccA_ * noccB_, ndf_ + 3);

    for (int b = 0; b < noccB_; b++) {
        C_DGEMM('N', 'N', noccA_, ndf_ + 3, nvirB_, 1.0, ssAS[0], nvirB_, B_p_BS[b * nvirB_], ndf_ + 3, 0.0, C_p_AB[b],
                noccB_ * (ndf_ + 3));
    }

    for (int a = 0; a < noccA_; a++) {
        C_DGEMM('N', 'N', noccB_, ndf_ + 3, nvirB_, 1.0, sBS[0], nvirB_, A_p_AS[a * nvirB_], ndf_ + 3, 0.0,
                D_p_AB[a * noccB_], ndf_ + 3);
    }

    energy += 2.0 * C_DDOT(noccA_ * noccB_ * (ndf_ + 3), C_p_AB[0], 1, D_p_AB[0], 1);

    free_block(C_p_AB);
    free_block(D_p_AB);

    double *X = init_array(ndf_ + 3);
    double *Y = init_array(ndf_ + 3);

    C_DGEMV('t', noccB_ * nvirB_, ndf_ + 3, 1.0, B_p_BS[0], ndf_ + 3, sBS[0], 1, 0.0, X, 1);
    C_DGEMV('t', noccA_ * nvirB_, ndf_ + 3, 1.0, A_p_AS[0], ndf_ + 3, ssAS[0], 1, 0.0, Y, 1);

    energy -= 4.0 * C_DDOT(ndf_ + 3, X, 1, Y, 1);

    double **xAB = block_matrix(noccA_, noccB_);
    double **xBS = block_matrix(noccB_, nvirB_);
    double **yBS = block_matrix(noccB_, nvirB_);

    C_DGEMM('N', 'T', noccA_, noccB_, nvirB_, 1.0, &(sAB_[0][noccB_]), nmoB_, sBS[0], nvirB_, 0.0, xAB[0], noccB_);

    C_DGEMM('T', 'N', noccB_, nvirB_, noccA_, 1.0, xAB[0], noccB_, ssAS[0], nvirB_, 0.0, xBS[0], nvirB_);

    C_DGEMV('n', noccB_ * nvirB_, ndf_ + 3, 1.0, B_p_BS[0], ndf_ + 3, diagAA_, 1, 0.0, yBS[0], 1);

    energy += 4.0 * C_DDOT(noccB_ * nvirB_, xBS[0], 1, yBS[0], 1);

    free_block(xBS);
    free_block(yBS);

    double **E_p_AB = block_matrix(noccA_ * noccB_, ndf_ + 3);
    double **C_p_AA = block_matrix(noccA_ * noccA_, ndf_ + 3);

    for (int b = 0; b < noccB_; b++) {
        C_DGEMM('N', 'N', noccA_, ndf_ + 3, nvirB_, 1.0, ssAS[0], nvirB_, B_p_BS[b * nvirB_], ndf_ + 3, 0.0,
                E_p_AB[b * noccA_], ndf_ + 3);
    }

    C_DGEMM('N', 'N', noccA_, noccA_ * (ndf_ + 3), noccB_, 1.0, xAB[0], noccB_, E_p_AB[0], noccA_ * (ndf_ + 3), 0.0,
            C_p_AA[0], noccA_ * (ndf_ + 3));

    energy -= 2.0 * C_DDOT(noccA_ * noccA_ * (ndf_ + 3), C_p_AA[0], 1, A_p_AA[0], 1);

    free_block(xAB);
    free_block(E_p_AB);
    free_block(C_p_AA);

    double **xAA = block_matrix(noccA_, noccA_);

    C_DGEMM('N', 'T', noccA_, noccA_, nvirB_, 1.0, ssAS[0], nvirB_, &(sAB_[0][noccB_]), nmoB_, 0.0, xAA[0], noccA_);

    C_DGEMV('t', noccA_ * noccA_, ndf_ + 3, 1.0, A_p_AA[0], ndf_ + 3, xAA[0], 1, 0.0, Y, 1);

    energy += 4.0 * C_DDOT(ndf_ + 3, X, 1, Y, 1);

    free_block(xAA);
    free(X);
    free(Y);

    free_block(ssAS);
    free_block(A_p_AS);
    free_block(A_p_AA);
    free_block(B_p_BS);

    return (energy);
}

// Compute the nonapproximated third-order exchange-induction //
// Konrad Patkowski, based on Jonathan Waldrop's Psi4NumPy code //
// September 2021 //
void SAPT2p3::sinf_e30ind() {

    if (print_) {
        outfile->Printf("  ==> Nonapproximated third-order induction <==\n\n");
    }

    // => Sizing <= //

    int nn = basisset_->nbf();
    int na = noccA_;
    int nb = noccB_;
    int nr = nvirA_;
    int ns = nvirB_;
    int nQ = ribasis_->nbf();

    int nT = 1;
#ifdef _OPENMP
    nT = Process::environment.get_n_threads();
#endif

    double **uAR = block_matrix(na, nr);
    double **uBS = block_matrix(nb, ns);
    psio_->read_entry(PSIF_SAPT_AMPS, "Ind30 uAR Amplitudes", (char *)uAR[0], sizeof(double) * na * nr);
    psio_->read_entry(PSIF_SAPT_AMPS, "Ind30 uBS Amplitudes", (char *)uBS[0], sizeof(double) * nb * ns);

    outfile->Printf("%d %d\n",CoccA_->rowdim(),CoccA_->coldim());
    outfile->Printf("%d %d\n",CoccB_->rowdim(),CoccB_->coldim());
    outfile->Printf("%d %d\n",Smat_->rowdim(),Smat_->coldim());
   
    // => Intermolecular overlap matrix and inverse <= //
    std::shared_ptr<Matrix> Sab = linalg::triplet(CoccA_, Smat_, CoccB_, true, false, false);

    double** Sabp = Sab->pointer();
    auto D = std::make_shared<Matrix>("D", na + nb, na + nb);
    D->identity();
    double** Dp = D->pointer();
    for (int a = 0; a < na; a++) {
        for (int b = 0; b < nb; b++) {
            Dp[a][b + na] = Dp[b + na][a] = Sabp[a][b];
        }
    }
    D->power(-1.0, 1.0E-12);
    Dp = D->pointer();

    // => New Stuff <= //
    // Start with T's
    std::shared_ptr<Matrix> Sbr = linalg::triplet(CoccB_, Smat_, CvirA_, true, false, false);
    std::shared_ptr<Matrix> Sas = linalg::triplet(CoccA_, Smat_, CvirB_, true, false, false);
    auto Tar = std::make_shared<Matrix>("Tar", na, nr);
    auto Tbr = std::make_shared<Matrix>("Tbr", nb, nr);
    auto Tas = std::make_shared<Matrix>("Tas", na, ns);
    auto Tbs = std::make_shared<Matrix>("Tbs", nb, ns);

    C_DGEMM('N', 'N', na, nr, nb, 1.0, &Dp[0][na], na + nb, Sbr->pointer()[0], nr, 0.0,
            Tar->pointer()[0], nr);
    C_DGEMM('N', 'N', nb, nr, nb, 1.0, &Dp[na][na], na + nb, Sbr->pointer()[0], nr, 0.0,
            Tbr->pointer()[0], nr);
    C_DGEMM('N', 'N', na, ns, na, 1.0, &Dp[0][0], na + nb, Sas->pointer()[0], ns, 0.0,
            Tas->pointer()[0], ns);
    C_DGEMM('N', 'N', nb, ns, na, 1.0, &Dp[na][0], na + nb, Sas->pointer()[0], ns, 0.0,
            Tbs->pointer()[0], ns);

    // C1's and C2's from D's and C's.
    // C1's are C times D diagonal blocks.
    // C2's are times off-diagonal blocks.
    auto C1a = std::make_shared<Matrix>("C1a", nn, na);
    auto C1b = std::make_shared<Matrix>("C1b", nn, nb);
    auto C2a = std::make_shared<Matrix>("C2a", nn, na);
    auto C2b = std::make_shared<Matrix>("C2b", nn, nb);

    C_DGEMM('N', 'N', nn, na, na, 1.0, CoccA_->pointer()[0], na, &Dp[0][0], na + nb, 0.0,
            C1a->pointer()[0], na);
    C_DGEMM('N', 'N', nn, nb, nb, 1.0, CoccB_->pointer()[0], nb, &Dp[na][na], na + nb, 0.0,
            C1b->pointer()[0], nb);
    C_DGEMM('N', 'N', nn, na, nb, 1.0, CoccB_->pointer()[0], nb, &Dp[na][0], na + nb, 0.0,
            C2a->pointer()[0], na);
    C_DGEMM('N', 'N', nn, nb, na, 1.0, CoccA_->pointer()[0], na, &Dp[0][na], na + nb, 0.0,
            C2b->pointer()[0], nb);

    // Coeffs for all occupied
    std::vector<std::shared_ptr<Matrix> > hold_these;
    hold_these.push_back(CoccA_);
    hold_these.push_back(CoccB_);

    auto Cocc0AB = linalg::horzcat(hold_these);
    hold_these.clear();

    // Half transform D_ia and D_ib for JK
    auto D_Ni_a = std::make_shared<Matrix>("D_Ni_a", nn, na + nb);
    auto D_Ni_b = std::make_shared<Matrix>("D_Ni_b", nn, na + nb);

    C_DGEMM('N', 'N', nn, na + nb, na, 1.0, CoccA_->pointer()[0], na, &Dp[0][0], na + nb, 0.0,
            D_Ni_a->pointer()[0], na + nb);
    C_DGEMM('N', 'N', nn, na + nb, nb, 1.0, CoccB_->pointer()[0], nb, &Dp[na][0], na + nb, 0.0,
            D_Ni_b->pointer()[0], na + nb);

    // Global JK object
    std::shared_ptr<JK> jk_;

    // Make JK's
    jk_ = JK::build_JK(basisset_, get_basisset("DF_BASIS_SCF"), options_, false, mem_);
    jk_->set_memory(mem_);
    jk_->set_do_J(true);
    jk_->set_do_K(true);
    jk_->initialize();
    jk_->print_header();

    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    Cl.clear();
    Cr.clear();
    Cl.push_back(Cocc0AB);
    Cr.push_back(D_Ni_a);
    Cl.push_back(Cocc0AB);
    Cr.push_back(D_Ni_b);
    jk_->compute();

    std::shared_ptr<Matrix> J_D_ia = J[0];
    std::shared_ptr<Matrix> K_D_ia = K[0];

    std::shared_ptr<Matrix> J_D_ib = J[1];
    std::shared_ptr<Matrix> K_D_ib = K[1];

    // Finish D_ia and D_ib transformation to make tilded C's
    auto D_ia = linalg::doublet(Cocc0AB, D_Ni_a, false, true);
    auto D_ib = linalg::doublet(Cocc0AB, D_Ni_b, false, true);

    auto Ct_Kr = linalg::triplet(D_ib, Smat_, CvirA_, false, false, false);
    Ct_Kr->scale(-1);
    Ct_Kr->add(CvirA_);
    auto Ct_Ks = linalg::triplet(D_ia, Smat_, CvirB_, false, false, false);
    Ct_Ks->scale(-1);
    Ct_Ks->add(CvirB_);

    // Make omega cores
    std::shared_ptr<Matrix> AJK(J_D_ia->clone());
    AJK->zero();
    AJK->add(J_D_ia);
    AJK->scale(2);
    AJK->add(VAmat_);
    AJK->subtract(K_D_ia);

    auto AJK_ar = linalg::triplet(C2a, AJK, Ct_Kr, true, false, false);
    auto AJK_as = linalg::triplet(C2a, AJK, Ct_Ks, true, false, false);
    auto AJK_br = linalg::triplet(C1b, AJK, Ct_Kr, true, false, false);
    auto AJK_bs = linalg::triplet(C1b, AJK, Ct_Ks, true, false, false);

    std::shared_ptr<Matrix> BJK(J_D_ib->clone());
    BJK->zero();
    BJK->add(J_D_ib);
    BJK->scale(2);
    BJK->add(VBmat_);
    BJK->subtract(K_D_ib);

    auto BJK_ar = linalg::triplet(C1a, BJK, Ct_Kr, true, false, false);
    auto BJK_as = linalg::triplet(C1a, BJK, Ct_Ks, true, false, false);
    auto BJK_br = linalg::triplet(C2b, BJK, Ct_Kr, true, false, false);
    auto BJK_bs = linalg::triplet(C2b, BJK, Ct_Ks, true, false, false);

    // Finish omega terms
    std::shared_ptr<Matrix> omega_ar(AJK_ar->clone());
    omega_ar->zero();
    omega_ar->add(AJK_ar);
    omega_ar->add(BJK_ar);
    omega_ar->scale(2);

    std::shared_ptr<Matrix> omega_as(AJK_as->clone());
    omega_as->zero();
    omega_as->add(AJK_as);
    omega_as->add(BJK_as);
    omega_as->scale(2);

    std::shared_ptr<Matrix> omega_br(AJK_br->clone());
    omega_br->zero();
    omega_br->add(AJK_br);
    omega_br->add(BJK_br);
    omega_br->scale(2);

    std::shared_ptr<Matrix> omega_bs(AJK_bs->clone());
    omega_bs->zero();
    omega_bs->add(AJK_bs);
    omega_bs->add(BJK_bs);
    omega_bs->scale(2);

    Cocc0AB.reset();
    D_Ni_a.reset();
    D_Ni_b.reset();
    J_D_ia.reset();
    K_D_ia.reset();
    J_D_ib.reset();
    K_D_ib.reset();
    AJK.reset();
    BJK.reset();
    AJK_ar.reset();
    AJK_as.reset();
    AJK_br.reset();
    AJK_bs.reset();
    BJK_ar.reset();
    BJK_as.reset();
    BJK_br.reset();
    BJK_bs.reset();

    //Konrad - all the stuff above still needed but I adjusted the scale factors in \omega's
    //we don't need the DF stuff that was below
    
    double CompleteInd30 = 0.0;
  
    auto sAR = std::make_shared<Matrix>("sAR", na, nr); 
    double **sARp = sAR->pointer();

    for (int a = 0; a < na; a++) {
        for (int r = 0; r < nr; r++) {
            sARp[a][r] = wBAR_[a][r] / (evalsA_[a] - evalsA_[na+r]);
        }
    }

    auto sBS = std::make_shared<Matrix>("sBS", nb, ns); 
    double **sBSp = sBS->pointer();

    for (int b = 0; b < nb; b++) {
        for (int s = 0; s < ns; s++) {
            sBSp[b][s] = wABS_[b][s] / (evalsB_[b] - evalsB_[nb+s]);
        }
    }

    std::shared_ptr<Matrix> STS_br = linalg::triplet(sBS, Tas, sAR, false, true, false);
    std::shared_ptr<Matrix> STS_as = linalg::triplet(sAR, Tbr, sBS, false, true, false);
    std::shared_ptr<Matrix> STS_ar = linalg::triplet(sAR, Tar, sAR, false, true, false);
    std::shared_ptr<Matrix> STS_bs = linalg::triplet(sBS, Tbs, sBS, false, true, false);

    CompleteInd30 += C_DDOT(na * nr, uAR[0], 1, omega_ar->pointer()[0], 1);
    CompleteInd30 += C_DDOT(nb * ns, uBS[0], 1, omega_bs->pointer()[0], 1);
    CompleteInd30 -= C_DDOT(nb * nr, STS_br->pointer()[0], 1, omega_br->pointer()[0], 1);
    CompleteInd30 -= C_DDOT(na * ns, STS_as->pointer()[0], 1, omega_as->pointer()[0], 1);
    CompleteInd30 -= C_DDOT(na * nr, STS_ar->pointer()[0], 1, omega_ar->pointer()[0], 1);
    CompleteInd30 -= C_DDOT(nb * ns, STS_bs->pointer()[0], 1, omega_bs->pointer()[0], 1);

//All Omega-dependent contributions have been completed, now Xi-dependent (J/K) contributions

    auto SCt_Na = std::make_shared<Matrix>("SCt_Na", nn, na);
    auto SCt_Nb = std::make_shared<Matrix>("SCt_Nb", nn, nb);
    C_DGEMM('N', 'T', nn, na, nr, 1.0, Ct_Kr->pointer()[0], nr, sAR->pointer()[0], nr, 0.0,
            SCt_Na->pointer()[0], na);
    C_DGEMM('N', 'T', nn, nb, ns, 1.0, Ct_Ks->pointer()[0], ns, sBS->pointer()[0], ns, 0.0,
            SCt_Nb->pointer()[0], nb);

    auto preXiAA = std::make_shared<Matrix>("preXiAA", nn, na);
    auto preXiBA = std::make_shared<Matrix>("preXiBA", nn, nb);
    auto preXiAB = std::make_shared<Matrix>("preXiAB", nn, na);
    auto preXiBB = std::make_shared<Matrix>("preXiBB", nn, nb);
    C_DGEMM('N', 'N', nn, na, na, 1.0, SCt_Na->pointer()[0], na, &Dp[0][0], na + nb, 0.0,
            preXiAA->pointer()[0], na);
    C_DGEMM('N', 'N', nn, nb, na, 1.0, SCt_Na->pointer()[0], na, &Dp[0][na], na + nb, 0.0,
            preXiBA->pointer()[0], nb);
    C_DGEMM('N', 'N', nn, na, nb, 1.0, SCt_Nb->pointer()[0], nb, &Dp[na][0], na + nb, 0.0,
            preXiAB->pointer()[0], na);
    C_DGEMM('N', 'N', nn, nb, nb, 1.0, SCt_Nb->pointer()[0], nb, &Dp[na][na], na + nb, 0.0,
            preXiBB->pointer()[0], nb);

    Cl.clear();
    Cr.clear();
    Cl.push_back(CoccB_);
    Cr.push_back(preXiBA);
    Cl.push_back(CoccB_);
    Cr.push_back(preXiBB);
    jk_->compute();

    std::shared_ptr<Matrix> J_XiBA = J[0];
    std::shared_ptr<Matrix> K_XiBA = K[0]->transpose();
    std::shared_ptr<Matrix> J_XiBB = J[1];
    std::shared_ptr<Matrix> K_XiBB = K[1]->transpose();

    auto XiAA = std::make_shared<Matrix>("XiAA", nn, nn);
    auto XiAB = std::make_shared<Matrix>("XiAB", nn, nn);
    auto XiBA = std::make_shared<Matrix>("XiBA", nn, nn);
    auto XiBB = std::make_shared<Matrix>("XiBB", nn, nn);
    C_DGEMM('N', 'T', nn, nn, na, 1.0, CoccA_->pointer()[0], na, preXiAA->pointer()[0], na, 0.0,
            XiAA->pointer()[0], nn);
    C_DGEMM('N', 'T', nn, nn, na, 1.0, CoccA_->pointer()[0], na, preXiAB->pointer()[0], na, 0.0,
            XiAB->pointer()[0], nn);
    C_DGEMM('N', 'T', nn, nn, nb, 1.0, CoccB_->pointer()[0], nb, preXiBA->pointer()[0], nb, 0.0,
            XiBA->pointer()[0], nn);
    C_DGEMM('N', 'T', nn, nn, nb, 1.0, CoccB_->pointer()[0], nb, preXiBB->pointer()[0], nb, 0.0,
            XiBB->pointer()[0], nn);

    CompleteInd30 += 4.0 * C_DDOT(nn * nn, XiAA->pointer()[0], 1, J_XiBB->pointer()[0], 1);
    CompleteInd30 -= 2.0 * C_DDOT(nn * nn, XiAA->pointer()[0], 1, K_XiBB->pointer()[0], 1);
    CompleteInd30 += 4.0 * C_DDOT(nn * nn, XiAB->pointer()[0], 1, J_XiBA->pointer()[0], 1);
    CompleteInd30 -= 2.0 * C_DDOT(nn * nn, XiAB->pointer()[0], 1, K_XiBA->pointer()[0], 1);
    CompleteInd30 += 4.0 * C_DDOT(nn * nn, XiAA->pointer()[0], 1, J_XiBA->pointer()[0], 1);
    CompleteInd30 -= 2.0 * C_DDOT(nn * nn, XiAA->pointer()[0], 1, K_XiBA->pointer()[0], 1);
    CompleteInd30 += 4.0 * C_DDOT(nn * nn, XiAB->pointer()[0], 1, J_XiBB->pointer()[0], 1);
    CompleteInd30 -= 2.0 * C_DDOT(nn * nn, XiAB->pointer()[0], 1, K_XiBB->pointer()[0], 1);

    e_exch_ind30_sinf_ = CompleteInd30 - e_ind30_;

    Process::environment.globals["SAPT EXCH-IND30(S^INF) ENERGY"] = e_exch_ind30_sinf_;
    if (print_) {
        outfile->Printf("    Exch-Ind30  (S^inf) = %18.12lf [Eh]\n", e_exch_ind30_sinf_);
        outfile->Printf("\n");
    }
}

}  // namespace sapt
}  // namespace psi
