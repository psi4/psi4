/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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
#include "psi4/libfock/jk.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/process.h"

namespace psi {
namespace sapt {

void SAPT2p3::exch_ind30() {
    auto tAR = std::make_unique<Matrix>("Ind30 uAR Amplitudes", noccA_, nvirA_);
    double **vAR = block_matrix(noccA_, nvirA_);

    tAR->load(psio_, PSIF_SAPT_AMPS, Matrix::SaveType::SubBlocks);
    psio_->read_entry(PSIF_SAPT_AMPS, "AR Exch-Ind Integrals", (char *)vAR[0], sizeof(double) * noccA_ * nvirA_);

    double ex_1 = -2.0 * C_DDOT(noccA_ * nvirA_, tAR->get_pointer(), 1, vAR[0], 1);

    tAR.reset();
    free_block(vAR);

    auto tBS = std::make_unique<Matrix>("Ind30 uBS Amplitudes", noccB_, nvirB_);
    double **vBS = block_matrix(noccB_, nvirB_);

    tBS->load(psio_, PSIF_SAPT_AMPS, Matrix::SaveType::SubBlocks);
    psio_->read_entry(PSIF_SAPT_AMPS, "BS Exch-Ind Integrals", (char *)vBS[0], sizeof(double) * noccB_ * nvirB_);

    double ex_2 = -2.0 * C_DDOT(noccB_ * nvirB_, tBS->get_pointer(), 1, vBS[0], 1);

    tBS.reset();
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

    auto uAR = std::make_unique<Matrix>("Ind30 uAR Amplitudes", na, nr);
    auto uBS = std::make_unique<Matrix>("Ind30 uBS Amplitudes", nb, ns);
    uAR->load(psio_, PSIF_SAPT_AMPS, Matrix::SaveType::SubBlocks);
    uBS->load(psio_, PSIF_SAPT_AMPS, Matrix::SaveType::SubBlocks);

    // => Intermolecular overlap matrix and inverse <= //
    auto Sab = linalg::triplet(CoccA_, Smat_, CoccB_, true, false, false);

    double** Sabp = Sab->pointer();
    auto D = Matrix("D", na + nb, na + nb);
    D.identity();
    double** Dp = D.pointer();
    for (int a = 0; a < na; a++) {
        for (int b = 0; b < nb; b++) {
            Dp[a][b + na] = Dp[b + na][a] = Sabp[a][b];
        }
    }
    D.power(-1.0, 1.0E-12);

    Dimension zero(1);
    Dimension nadim(1);
    Dimension nabdim(1);
    zero[0] = 0;
    nadim[0] = na;
    nabdim[0] = na + nb;
    auto Daa = D.get_block(Slice(zero, nadim), Slice(zero, nadim));
    auto Dab = D.get_block(Slice(zero, nadim), Slice(nadim, nabdim));
    auto Dbb = D.get_block(Slice(nadim, nabdim), Slice(nadim, nabdim));

    // => New Stuff <= //
    // Start with T's
    auto Sbr = linalg::triplet(CoccB_, Smat_, CvirA_, true, false, false);
    auto Sas = linalg::triplet(CoccA_, Smat_, CvirB_, true, false, false);
    auto Tar = linalg::doublet(Dab, Sbr, false, false);
    auto Tbr = linalg::doublet(Dbb, Sbr, false, false);
    auto Tas = linalg::doublet(Daa, Sas, false, false);
    auto Tbs = linalg::doublet(Dab, Sas, true, false);

    // C1's and C2's from D's and C's.
    // C1's are C times D diagonal blocks.
    // C2's are times off-diagonal blocks.
    auto C1a = linalg::doublet(CoccA_, Daa, false, false);
    auto C1b = linalg::doublet(CoccB_, Dbb, false, false);
    auto C2a = linalg::doublet(CoccB_, Dab, false, true);
    auto C2b = linalg::doublet(CoccA_, Dab, false, false);

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

    auto& Cl = jk_->C_left();
    auto& Cr = jk_->C_right();
    const auto& J = jk_->J();
    const auto& K = jk_->K();

    Cl.clear();
    Cr.clear();
    Cl.push_back(Cocc0AB);
    Cr.push_back(D_Ni_a);
    Cl.push_back(Cocc0AB);
    Cr.push_back(D_Ni_b);
    jk_->compute();

    auto J_D_ia = J[0];
    auto K_D_ia = K[0];

    auto J_D_ib = J[1];
    auto K_D_ib = K[1];

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
    AJK->scale(2);
    AJK->add(VAmat_);
    AJK->subtract(K_D_ia);

    auto AJK_ar = linalg::triplet(C2a, AJK, Ct_Kr, true, false, false);
    auto AJK_as = linalg::triplet(C2a, AJK, Ct_Ks, true, false, false);
    auto AJK_br = linalg::triplet(C1b, AJK, Ct_Kr, true, false, false);
    auto AJK_bs = linalg::triplet(C1b, AJK, Ct_Ks, true, false, false);

    std::shared_ptr<Matrix> BJK(J_D_ib->clone());
    BJK->scale(2);
    BJK->add(VBmat_);
    BJK->subtract(K_D_ib);

    auto BJK_ar = linalg::triplet(C1a, BJK, Ct_Kr, true, false, false);
    auto BJK_as = linalg::triplet(C1a, BJK, Ct_Ks, true, false, false);
    auto BJK_br = linalg::triplet(C2b, BJK, Ct_Kr, true, false, false);
    auto BJK_bs = linalg::triplet(C2b, BJK, Ct_Ks, true, false, false);

    // Finish omega terms
    Matrix omega_ar(AJK_ar->clone());
    omega_ar.add(BJK_ar);
    omega_ar.scale(2);

    Matrix omega_as(AJK_as->clone());
    omega_as.add(BJK_as);
    omega_as.scale(2);

    Matrix omega_br(AJK_br->clone());
    omega_br.add(BJK_br);
    omega_br.scale(2);

    Matrix omega_bs(AJK_bs->clone());
    omega_bs.add(BJK_bs);
    omega_bs.scale(2);

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

    auto STS_br = linalg::triplet(sBS, Tas, sAR, false, true, false);
    auto STS_as = linalg::triplet(sAR, Tbr, sBS, false, true, false);
    auto STS_ar = linalg::triplet(sAR, Tar, sAR, false, true, false);
    auto STS_bs = linalg::triplet(sBS, Tbs, sBS, false, true, false);

    CompleteInd30 += uAR->vector_dot(omega_ar);
    CompleteInd30 += uBS->vector_dot(omega_bs);
    CompleteInd30 -= STS_br->vector_dot(omega_br);
    CompleteInd30 -= STS_as->vector_dot(omega_as);
    CompleteInd30 -= STS_ar->vector_dot(omega_ar);
    CompleteInd30 -= STS_bs->vector_dot(omega_bs);
    uAR.reset(); uBS.reset(); STS_br.reset(); STS_as.reset(); STS_ar.reset(); STS_bs.reset();

//All Omega-dependent contributions have been completed, now Xi-dependent (J/K) contributions

    auto SCt_Na = linalg::doublet(Ct_Kr, sAR, false, true);
    auto SCt_Nb = linalg::doublet(Ct_Ks, sBS, false, true);

    auto preXiAA = linalg::doublet(SCt_Na, Daa, false, false);
    auto preXiBA = linalg::doublet(SCt_Na, Dab, false, false);
    auto preXiAB = linalg::doublet(SCt_Nb, Dab, false, true);
    auto preXiBB = linalg::doublet(SCt_Nb, Dbb, false, false);

    auto XiAA = linalg::doublet(CoccA_, preXiAA, false, true);
    auto XiAB = linalg::doublet(CoccA_, preXiAB, false, true);
    auto XiBA = linalg::doublet(CoccB_, preXiBA, false, true);
    auto XiBB = linalg::doublet(CoccB_, preXiBB, false, true);

    preXiBB->add(preXiBA);  //enough to calculate JK for the sum of the two

    Cl.clear();
    Cr.clear();
    Cl.push_back(CoccB_);
    Cr.push_back(preXiBB);
    jk_->compute();

    std::shared_ptr<Matrix> J_XiBB = J[0];
    std::shared_ptr<Matrix> K_XiBB = K[0]->transpose();

    CompleteInd30 += 4.0 * XiAA->vector_dot(J_XiBB);
    CompleteInd30 -= 2.0 * XiAA->vector_dot(K_XiBB);
    CompleteInd30 += 4.0 * XiAB->vector_dot(J_XiBB);
    CompleteInd30 -= 2.0 * XiAB->vector_dot(K_XiBB);

    e_exch_ind30_sinf_ = CompleteInd30 - e_ind30_;

    Process::environment.globals["SAPT EXCH-IND30(S^INF) ENERGY"] = e_exch_ind30_sinf_;
    if (print_) {
        outfile->Printf("    Exch-Ind30  (S^inf) = %18.12lf [Eh]\n", e_exch_ind30_sinf_);
        outfile->Printf("\n");
    }
}

}  // namespace sapt
}  // namespace psi
