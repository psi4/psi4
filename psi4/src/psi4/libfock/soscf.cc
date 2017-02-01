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


#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"

#include "psi4/libthce/thce.h"
#include "psi4/libthce/lreri.h"

#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "soscf.h"
#include "jk.h"

namespace psi {

/// SOMCSCF class

SOMCSCF::SOMCSCF(std::shared_ptr<JK> jk, SharedMatrix AOTOSO, SharedMatrix H) :
    jk_(jk)
{
    matrices_["H"] = H;
    matrices_["AOTOSO"] = AOTOSO;
    nao_ = AOTOSO->rowspi()[0];
    casscf_ = true;
    has_fzc_ = false;
    compute_IFock_ = true;
    energy_drc_ = 0.0;
    energy_ci_ = 0.0;
}
SOMCSCF::~SOMCSCF()
{
}
void SOMCSCF::transform(bool approx_only)
{
    throw PSIEXCEPTION("The SOMCSCF object must be initilized as a DF or Disk object.");
}
SharedMatrix SOMCSCF::compute_Q(SharedMatrix TPDM)
{
    throw PSIEXCEPTION("The SOMCSCF object must be initilized as a DF or Disk object.");
}
SharedMatrix SOMCSCF::compute_Qk(SharedMatrix TPDM, SharedMatrix U, SharedMatrix Uact)
{
    throw PSIEXCEPTION("The SOMCSCF object must be initilized as a DF or Disk object.");
}
void SOMCSCF::set_act_MO()
{
    throw PSIEXCEPTION("The SOMCSCF object must be initilized as a DF or Disk object.");
}
void SOMCSCF::set_ras(std::vector<Dimension> ras_spaces)
{
    ras_spaces_ = ras_spaces;
    casscf_ = false;
}
void SOMCSCF::check_ras(void)
{
    // Check the sum of ras equals the act size
    Dimension tot_ras = Dimension(ras_spaces_[0].n(), "Total ras count.");
    for (int i=0; i<ras_spaces_.size(); i++){
        tot_ras += ras_spaces_[i];
    }
    if (tot_ras != nactpi_){
        throw PSIEXCEPTION("SOMSCF: RAS Spaces do not sum up to the total of active spaces\n");
    }

}
void SOMCSCF::set_frozen_orbitals(SharedMatrix Cfzc)
{
    if (Cfzc->ncol()){
        std::vector<SharedMatrix>& Cl = jk_->C_left();
        Cl.clear();
        Cl.push_back(Cfzc);
        jk_->compute();
        Cl.clear();

        const std::vector<SharedMatrix>& J = jk_->J();
        const std::vector<SharedMatrix>& K = jk_->K();

        J[0]->scale(2.0);
        J[0]->subtract(K[0]);

        matrices_["FZC_JK_AO"] = J[0]->clone();
        matrices_["Cfzc"] = Cfzc;

        has_fzc_ = true;
    }
}
void SOMCSCF::set_AO_IFock(SharedMatrix IFock)
{
    matrices_["AO_IFock"] = IFock->clone();
    compute_IFock_ = false;
}
double SOMCSCF::rhf_energy(SharedMatrix C)
{
    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();
    Cl.push_back(C);
    jk_->compute();
    Cl.clear();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    J[0]->scale(2.0);
    J[0]->subtract(K[0]);

    J[0]->add(matrices_["H"]);
    J[0]->add(matrices_["H"]);

    SharedMatrix D = Matrix::doublet(C, C, false, true);
    double erhf = J[0]->vector_dot(D);
    D.reset();
    return erhf;
}
SharedMatrix SOMCSCF::form_rotation_matrix(SharedMatrix x, size_t order) {
    SharedMatrix U(new Matrix("Ck", nirrep_, nmopi_, nmopi_));

    // Form full antisymmetric matrix
    for (size_t h = 0; h < nirrep_; h++) {
        if (!noapi_[h] || !navpi_[h]) continue;
        double** Up = U->pointer(h);
        double** xp = x->pointer(h);

        // Matrix::schmidt orthogonalizes rows not columns so we need to transpose
        for (size_t i = 0, target = 0; i < noapi_[h]; i++) {
            for (size_t a = fmax(noccpi_[h], i); a < nmopi_[h]; a++) {
                Up[i][a] = xp[i][a - noccpi_[h]];
                Up[a][i] = -1.0 * xp[i][a - noccpi_[h]];
            }
        }
    }

    // Build exp(U)
    U->expm(order, true);
    return U;
}

SharedMatrix SOMCSCF::Ck(SharedMatrix C, SharedMatrix x) {
    // C' = C U
    SharedMatrix U = form_rotation_matrix(x);
    SharedMatrix Cp = Matrix::doublet(C, U);

    return Cp;
}
void SOMCSCF::update(SharedMatrix Cocc, SharedMatrix Cact, SharedMatrix Cvir,
            SharedMatrix OPDM, SharedMatrix TPDM)
{
    // => Update orbitals and density matrices <= //
    std::vector<std::shared_ptr<Matrix> > fullC;
    nocc_ = Cocc->ncol();
    noccpi_ = Cocc->colspi();
    matrices_["Cocc"] = Cocc;
    fullC.push_back(Cocc);

    nact_ = Cact->ncol();
    nactpi_ = Cact->colspi();
    matrices_["Cact"] = Cact;
    fullC.push_back(Cact);

    nvir_ = Cvir->ncol();
    nvirpi_ = Cvir->colspi();
    matrices_["Cvir"] = Cvir;
    fullC.push_back(Cvir);

    matrices_["C"] = Matrix::horzcat(fullC);
    matrices_["C"]->set_name("C");

    nirrep_ = matrices_["C"]->nirrep();
    nso_    = matrices_["C"]->nrow();
    nsopi_  = matrices_["C"]->rowspi();
    nmo_    = matrices_["C"]->ncol();
    nmopi_  = matrices_["C"]->colspi();

    noapi_ = noccpi_ + nactpi_;
    navpi_ = nactpi_ + nvirpi_;

    matrices_["OPDM"] = OPDM;
    matrices_["TPDM"] = TPDM;
    set_act_MO();

    // Make sure our ras and active spaces align
    if (!casscf_) check_ras();
    timer_on("SOMCSCF: Compute Gradient");

    // => Build generalized inactive and active Fock matrices <= //
    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();

    Cl.clear();
    Cr.clear();

    // For active Fock
    SharedMatrix CL_COPDM = Matrix::doublet(matrices_["Cact"], matrices_["OPDM"]);
    Cl.push_back(CL_COPDM);
    Cr.push_back(matrices_["Cact"]);

    if (compute_IFock_){
        // For inactive Fock
        Cl.push_back(matrices_["Cocc"]);
        Cr.push_back(matrices_["Cocc"]);
    }

    jk_->compute();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();


    // AFock build
    K[0]->scale(0.5);
    J[0]->subtract(K[0]);
    matrices_["AFock"] = Matrix::triplet(matrices_["C"], J[0], matrices_["C"], true, false, false);
    matrices_["AFock"]->set_name("AFock");

    // IFock build
    if (compute_IFock_) {
        J[1]->scale(2.0);
        J[1]->subtract(K[1]);
        J[1]->add(matrices_["H"]);
        if (has_fzc_) {
            J[1]->add(matrices_["FZC_JK_AO"]);
        }
        matrices_["IFock"] = Matrix::triplet(matrices_["C"], J[1], matrices_["C"],
                                             true, false, false);
        matrices_["IFock"]->set_name("IFock");
    } else {
        matrices_["IFock"] = Matrix::triplet(matrices_["C"], matrices_["AO_IFock"],
                                             matrices_["C"], true, false, false);
        matrices_["IFock"]->set_name("IFock");
    }

    matrices_["Q"] = compute_Q(matrices_["TPDM"]);
    // matrices_["Q"]->print();

    // => Generalized Fock matrix <= //
    matrices_["Fock"] = SharedMatrix(new Matrix("Generalized Fock", nirrep_, nmopi_, nmopi_));
    double *Fp, *IFp, *AFp, *Qp, *OPDMp;
    for (int h=0; h<nirrep_; h++){
        int on = noccpi_[h] * nmopi_[h];
        int an = nactpi_[h] * nmopi_[h];
        if (!on && !an) continue;

        // First index occupied
        if (on){
            IFp = matrices_["IFock"]->pointer(h)[0];
            AFp = matrices_["AFock"]->pointer(h)[0];
            Fp  = matrices_["Fock"] ->pointer(h)[0];
            for (int i=0; i<on; i++){
                Fp[i] = 2.0 * (IFp[i] + AFp[i]);
            }
        }
        // First index active
        if (an){
            IFp = matrices_["IFock"]->pointer(h)[0] + on;
            Fp  = matrices_["Fock"] ->pointer(h)[0] + on;
            OPDMp = matrices_["OPDM"]->pointer(h)[0];

            // OPDM_vw IF_wn => F_vn
            C_DGEMM('N','N',nactpi_[h],nmopi_[h],nactpi_[h],1.0,OPDMp,nactpi_[h],IFp,nmopi_[h],1.0,Fp,nmopi_[h]);

            Qp = matrices_["Q"]->pointer(h)[0];
            for (int i=0; i<an; i++){
                Fp[i] += Qp[i];
            }
        }
    }
    // matrices_["Fock"]->print();

    // => Orbtial Gradient <= //
    matrices_["Gradient"] = SharedMatrix(new Matrix("Gradient", nirrep_, noapi_, navpi_));
    for (int h=0; h<nirrep_; h++){
        if (!noapi_[h] || !navpi_[h]) continue;

        double** dFp = matrices_["Fock"]->pointer(h);
        double** Gp = matrices_["Gradient"]->pointer(h);

        for (int i=0; i<noapi_[h]; i++){
            for (int j=0; j<navpi_[h]; j++){
                int nj = noccpi_[h] + j;
                Gp[i][j] = 2.0 * (dFp[i][nj] - dFp[nj][i]);

                // Ensure the gradient is zero for the active-active diagonal block
                if (nj == i){
                    Gp[i][j] = 0.0;
                }
            }
        }
    }

    // Compute DRC and CI energy energy
    energy_drc_ = 0.0;
    energy_ci_ = 0.0;
    J[0]->add(matrices_["H"]);

    SharedMatrix tmpD = Matrix::doublet(matrices_["Cocc"], matrices_["Cocc"], false, true);
    if (has_fzc_){
        SharedMatrix tmpDf = Matrix::doublet(matrices_["Cfzc"], matrices_["Cfzc"], false, true);
        tmpD->add(tmpDf);
    }
    energy_drc_ = J[0]->vector_dot(tmpD);

    // Compute CI energy
    for (int h = 0; h < nirrep_; h++) {
        for (int t = 0; t < nactpi_[h]; t++) {
            for (int v = 0; v < nactpi_[h]; v++) {
                energy_ci_ +=
                    matrices_["OPDM"]->get(h, t, v) *
                    matrices_["IFock"]->get(h, t + noccpi_[h], v + noccpi_[h]);
            }
        }
    }
    energy_ci_ += 0.5 * matrices_["actMO"]->vector_dot(matrices_["TPDM"]);

    // outfile->Printf("Frozen energy:     %18.12lf\n", energy_fzc_);
    // outfile->Printf("Restricted energy: %18.12lf\n", energy_drc_);
    // outfile->Printf("CI energy:         %18.12lf\n", energy_ci_);

    zero_redundant(matrices_["Gradient"]);
    timer_off("SOMCSCF: Compute Gradient");
    matrices_["Precon"] = H_approx_diag();
    // matrices_["Gradient"]->print();
    // matrices_["Precon"]->print();
}
SharedMatrix SOMCSCF::H_approx_diag()
{

    timer_on("SOMCSCF: Approximate hessian");
    // d_tuvw I_tuvw -> dI_t, D_tu IF_tu -> DIF_t
    double** actMOp = matrices_["actMO"]->pointer();
    double** TPDMp = matrices_["TPDM"]->pointer();
    int relact = 0;
    int nact3 = nact_*nact_*nact_;

    SharedVector dI(new Vector("dI", nactpi_));
    SharedVector DIF(new Vector("IF * OPDM", nactpi_));
    for (int h=0; h<nirrep_; h++){
        if (!nactpi_[h]) continue;

        double* dIp = dI->pointer(h);
        double* DIFp = DIF->pointer(h);
        double** IFp = matrices_["IFock"]->pointer(h);
        double** OPDMp = matrices_["OPDM"]->pointer(h);
        for (int a=0; a<nactpi_[h]; a++){
            size_t shift = nact3 * (a + relact);
            dIp[a] = C_DDOT(nact3, TPDMp[0]+shift, 1, actMOp[0]+shift, 1);
            DIFp[a] = C_DDOT(nactpi_[h], OPDMp[a], 1, IFp[noccpi_[h] + a]+noccpi_[h], 1);
        }
        relact += nactpi_[h];
    }

    SharedMatrix H(new Matrix("Approximate diag hessian", nirrep_, noapi_, navpi_));
    int offset_act = 0;
    for (int h=0; h<nirrep_; h++){
        if (!noapi_[h] || !navpi_[h]) continue;
        double** Hp = H->pointer(h);
        double** IFp = matrices_["IFock"]->pointer(h);
        double** AFp = matrices_["AFock"]->pointer(h);

        // iv block
        if (noccpi_[h] && nvirpi_[h]){
            for (int i=0; i<noccpi_[h]; i++){
                double tmp_i = IFp[i][i] + AFp[i][i];
                for (int v=0; v<nvirpi_[h]; v++){
                    int nv = noapi_[h] + v;
                    Hp[i][nactpi_[h] + v] = 4.0 * (IFp[nv][nv] + AFp[nv][nv] - tmp_i);
                }
            }
        } // end iv block

        // av block
        if (nactpi_[h] && nvirpi_[h]){
            double** OPDMp = matrices_["OPDM"]->pointer(h);
            double* dIp = dI->pointer(h);
            double* DIFp = DIF->pointer(h);
            for (int a=0; a<nactpi_[h]; a++){
                int oa = noccpi_[h] + a;
                for (int v=0; v<nvirpi_[h]; v++){
                    int oav = noapi_[h] + v;
                    int av = nactpi_[h] + v;
                    Hp[oa][av]  = OPDMp[a][a] * IFp[oav][oav];
                    Hp[oa][av] -= dIp[a];
                    Hp[oa][av] -= DIFp[a];
                    Hp[oa][av] += OPDMp[a][a] * AFp[oav][oav];
                    Hp[oa][av] *= 2.0;
                }
            }
        } // end av block

        // ia block
        if (nactpi_[h] && noccpi_[h]){
            double** OPDMp = matrices_["OPDM"]->pointer(h);
            double* dIp = dI->pointer(h);
            double* DIFp = DIF->pointer(h);
            for (int i=0; i<noccpi_[h]; i++){
                for (int a=0; a<nactpi_[h]; a++){
                    int oa = noccpi_[h] + a;
                    Hp[i][a]  = 2.0 * (IFp[oa][oa] + AFp[oa][oa]);
                    Hp[i][a] -= 2.0 * (IFp[i][i] + AFp[i][i]);
                    Hp[i][a] += (OPDMp[a][a] * IFp[i][i]);
                    Hp[i][a] += (OPDMp[a][a] * AFp[i][i]);
                    Hp[i][a] -= dIp[a];
                    Hp[i][a] -= DIFp[a];
                    Hp[i][a] *= 2.0;
                }
            }
        } // end ia block

        // aa block is one to prevent divide by zero
        if (nactpi_[h]){
            for (int i=0; i<nactpi_[h]; i++){
                for (int a=0; a<nactpi_[h]; a++){
                    Hp[i + noccpi_[h]][a] = 1.0;
                }
            }
        }// end aa block

        // aa block for RASSCF
        if (nactpi_[h] && !casscf_){
            double* dIp = dI->pointer(h);
            double* DIFp = DIF->pointer(h);
            double** IFp = matrices_["IFock"]->pointer(h);
            double** OPDMp = matrices_["OPDM"]->pointer(h);

            // We need to do these completely flat, some compilers vectorize this incorrectly
            double* actMO_fp = matrices_["actMO"]->pointer()[0];
            double* TPDM_fp = matrices_["TPDM"]->pointer()[0];

            int nact2 = nact_ * nact_;
            int offset_col = 0;
            int offset_row = 0;

            // Loop over spaces, last space will have no rotations
            for (int nras = 0; nras < ras_spaces_.size()-1; nras++){
                int ras_size = ras_spaces_[nras][h];
                offset_col += ras_size;

                // Loop over pairs
                for (int i=offset_row; i<(offset_row + ras_size); i++){
                for (int a=offset_col; a<nactpi_[h]; a++){
                    int oa = a + noccpi_[h];
                    int oi = i + noccpi_[h];

                    double value = 0.0;
                    value += 2.0 * OPDMp[a][a] * IFp[oi][oi];
                    value += 2.0 * OPDMp[i][i] * IFp[oa][oa];
                    value -= 4.0 * OPDMp[i][a] * IFp[oi][oa];
                    Hp[oi][a] = value;

                    value = 0.0;
                    value += (dIp[i] + dIp[a]);
                    value += DIFp[i];
                    value += DIFp[a];
                    Hp[oi][a] -= 2.0 * value;


                    value = 0.0;
                    int p = i + offset_act;
                    int q = a + offset_act;
                    int pp = p * nact_ + p;
                    int pq = p * nact_ + q;
                    int qq = q * nact_ + q;
                    for (int u=0; u<nact_; u++){
                        int pu = p * nact_ + u;
                        int qu = q * nact_ + u;

                        // Unrolling this to prevent some strange vectorization error
                        for (int v=0; v<nact_; v++){
                            int pv = p * nact_ + v;
                            int qv = q * nact_ + v;
                            int uv = u * nact_ + v;
                            value += 2.0 * TPDM_fp[pu * nact2 + pv] * actMO_fp[qu * nact2 + qv];
                            value += 2.0 * TPDM_fp[qu * nact2 + qv] * actMO_fp[pu * nact2 + pv];
                            value +=       TPDM_fp[pp * nact2 + uv] * actMO_fp[qq * nact2 + uv];
                            value +=       TPDM_fp[qq * nact2 + uv] * actMO_fp[pp * nact2 + uv];
                            value -= 4.0 * TPDM_fp[pv * nact2 + qu] * actMO_fp[pu * nact2 + qv];
                            value -= 2.0 * TPDM_fp[pq * nact2 + uv] * actMO_fp[pq * nact2 + uv];
                        }
                    }
                    // outfile->Printf("%d %d | %d %d : %lf %lf\n", h, nras, i, a, Hp[oi][a], value);
                    Hp[oi][a] += 2.0 * value;

                } // End i loop
                } // End a loop
                offset_row += ras_size;
            } // End ras loop
            offset_act += nactpi_[h];
        } // End active-active RASSCF block

    } // End irrep loop
    timer_off("SOMCSCF: Approximate hessian");
    return H;
}
SharedMatrix SOMCSCF::compute_AFock(SharedMatrix OPDM){
    // => Build generalized inactive and active Fock matrices <= //
    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();

    Cl.clear();
    Cr.clear();

    // For active Fock
    SharedMatrix CL_COPDM = Matrix::doublet(matrices_["Cact"], OPDM);
    Cl.push_back(CL_COPDM);
    Cr.push_back(matrices_["Cact"]);

    jk_->compute();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    // AFock build
    K[0]->scale(0.5);
    J[0]->subtract(K[0]);
    SharedMatrix AFock = Matrix::triplet(matrices_["C"], J[0], matrices_["C"], true, false, false);
    AFock->set_name("AFock");
    return AFock;
}
SharedMatrix SOMCSCF::Hk(SharedMatrix x)
{
    timer_on("SOMCSCF: Rotated fock");

    // => Antisymmetric rotation matrix <= //
    SharedMatrix U(new Matrix("U", nirrep_, nmopi_, nmopi_));
    SharedMatrix Uocc(new Matrix("Uocc", nirrep_, noccpi_, nmopi_));
    SharedMatrix Uact(new Matrix("Uact", nirrep_, nactpi_, nmopi_));
    for (int h=0; h<nirrep_; h++){

        if (!noapi_[h] || !navpi_[h]) continue;
        double** Up = U->pointer(h);
        double** xp = x->pointer(h);

        for (int i=0; i<noapi_[h]; i++){
            for (int a=0; a < navpi_[h]; a++){
                int offa = noccpi_[h] + a;
                Up[i][offa] = xp[i][a];
                Up[offa][i] = -1.0 * xp[i][a];
            }
        }
        // Fill Uocc
        if (noccpi_[h]){
            double** Uoccp = Uocc->pointer(h);
            for (int i=0; i<noccpi_[h]; i++){
                for (int j=0; j<nmopi_[h]; j++){
                    Uoccp[i][j] = Up[i][j];
                }
            }
        }
        // Fill Ua
        if (nactpi_[h]){
            double** Uactp = Uact->pointer(h);
            for (int i=0; i<nactpi_[h]; i++){
                for (int j=0; j<nmopi_[h]; j++){
                    Uactp[i][j] = Up[i + noccpi_[h]][j];
                }
            }
        }
    }

    // => Rotated inactive and active Fock matrices <= //
    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();

    // For inactive Fock

    SharedMatrix CLUocc = Matrix::doublet(matrices_["C"], Uocc, false, true);
    Cl.push_back(CLUocc);
    Cr.push_back(matrices_["Cocc"]);

    // For active Fock
    SharedMatrix CLUact = Matrix::triplet(matrices_["C"], Uact, matrices_["OPDM"], false, true, true);
    Cr.push_back(CLUact);
    Cl.push_back(matrices_["Cact"]);

    jk_->compute();
    Cl.clear();
    Cr.clear();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    // Rotated inactive fock
    SharedMatrix IFk = Matrix::doublet(matrices_["IFock"], U, false, true);
    IFk->gemm(false, false, 1.0, U, matrices_["IFock"], 1.0);

    J[0]->scale(4.0);
    J[0]->subtract(K[0]);
    J[0]->subtract(K[0]->transpose());
    SharedMatrix trans_half = Matrix::doublet(J[0], matrices_["C"]);
    IFk->gemm(true, false, 1.0, matrices_["C"], trans_half, 1.0);

    // Rotated active fock
    SharedMatrix ret = Matrix::doublet(matrices_["AFock"], U, false, true);
    ret->gemm(false, false, 1.0, U, matrices_["AFock"], 1.0);

    J[1]->scale(2.0);
    K[1]->scale(0.5);
    J[1]->subtract(K[1]);
    J[1]->subtract(K[1]->transpose());
    trans_half = Matrix::doublet(J[1], matrices_["C"]);
    ret->gemm(true, false, 1.0, matrices_["C"], trans_half, 1.0);

    trans_half.reset();
    ret->add(IFk);
    ret->scale(2.0);

    /// Build Qk
    matrices_["Qk"] = compute_Qk(matrices_["TPDM"], U, Uact);
    // outfile->Printf("Active Fock\n");
    // matrices_["Qk"]->print();

    // Add in Q and zero out virtual
    for (int h=0; h<nirrep_; h++){
        if (nactpi_[h]) {
            double** Fkp = ret->pointer(h);
            double** Qkp = matrices_["Qk"]->pointer(h);

            // OPDM_vw IF_wn->vn
            C_DGEMM('N','N',nactpi_[h],nmopi_[h],nactpi_[h],1.0,
                    matrices_["OPDM"]->pointer(h)[0],nactpi_[h],
                    IFk->pointer(h)[noccpi_[h]],nmopi_[h],
                    0.0,Fkp[noccpi_[h]],nmopi_[h]);

            // OPDM_vw += Qk
            C_DAXPY(nmopi_[h]*nactpi_[h], 1.0, Qkp[0], 1, Fkp[noccpi_[h]], 1);
        }

        if (nvirpi_[h]){
            double** Fkp = ret->pointer(h);
            // Zero out the Fk[a,n] part
            for (int i=noapi_[h]; i<nmopi_[h]; i++){
                for (int j=0; j<nmopi_[h]; j++){
                    Fkp[i][j] = 0.0;
                }
            }
        }
    }


    // => Hessian <= //
    SharedMatrix hessx(new Matrix("Hessian x", nirrep_, noapi_, navpi_));

    for (int h=0; h<nirrep_; h++){
        if (!noapi_[h] || !navpi_[h]) continue;

        double** Hxp = hessx->pointer(h);
        double** Fkp = ret->pointer(h);

        for (int i=0; i<noapi_[h]; i++){
            for (int j=0; j<navpi_[h]; j++){
                int nj = noccpi_[h] + j;
                Hxp[i][j] = 2.0 * (Fkp[i][nj] - Fkp[nj][i]);
            }
        }
    }

    zero_redundant(hessx);
    timer_off("SOMCSCF: Rotated fock");
    return hessx;

}
SharedMatrix SOMCSCF::approx_solve()
{
    // outfile->Printf("In approx solve\n");

    SharedMatrix ret = matrices_["Gradient"]->clone();
    ret->apply_denominator(matrices_["Precon"]);
    zero_redundant(ret);
    return ret;
}
SharedMatrix SOMCSCF::solve(int max_iter, double conv, bool print)
{
    if (print){
        outfile->Printf("\n");
        outfile->Printf("    ==> SOMCSCF Iterations <==\n");
        outfile->Printf("    Maxiter     = %11d\n", max_iter);
        outfile->Printf("    Convergence = %11.3E\n", conv);
        outfile->Printf("    ---------------------------------------\n");
        outfile->Printf("    %-4s   %11s     %10s\n", "Iter", "Residual RMS", "Time [s]");
        outfile->Printf("    ---------------------------------------\n");
    }

    time_t start;
    time_t stop;
    start = time(NULL);

    // Initial guess
    SharedMatrix x = matrices_["Gradient"]->clone();
    x->set_name("Trial Vector x");
    x->apply_denominator(matrices_["Precon"]);

    // Calc hessian vector product, find residual and conditioned residual
    SharedMatrix r = matrices_["Gradient"]->clone();


    SharedMatrix Ap = Hk(x);
    // outfile->Printf("Ap\n");
    // Ap->print();

    // outfile->Printf("Gradient\n");
    // r->print();
    // outfile->Printf("Ap guess\n");
    // Ap->print();

    r->subtract(Ap);
    if (print){
        double rconv = r->rms();
        stop = time(NULL);
        outfile->Printf("    %-4d %11.3E %10ld\n", 0, rconv, stop-start);
    }

    SharedMatrix z = r->clone();
    z->apply_denominator(matrices_["Precon"]);

    SharedMatrix p = z->clone();

    SharedMatrix best = x->clone();
    double best_conv = r->rms();
    for (int iter = 0; iter < max_iter; iter++) {
        // Calc hessian vector product
        Ap = Hk(p);

        // Find factors and scale
        double rzpre = r->vector_dot(z);
        double alpha = rzpre / p->vector_dot(Ap);

        x->axpy(alpha, p);
        r->axpy(-alpha, Ap);

        // Get residual
        double rconv = r->rms();
        stop = time(NULL);
        if (print){
            outfile->Printf("    %-4d %11.3E %10ld\n", iter+1, rconv, stop-start);
        }

        // Convergence may not be monotonic
        if (rconv < best_conv){
            best_conv = rconv;
            best->copy(x);
        }

        // Check convergence
        if (rconv < conv){
            break;
        }

        // Update p and z
        z->copy(r);
        z->apply_denominator(matrices_["Precon"]);

        double beta = r->vector_dot(z) / rzpre;
        // if (beta > 0.15) beta = 0.15;
        p->scale(beta);
        p->add(z);

    }// End iterations

    if (print){
        outfile->Printf("    %-4s %11.3E %10s\n", "Best", best_conv, "--");
        outfile->Printf("    ---------------------------------------\n");
        outfile->Printf("\n");
    }
    zero_redundant(best);
    return best;
}
SharedMatrix SOMCSCF::gradient()
{
   return matrices_["Gradient"];
}
double SOMCSCF::gradient_rms()
{
   return matrices_["Gradient"]->rms();
}
void SOMCSCF::zero_redundant(SharedMatrix vector){
    if (casscf_){
        zero_act(vector);
    }
    else{
        zero_ras(vector);
    }
}
void SOMCSCF::zero_act(SharedMatrix vector){
    for (int h=0; h<nirrep_; h++){
        if (!nactpi_[h]) continue;

        double** vp = vector->pointer(h);
        for (int i=0; i<nactpi_[h]; i++){
            for (int j=0; j<nactpi_[h]; j++){
                int io = noccpi_[h] + i;
                vp[io][j] = 0.0;
            }
        }
    }
}
void SOMCSCF::zero_ras(SharedMatrix vector){
    for (int h=0; h<nirrep_; h++){
        if (!nactpi_[h]) continue;

        double** vp = vector->pointer(h);

        // Loop over spaces last space will have no rotations
        int offset_row = 0;
        int offset_col = 0;
        for (int nras = 0; nras < ras_spaces_.size(); nras++){
            int ras_size = ras_spaces_[nras][h];
            offset_col += ras_size;

            // Loop over pairs
            for (int i=offset_row; i<(offset_row + ras_size); i++){
                for (int a=0; a<offset_col; a++){
                    vp[noccpi_[h]+i][a] = 0.0;
                }
            }
            offset_row += ras_size;
        } // End ras loop
    }
}
/// End SOMCSCF class

/// DFSOMCSCF class
DFSOMCSCF::DFSOMCSCF(std::shared_ptr<JK> jk, std::shared_ptr<DFERI> df, SharedMatrix AOTOSO,
            SharedMatrix H) :
            SOMCSCF(jk, AOTOSO, H)
{
    dferi_ = df;
}
DFSOMCSCF::~DFSOMCSCF()
{
}
void DFSOMCSCF::transform(bool approx_only)
{
    // => AO C matrices <= //
    // We want a pitzer order C matrix with appended Cact
    SharedMatrix Cocc = matrices_["Cocc"];
    SharedMatrix Cact = matrices_["Cact"];
    SharedMatrix Cvir = matrices_["Cvir"];

    int nrot = Cocc->ncol() + Cact->ncol() + Cvir->ncol();
    int aoc_rowdim =  nrot + Cact->ncol();
    SharedMatrix AO_C = SharedMatrix(new Matrix("AO_C", nao_, aoc_rowdim));

    double** AO_Cp = AO_C->pointer();
    for (int h=0, offset=0, offset_act=0; h < nirrep_; h++){
        int hnso = nsopi_[h];
        if (hnso == 0) continue;
        double** Up = matrices_["AOTOSO"]->pointer(h);

        int noccpih = Cocc->colspi()[h];
        int nactpih = Cact->colspi()[h];
        int nvirpih = Cvir->colspi()[h];
        // occupied
        if (noccpih){
            double** CSOp = Cocc->pointer(h);
            C_DGEMM('N','N',nao_,noccpih,hnso,1.0,Up[0],hnso,CSOp[0],noccpih,0.0,&AO_Cp[0][offset],aoc_rowdim);
            offset += noccpih;
        }
        // active
        if (nactpih){
            double** CSOp = Cact->pointer(h);
            C_DGEMM('N','N',nao_,nactpih,hnso,1.0,Up[0],hnso,CSOp[0],nactpih,0.0,&AO_Cp[0][offset],aoc_rowdim);
            offset += nactpih;

            C_DGEMM('N','N',nao_,nactpih,hnso,1.0,Up[0],hnso,CSOp[0],nactpih,0.0,&AO_Cp[0][offset_act + nrot],aoc_rowdim);
            offset_act += nactpih;
        }
        // virtual
        if (nvirpih){
            double** CSOp = Cvir->pointer(h);
            C_DGEMM('N','N',nao_,nvirpih,hnso,1.0,Up[0],hnso,CSOp[0],nvirpih,0.0,&AO_Cp[0][offset],aoc_rowdim);
            offset += nvirpih;
        }
    }


    // => Compute DF ints <= //
    dferi_->clear();
    dferi_->set_C(AO_C);
    dferi_->add_space("R", 0, nrot);
    dferi_->add_space("a", nrot, aoc_rowdim);
    dferi_->add_space("F", 0, aoc_rowdim);

    if (approx_only){
      dferi_->add_pair_space("aaQ", "a", "a");
      dferi_->add_pair_space("RaQ", "a", "R", -1.0/2.0, true);
    }
    else{
      dferi_->add_pair_space("aaQ", "a", "a");
      dferi_->add_pair_space("RaQ", "a", "R", -1.0/2.0, true);
      dferi_->add_pair_space("RRQ", "R", "R");
    }
    dferi_->compute();
}
void DFSOMCSCF::set_act_MO()
{
    // Build (aa|aa)
    std::map<std::string, std::shared_ptr<Tensor> >& dfints = dferi_->ints();
    int nQ = dferi_->size_Q();
    std::shared_ptr<Tensor> aaQT = dfints["aaQ"];
    SharedMatrix aaQ(new Matrix("aaQ", nact_ * nact_, nQ));
    double* aaQp = aaQ->pointer()[0];
    FILE* aaQF = aaQT->file_pointer();
    fseek(aaQF,0L,SEEK_SET);
    fread(aaQp, sizeof(double), nact_ * nact_ * nQ, aaQF);
    matrices_["actMO"] = Matrix::doublet(aaQ, aaQ, false, true);
    aaQ.reset();
}
SharedMatrix DFSOMCSCF::compute_Q(SharedMatrix TPDM){

    timer_on("SOMCSCF: DF-Q matrix");
    std::map<std::string, std::shared_ptr<Tensor> >& dfints = dferi_->ints();

    int nQ = dferi_->size_Q();
    int nact2 = nact_ * nact_;
    double* TPDMp = TPDM->pointer()[0];

    // Load aaQ
    std::shared_ptr<Tensor> aaQT = dfints["aaQ"];
    SharedMatrix aaQ(new Matrix("aaQ", nact_ * nact_, nQ));
    double* aaQp = aaQ->pointer()[0];
    FILE* aaQF = aaQT->file_pointer();
    fseek(aaQF,0L,SEEK_SET);
    fread(aaQp, sizeof(double), nact_ * nact_ * nQ, aaQF);

    // d_vwxy I_xyQ -> d_vwQ (Qa^4)
    SharedMatrix vwQ(new Matrix("vwQ", nact_ * nact_, nQ));
    double* vwQp = vwQ->pointer()[0];
    C_DGEMM('N','N',nact2,nQ,nact2,1.0,TPDMp,nact2,aaQp,nQ,0.0,vwQp,nQ);
    aaQ.reset();

    // Load NaQ
    std::shared_ptr<Tensor> NaQT = dfints["RaQ"];
    SharedMatrix NaQ(new Matrix("NaQ", nmo_ * nact_, nQ));
    double* NaQp = NaQ->pointer()[0];
    FILE* NaQF = NaQT->file_pointer();
    fseek(NaQF,0L,SEEK_SET);
    fread(NaQp, sizeof(double), nmo_ * nact_ * nQ, NaQF);

    // d_vwQ I_NwQ -> Q_vN (NQa^2)
    SharedMatrix denQ(new Matrix("Dense Qvn", nact_, nmo_));
    double** denQp = denQ->pointer();
    C_DGEMM('N','T',nact_,nmo_,nQ*nact_,1.0,vwQp,nQ*nact_,NaQp,nQ*nact_,0.0,denQp[0],nmo_);
    NaQ.reset();

    // Symmetry block Q
    SharedMatrix Q(new Matrix("Qvn", nirrep_, nactpi_, nmopi_));

    int offset_act = 0;
    int offset_nmo = 0;
    for (int h=0; h<nirrep_; h++){
        if (!nactpi_[h] || !nmopi_[h]){
            offset_nmo += nmopi_[h];
            continue;
        }

        double* Qp = Q->pointer(h)[0];
        for (int i=0, target=0; i<nactpi_[h]; i++){
            for (int j=0; j<nmopi_[h]; j++){
                Qp[target++] = denQp[offset_act + i][offset_nmo + j];
            }
        }
        offset_act += nactpi_[h];
        offset_nmo += nmopi_[h];
    }
    timer_off("SOMCSCF: DF-Q matrix");
    // outfile->Printf("Printing the Q matrix\n");
    // matrices_["Q"]->print();
    return Q;
}
SharedMatrix DFSOMCSCF::compute_Qk(SharedMatrix TPDM, SharedMatrix U, SharedMatrix Uact)
{
    timer_on("SOMCSCF: DF-Qk matrix");
    // Remove U symmetry
    SharedMatrix dUact;
    if (nirrep_ == 1){
        dUact = Uact;
    }
    else{
        dUact = Uact->to_block_sharedmatrix();
    }

    int nQ = dferi_->size_Q();
    int nact2 = nact_*nact_;
    int nact3 = nact2*nact_;
    double* TPDMp = TPDM->pointer()[0];
    std::map<std::string, std::shared_ptr<Tensor> >& dfints = dferi_->ints();

    // Read NaQ
    std::shared_ptr<Tensor> NaQT = dfints["RaQ"];
    SharedMatrix NaQ(new Matrix("RaQ", nmo_ * nact_, nQ));
    double* NaQp = NaQ->pointer()[0];
    FILE* NaQF = NaQT->file_pointer();
    fseek(NaQF,0L,SEEK_SET);
    fread(NaQp, sizeof(double), nmo_ * nact_ * nQ, NaQF);

    SharedMatrix xyQ(new Matrix("xyQ", nact_ * nact_, nQ));
    double* xyQp = xyQ->pointer()[0];
    double** dUactp = dUact->pointer();

    // oyQ,xo->xyQ (NaQ, aU) QNa^2
    C_DGEMM('N','N',nact_,nQ*nact_,nmo_,1.0,dUactp[0],nmo_,NaQp,nQ*nact_,0.0,xyQp,nQ*nact_);

    // xyQ += yxQ
    for (int x=0; x<nact_; x++){
        for (int y=x; y<nact_; y++){
            for (int q=0; q<nQ; q++){
                size_t idx1 = x * nQ * nact_ + y * nQ + q;
                size_t idx2 = y * nQ * nact_ + x * nQ + q;
                double tmpn = xyQp[idx1] + xyQp[idx2];
                xyQp[idx1] = tmpn;
                xyQp[idx2] = tmpn;
            }
        }
    }

    // nwQ,xyQ => tmp_nwxy (NaQ, xyQ) NQa^3
    SharedMatrix Gnwxy = Matrix::doublet(NaQ, xyQ, false, true);
    double* Gnwxyp = Gnwxy->pointer()[0];
    NaQ.reset();

    // vwxy,nwxy => Qk_vm (TPDM, tmp_nwxy) Na^4
    SharedMatrix dQk(new Matrix("dQk", nact_, nmo_));
    double** dQkp = dQk->pointer();

    C_DGEMM('N','T',nact_,nmo_,nact3,1.0,TPDMp,nact3,Gnwxyp,nact3,0.0,dQkp[0],nmo_);

    // wo,onQ->wnQ (aU, NNQ) QN^2a^2 (rate determining)
    // Read and gemm in chunks, need to do blocking later
    int chunk_size = 200;
    std::shared_ptr<Tensor> NNQT = dfints["RRQ"];

    SharedMatrix wnQ(new Matrix("nwQ", nact_*nmo_, nQ));
    double** wnQp = wnQ->pointer();
    SharedMatrix NNQ(new Matrix("RRQ", chunk_size * nmo_, nQ));
    double** NNQp = NNQ->pointer();
    FILE* NNQF = NNQT->file_pointer();

    for (int start=0; start < nmo_; start += chunk_size){
        int block = (start + chunk_size > nmo_ ? nmo_ - start : chunk_size);
        size_t read_start = sizeof(double) * start * nmo_ * nQ;
        fseek(NNQF,read_start,SEEK_SET);
        fread(NNQp[0], sizeof(double), block * nmo_ * nQ, NNQF);
        C_DGEMM('N','N',nact_,nmo_*nQ,block,1.0,
                dUactp[0]+start,nmo_,
                NNQp[0],nmo_*nQ,
                1.0,wnQp[0],nmo_*nQ);
    }
    NNQ.reset();

    // Read aaQ
    std::shared_ptr<Tensor> aaQT = dfints["aaQ"];
    SharedMatrix aaQ(new Matrix("aaQ", nact_ * nact_, nQ));
    double* aaQp = aaQ->pointer()[0];
    FILE* aaQF = aaQT->file_pointer();
    fseek(aaQF,0L,SEEK_SET);
    fread(aaQp, sizeof(double), nact_ * nact_ * nQ, aaQF);

    // wnQ,xyQ => tmp_wnxy (wNQ, aaQ) NQa^3
    C_DGEMM('N','T',nmo_*nact_,nact2,nQ,1.0,wnQp[0],nQ,aaQp,nQ,0.0,Gnwxyp,nact2);
    aaQ.reset();

    // Should probably figure out an inplace algorithm
    SharedMatrix Gleft(new Matrix("Gnwxy", nmo_, nact3));
    double* Gleftp = Gleft->pointer()[0];

    // wnxy => nwxy
    size_t target = 0;
    for (int n=0; n<nmo_; n++){
    for (int w=0; w<nact_; w++){
    for (int x=0; x<nact_; x++){
    for (int y=0; y<nact_; y++){
        Gleftp[target++] = Gnwxyp[w*nmo_*nact2 + n*nact2 + x*nact_ + y];
    }}}}
    Gnwxy.reset();

    // vwxy,nwxy => Qk_vm (TPDM, tmp_nwxy) Na^4
    C_DGEMM('N','T',nact_,nmo_,nact3,1.0,TPDMp,nact3,Gleftp,nact3,1.0,dQkp[0],nmo_);

    // Symm block Qk
    SharedMatrix tQ = compute_Q(TPDM);
    SharedMatrix Qk = Matrix::doublet(tQ, U, false, true);
    // dQk->print();

    int offset_act = 0;
    int offset_nmo = 0;
    for (int h=0; h<nirrep_; h++){
        if (!nactpi_[h]){
            offset_nmo += nmopi_[h];
            continue;
        }
        double** Qkp = Qk->pointer(h);
        for (int a=0; a<nactpi_[h]; a++){
            C_DAXPY(nmopi_[h], 1.0, dQkp[offset_act+a]+offset_nmo, 1, Qkp[a], 1);
        }
        offset_act += nactpi_[h];
        offset_nmo += nmopi_[h];
    }
    // matrices_["Qk"]->print();
    timer_off("SOMCSCF: DF-Qk matrix");

    return Qk;

}// End DFSOMCSCF object

/// DiskSOMCSCF class
DiskSOMCSCF::DiskSOMCSCF(std::shared_ptr<JK> jk,
            std::shared_ptr<IntegralTransform> ints,
            SharedMatrix AOTOSO, SharedMatrix H) :
            SOMCSCF(jk, AOTOSO, H)
{
    ints_ = ints;
    psio_ = _default_psio_lib_;
}
DiskSOMCSCF::~DiskSOMCSCF()
{
}
void DiskSOMCSCF::transform(bool approx_only)
{
    throw PSIEXCEPTION("DiskSOMCSCF::transform is not supported for Disk integrals.");
}
void DiskSOMCSCF::set_act_MO()
{
    dpdbuf4 I;

    // => Read dense active MO <= //
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0,
                ints_->DPD_ID("[X>=X]+"), ints_->DPD_ID("[X>=X]+"),
                ints_->DPD_ID("[X>=X]+"), ints_->DPD_ID("[X>=X]+"), 0, "MO Ints (XX|XX)");

    matrices_["actMO"] = SharedMatrix(new Matrix("actMO", nact_*nact_, nact_*nact_));
    double** actMOp = matrices_["actMO"]->pointer();

    for (int h=0; h<nirrep_; h++){
        global_dpd_->buf4_mat_irrep_init(&I, h);
        global_dpd_->buf4_mat_irrep_rd(&I, h);
    }

    // 8 fold symmetry
    for(int p = 0; p < nact_; p++){
      int p_sym = I.params->psym[p];

      for(int q = 0; q <= p; q++){
        int q_sym = I.params->qsym[q];
        int pq = I.params->rowidx[p][q];
        int pq_sym = p_sym^q_sym;

        for(int r = 0; r <= p; r++){
          int r_sym = I.params->rsym[r];
          int smax = (p==r) ? q+1 : r+1;

          for(int s = 0; s < smax; s++){
            int s_sym = I.params->ssym[s];
            int rs_sym = r_sym^s_sym;
            int rs = I.params->colidx[r][s];

            if (pq_sym != rs_sym) continue;

            double value = I.matrix[pq_sym][pq][rs];

            actMOp[p*nact_ + q][r*nact_ + s] = value;
            actMOp[q*nact_ + p][r*nact_ + s] = value;
            actMOp[p*nact_ + q][s*nact_ + r] = value;
            actMOp[q*nact_ + p][s*nact_ + r] = value;

            actMOp[r*nact_ + s][p*nact_ + q] = value;
            actMOp[s*nact_ + r][p*nact_ + q] = value;
            actMOp[r*nact_ + s][q*nact_ + p] = value;
            actMOp[s*nact_ + r][q*nact_ + p] = value;

          }
        }
      }
    }

    // Close everything out
    for (int h=0; h<nirrep_; h++){
        global_dpd_->buf4_mat_irrep_close(&I, h);
    }
    global_dpd_->buf4_close(&I);
    psio_->close(PSIF_LIBTRANS_DPD, 1);
}
SharedMatrix DiskSOMCSCF::compute_Q(SharedMatrix TPDMmat)
{
    timer_on("SOMCSCF: Q matrix");

    // => Write active TPDM <= //
    dpdbuf4 G, TPDM;
    dpdfile2 Q;

    double** TPDMmatp = TPDMmat->pointer();
    psio_->open(PSIF_MCSCF, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&TPDM, PSIF_MCSCF, 0,
                ints_->DPD_ID("[X,X]"), ints_->DPD_ID("[X,X]"),
                ints_->DPD_ID("[X>=X]+"), ints_->DPD_ID("[X>=X]+"), 0, "CI TPDM (XX|XX)");

    for (int h=0; h<nirrep_; h++){
        global_dpd_->buf4_mat_irrep_init(&TPDM, h);
    }

    // 4 fold symmetry
    for(int p = 0; p < nact_; p++){
      int p_sym = TPDM.params->psym[p];

      for(int q = 0; q <= p; q++){
        int q_sym = TPDM.params->psym[q];
        int pq_sym = p_sym^q_sym;
        int pq = TPDM.params->rowidx[p][q];

        for(int r = 0; r < nact_; r++){
          int r_sym = TPDM.params->psym[r];

          for(int s = 0; s <= r; s++){
            int s_sym = TPDM.params->psym[s];
            int rs_sym = r_sym^s_sym;
            int rs = TPDM.params->colidx[r][s];

            if (pq_sym != rs_sym) continue;

            TPDM.matrix[pq_sym][pq][rs] = TPDMmatp[p*nact_ + q][r*nact_ + s];
    }}}}

    for (int h=0; h<nirrep_; h++){
        global_dpd_->buf4_mat_irrep_wrt(&TPDM, h);
        global_dpd_->buf4_mat_irrep_close(&TPDM, h);
    }

    // G_mwxy TPDM_vwxy -> Q_mv
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // We init R then A so we need 1, 0 here for index types
    global_dpd_->file2_init(&Q, PSIF_MCSCF, 0, 1, 0, "Q");

    global_dpd_->buf4_init(&G, PSIF_LIBTRANS_DPD, 0,
                ints_->DPD_ID("[X,X]"), ints_->DPD_ID("[X,R]"),
                ints_->DPD_ID("[X>=X]+"), ints_->DPD_ID("[X,R]"), 0, "MO Ints (XX|XR)");

    global_dpd_->contract442(&TPDM, &G, &Q, 3, 3, 1.0, 0.0);
    SharedMatrix Qmat(new Matrix(&Q));

    global_dpd_->file2_close(&Q);
    global_dpd_->buf4_close(&TPDM);
    global_dpd_->buf4_close(&G);

    psio_->close(PSIF_LIBTRANS_DPD, 1);
    psio_->close(PSIF_MCSCF, 1);

    timer_off("SOMCSCF: Q matrix");

    return Qmat;

}
SharedMatrix DiskSOMCSCF::compute_Qk(SharedMatrix TPDMmat, SharedMatrix U, SharedMatrix Uact)
{
    timer_on("SOMCSCF: Qk matrix");
    // \TPDM_{vwxy}\kappa_{mo}g_{owxy}
    // \TPDM_{vwxy}(\kappa_{wo}g_{moxy} +\kappa_{xo}g_{mwoy} + \kappa_{yo}g_{mwxo})

    bool debug = false;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    psio_->open(PSIF_MCSCF, PSIO_OPEN_OLD);

    dpdfile2 Qk, dpdUact;
    dpdbuf4 G, Gk, TPDM;

    // Write out the incoming TPDM
    double** TPDMmatp = TPDMmat->pointer();
    global_dpd_->buf4_init(&TPDM, PSIF_MCSCF, 0,
                ints_->DPD_ID("[X,X]"), ints_->DPD_ID("[X,X]"),
                ints_->DPD_ID("[X,X]"), ints_->DPD_ID("[X,X]"), 0, "CI Qk TPDM (XX|XX)");

    for (int h=0; h<nirrep_; h++){
        global_dpd_->buf4_mat_irrep_init(&TPDM, h);
    }

    // 4 fold symmetry
    for(int p = 0; p < nact_; p++){
      int p_sym = TPDM.params->psym[p];

      // for(int q = 0; q <= p; q++){
      for(int q = 0; q < nact_; q++){
        int q_sym = TPDM.params->psym[q];
        int pq_sym = p_sym^q_sym;
        int pq = TPDM.params->rowidx[p][q];

        for(int r = 0; r < nact_; r++){
          int r_sym = TPDM.params->psym[r];

          for(int s = 0; s < nact_; s++){
          // for(int s = 0; s <= r; s++){
            int s_sym = TPDM.params->psym[s];
            int rs_sym = r_sym^s_sym;
            int rs = TPDM.params->colidx[r][s];

            if (pq_sym != rs_sym) continue;

            TPDM.matrix[pq_sym][pq][rs] = TPDMmatp[p*nact_ + q][r*nact_ + s];
    }}}}

    for (int h=0; h<nirrep_; h++){
        global_dpd_->buf4_mat_irrep_wrt(&TPDM, h);
        global_dpd_->buf4_mat_irrep_close(&TPDM, h);
    }

    // We init R then a so we need 1, 0 here for index types
    global_dpd_->file2_init(&dpdUact, PSIF_MCSCF, 0, 1, 0, "Uact");

    // Copy SharedMatrix Uact into Ua
    global_dpd_->file2_mat_init(&dpdUact);

    for (int h=0; h<nirrep_; h++){
        if (!nactpi_[h]) continue;
        double** Uactp = Uact->pointer(h);
        int size = nactpi_[h] * nmopi_[h];
        C_DCOPY(size, Uactp[0], 1, dpdUact.matrix[h][0], 1);
    }

    global_dpd_->file2_mat_wrt(&dpdUact);
    global_dpd_->file2_mat_close(&dpdUact);


    // Rotate the two electron integrals
    global_dpd_->buf4_init(&Gk, PSIF_MCSCF, 0,
                ints_->DPD_ID("[X,X]"), ints_->DPD_ID("[X,R]"),
                ints_->DPD_ID("[X,X]"), ints_->DPD_ID("[X,R]"), 0, "Rotated MO Ints (XX|XR)");


    // \kappa_{uP}g_{tPvR} -> Gk_{tuvR}
    global_dpd_->buf4_init(&G, PSIF_LIBTRANS_DPD, 0,
                ints_->DPD_ID("[X,R]"), ints_->DPD_ID("[X,R]"),
                ints_->DPD_ID("[X,R]"), ints_->DPD_ID("[X,R]"), 0, "MO Ints (XR|XR)");

    global_dpd_->contract424(&G, &dpdUact, &Gk, 1, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&G);

    if (debug){
        outfile->Printf("First contraction done\n\n");
        global_dpd_->buf4_print(&Gk, "outfile", 1);
    }


    // \kappa_{tP}g_{PuvR} -> Gk_{tuvR}
    global_dpd_->buf4_copy(&Gk, PSIF_MCSCF, "Tran Copy Ints (XX|XR)");
    global_dpd_->buf4_close(&Gk);
    global_dpd_->buf4_init(&G, PSIF_MCSCF, 0,
                ints_->DPD_ID("[X,X]"), ints_->DPD_ID("[X,R]"),
                ints_->DPD_ID("[X,X]"), ints_->DPD_ID("[X,R]"), 0, "Tran Copy Ints (XX|XR)");
    global_dpd_->buf4_sort_axpy(&G, PSIF_MCSCF, qprs,
                            ints_->DPD_ID("[X,X]"), ints_->DPD_ID("[X,R]"),
                            "Rotated MO Ints (XX|XR)", 1.0);

    global_dpd_->buf4_init(&Gk, PSIF_MCSCF, 0,
                ints_->DPD_ID("[X,X]"), ints_->DPD_ID("[X,R]"),
                ints_->DPD_ID("[X,X]"), ints_->DPD_ID("[X,R]"), 0, "Rotated MO Ints (XX|XR)");
    global_dpd_->buf4_close(&G);
    if (debug){
        outfile->Printf("Second contraction done\n\n");
        global_dpd_->buf4_print(&Gk, "outfile", 1);
    }



    // \kappa_{wo}g_{moxy}
    global_dpd_->buf4_init(&G, PSIF_LIBTRANS_DPD, 0,
                ints_->DPD_ID("[X,X]"), ints_->DPD_ID("[R,R]"),
                ints_->DPD_ID("[X>=X]+"), ints_->DPD_ID("[R>=R]+"), 0, "MO Ints (XX|RR)");

    global_dpd_->contract244(&dpdUact, &G, &Gk, 1, 2, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->file2_close(&dpdUact);
    if (debug){
        outfile->Printf("Third contraction done\n\n");
        global_dpd_->buf4_print(&Gk, "outfile", 1);
    }

    // \Gamma_{tuvw}g_{tuvR} -> Qk_{wR}
    global_dpd_->file2_init(&Qk, PSIF_MCSCF, 0, 1, 0, "Qk");
    global_dpd_->contract442(&TPDM, &Gk, &Qk, 3, 3, 1.0, 0.0);

    global_dpd_->buf4_close(&TPDM);
    global_dpd_->buf4_close(&Gk);


    SharedMatrix Qkmat = SharedMatrix(new Matrix(&Qk));
    global_dpd_->file2_close(&Qk);

    psio_->close(PSIF_LIBTRANS_DPD, 1);
    psio_->close(PSIF_MCSCF, 1);

    // Qkmat->print();
    // Transform last index
    SharedMatrix Qmat = compute_Q(TPDMmat);
    Qkmat->gemm(false, false, -1.0, Qmat, U, 1.0);
    if (debug){
        Qkmat->print();

    }
    timer_off("SOMCSCF: Qk matrix");

    return Qkmat;

}// End DiskSOMCSCF object


/// IncoreSOMCSCF class
IncoreSOMCSCF::IncoreSOMCSCF(std::shared_ptr<JK> jk,
            SharedMatrix AOTOSO, SharedMatrix H) :
            SOMCSCF(jk, AOTOSO, H)
{
    eri_tensor_set_ = false;
}
IncoreSOMCSCF::~IncoreSOMCSCF()
{
}
SharedMatrix IncoreSOMCSCF::compute_Q(SharedMatrix TPDM)
{
    if (!eri_tensor_set_){
        throw PSIEXCEPTION("IncoreSOMCSCF: Eri tensors were not set!");
    }

    timer_on("SOMCSCF: Q matrix");

    // G_mwxy TPDM_vwxy -> Q_mv
    SharedMatrix denQ(new Matrix("Dense Qvn", nact_, nmo_));
    double** denQp = denQ->pointer();

    int nact3 = nact_ * nact_ * nact_;
    double** TPDMp = TPDM->pointer();
    double** aaaRp = mo_aaar_->pointer();
    ///KPH found that this didn't work for my purpose.  Not sure if it works for anyone else.  No test cases for this.
    //C_DGEMM('N','N',nact_,nmo_,nact3,1.0,TPDMp[0],nact3,aaaRp[0],nact3,1.0,denQp[0],nmo_);
    ///TPDM_vwxy G_mwxy -> Q_vm
    C_DGEMM('N','T',nact_,nmo_,nact3,1.0,TPDMp[0],nact3,aaaRp[0],nact3,1.0,denQp[0],nmo_);

    // Symmetry block Q
    SharedMatrix Q(new Matrix("Qvn", nirrep_, nactpi_, nmopi_));

    int offset_act = 0;
    int offset_nmo = 0;
    for (int h=0; h<nirrep_; h++){
        if (!nactpi_[h] || !nmopi_[h]){
            offset_nmo += nmopi_[h];
            continue;
        }

        double* Qp = Q->pointer(h)[0];
        for (int i=0, target=0; i<nactpi_[h]; i++){
            for (int j=0; j<nmopi_[h]; j++){
                Qp[target++] = denQp[offset_act + i][offset_nmo + j];
            }
        }
        offset_act += nactpi_[h];
        offset_nmo += nmopi_[h];
    }

    timer_off("SOMCSCF: Q matrix");
    return Q;
}
void IncoreSOMCSCF::set_act_MO(void)
{
    if (eri_tensor_set_){
        matrices_["actMO"] = mo_aaaa_;
    }
    else{
        throw PSIEXCEPTION("IncoreSOMCSCF: ERI tensors were not set!");
    }
}
SharedMatrix IncoreSOMCSCF::compute_Qk(SharedMatrix TPDM, SharedMatrix U, SharedMatrix Uact)
{
    throw PSIEXCEPTION("IncoreSOMCSCF::Qk: Qk does not yet.");
}
void IncoreSOMCSCF::set_eri_tensors(SharedMatrix aaaa, SharedMatrix aaar)
{
    mo_aaaa_ = aaaa;
    mo_aaar_ = aaar;
    eri_tensor_set_ = true;
}// End DiskSOMCSCF object

} // Namespace psi
