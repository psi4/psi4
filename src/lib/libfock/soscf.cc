/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include <libmints/mints.h>
#include <libqt/qt.h>
#include <psi4-dec.h>
#include <libthce/thce.h>
#include <libthce/lreri.h>


#include "soscf.h"
#include "jk.h"

namespace psi {

/// SOMCSCF class

SOMCSCF::SOMCSCF(boost::shared_ptr<JK> jk, boost::shared_ptr<DFERI> df, SharedMatrix AOTOSO,
                SharedMatrix H)
{
    jk_ = jk;
    dferi_ = df;
    matrices_["H"] = H;
    matrices_["AOTOSO"] = AOTOSO;
    nao_ = AOTOSO->rowspi()[0];
    casscf_ = true;
}
SOMCSCF::~SOMCSCF()
{
}
void SOMCSCF::set_ras(std::vector<Dimension> ras_spaces)
{
    ras_spaces_ = ras_spaces;
    casscf_ = false;
    // Check the sum of ras equals the act size
}
SharedMatrix SOMCSCF::Ck(SharedMatrix x)
{

    SharedMatrix tmp(new Matrix("Ck", nirrep_, nsopi_, nmopi_));

    // Form full antisymmetric matrix
    for (size_t h=0; h<nirrep_; h++){

        if (!noapi_[h] || !navpi_[h]) continue;
        double** tp = tmp->pointer(h);
        double*  xp = x->pointer(h)[0];

        // Matrix::schmidt orthogonalizes rows not columns so we need to transpose
        for (size_t i=0, target=0; i<noapi_[h]; i++){
            for (size_t a=noccpi_[h]; a<nmopi_[h]; a++){
                // outfile->Printf("%d %d %d | %lf\n", h, i, a, xp[target]);
                tp[i][a] = xp[target];
                tp[a][i] = -1.0 * xp[target++];
            }
        }
    }

    // Build exp(U) = 1 + U + 0.5 U U
    SharedMatrix U = tmp->clone();
    for (size_t h=0; h<nirrep_; h++){
        double** Up = U->pointer(h);
        for (size_t i=0; i<U->rowspi()[h]; i++){
            Up[i][i] += 1.0;
        }
    }
    U->gemm(false, false, 0.5, tmp, tmp, 1.0);

    // We did not fully exponentiate the matrix, need to orthogonalize
    U->schmidt();

    // C' = C U
    tmp->gemm(false, false, 1.0, matrices_["C"], U, 0.0);

    return tmp;
}
void SOMCSCF::update(SharedMatrix Cocc, SharedMatrix Cact, SharedMatrix Cvir,
            SharedMatrix OPDM, SharedMatrix TPDM)
{
    // => Update orbitals and density matrices <= //
    std::vector<boost::shared_ptr<Matrix> > fullC;
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

    // => AO C matrices <= //
    // Only need for DF SOMCSCF
    matrices_["AO_Cact"] = SharedMatrix(new Matrix("AO_Cact", nao_, nact_));
    matrices_["AO_Cocc"] = SharedMatrix(new Matrix("AO_Cocc", nao_, nocc_));
    matrices_["AO_Cvir"] = SharedMatrix(new Matrix("AO_Cvir", nao_, nvir_));

    /// Not a lot of checks in case someone did something stupid, need to think about that
    int offset_act = 0;
    int offset_occ = 0;
    int offset_vir = 0;
    for (int h = 0; h < nirrep_; ++h) {
        int hnso = matrices_["AOTOSO"]->colspi()[h];
        if (hnso == 0) continue;
        double** Up = matrices_["AOTOSO"]->pointer(h);

        // occupied
        if (noccpi_[h]){
            double** CAOp = matrices_["AO_Cocc"]->pointer();
            double** CSOp = matrices_["Cocc"]->pointer(h);
            C_DGEMM('N','N',nao_,noccpi_[h],hnso,1.0,Up[0],hnso,CSOp[0],noccpi_[h],0.0,&CAOp[0][offset_occ],nocc_);
            offset_occ += noccpi_[h];
        }
        // active
        if (nactpi_[h]){
            double** CAOp = matrices_["AO_Cact"]->pointer();
            double** CSOp = matrices_["Cact"]->pointer(h);
            C_DGEMM('N','N',nao_,nactpi_[h],hnso,1.0,Up[0],hnso,CSOp[0],nactpi_[h],0.0,&CAOp[0][offset_act],nact_);
            offset_act += nactpi_[h];
        }
        // virtual
        if (nvirpi_[h]){
            double** CAOp = matrices_["AO_Cvir"]->pointer();
            double** CSOp = matrices_["Cvir"]->pointer(h);
            C_DGEMM('N','N',nao_,nvirpi_[h],hnso,1.0,Up[0],hnso,CSOp[0],nvirpi_[h],0.0,&CAOp[0][offset_vir],nvir_);
            offset_vir += nvirpi_[h];
        }
    }
    std::vector<boost::shared_ptr<Matrix> > full_AOC;
    full_AOC.push_back(matrices_["AO_Cocc"]);
    full_AOC.push_back(matrices_["AO_Cact"]);
    full_AOC.push_back(matrices_["AO_Cvir"]);

    matrices_["AO_C"] = Matrix::horzcat(full_AOC);
    matrices_["AO_C"]->set_name("AO_C");
    // matrices_["C"]->print();
    // matrices_["AO_C"]->print();


    // => Build generalized inactive and active Fock matrices <= //
    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();

    // For inactive Fock
    Cl.push_back(matrices_["Cocc"]);
    Cr.push_back(matrices_["Cocc"]);

    // For active Fock
    SharedMatrix CL_COPDM = Matrix::doublet(matrices_["Cact"], matrices_["OPDM"]);
    Cl.push_back(CL_COPDM);
    Cr.push_back(matrices_["Cact"]);

    jk_->compute();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    J[0]->scale(2.0);
    J[0]->subtract(K[0]);
    J[0]->add(matrices_["H"]);
    K[0]->gemm(false, false, 1.0, J[0], matrices_["C"], 0.0);
    matrices_["IFock"] = Matrix::doublet(matrices_["C"], K[0], true, false);
    matrices_["IFock"]->set_name("IFock");
    // matrices_["IFock"]->print();

    K[1]->scale(0.5);
    J[1]->subtract(K[1]);
    K[1]->gemm(false, false, 1.0, J[1], matrices_["C"], 0.0);
    matrices_["AFock"] = Matrix::doublet(matrices_["C"], K[1], true, false);
    matrices_["AFock"]->set_name("AFock");
    // matrices_["AFock"]->print();
    // matrices_["OPDM"]->print();
    // matrices_["C"]->print();
    // matrices_["Cact"]->print();

    // => Q matrix <= //
    // This will be moved out to the DFSOMCSCF and PKSOMCSCF objects
    // For now we will implement DFSOMCSCF
    // This is done simply
    // d_vwxy I_nwxy => Q_vn
    dferi_->clear();
    dferi_->set_C(matrices_["AO_C"]);
    dferi_->add_space("o", 0, nocc_);
    dferi_->add_space("a", nocc_, nocc_ + nact_);
    dferi_->add_space("v", nocc_ + nact_, nocc_ + nact_ + nvir_);
    dferi_->add_space("N", 0, nmo_);

    // Is it smart enough to order then untranspose?
    // In the future build this once then slice it
    dferi_->add_pair_space("aaQ", "a", "a");
    dferi_->add_pair_space("NaQ", "N", "a");
    dferi_->add_pair_space("NNQ", "N", "N");

    // dferi_->print_header(2);
    dferi_->compute();
    std::map<std::string, boost::shared_ptr<Tensor> >& dfints = dferi_->ints();

    int nQ = dferi_->size_Q();
    int nact2 = nact_ * nact_;
    double* TPDMp = TPDM->pointer()[0];

    // Load aaQ
    boost::shared_ptr<Tensor> aaQT = dfints["aaQ"];
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
    boost::shared_ptr<Tensor> NaQT = dfints["NaQ"];
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
    matrices_["Q"] = SharedMatrix(new Matrix("Qvn", nirrep_, nactpi_, nmopi_));

    offset_act = 0;
    int offset_nmo = 0;
    for (int h=0; h<nirrep_; h++){
        if (!nactpi_[h] || !nmopi_[h]) continue;

        double* Qp = matrices_["Q"]->pointer(h)[0];
        for (int i=0, target=0; i<nactpi_[h]; i++){
            for (int j=0; j<nmopi_[h]; j++){
                Qp[target++] = denQp[offset_act + i][offset_nmo + j];
            }
        }
        offset_act += nactpi_[h];
        offset_nmo += nmopi_[h];
    }
    // matrices_["Q"]->print();
    // outfile->Printf("Finished building Q\n");


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
    zero_redundant(matrices_["Gradient"]);
    matrices_["Precon"] = H_approx_diag();
}
SharedMatrix SOMCSCF::H_approx_diag()
{
    // Build (aa|aa)
    std::map<std::string, boost::shared_ptr<Tensor> >& dfints = dferi_->ints();
    int nQ = dferi_->size_Q();
    boost::shared_ptr<Tensor> aaQT = dfints["aaQ"];
    SharedMatrix aaQ(new Matrix("aaQ", nact_ * nact_, nQ));
    double* aaQp = aaQ->pointer()[0];
    FILE* aaQF = aaQT->file_pointer();
    fseek(aaQF,0L,SEEK_SET);
    fread(aaQp, sizeof(double), nact_ * nact_ * nQ, aaQF);
    SharedMatrix actMO = Matrix::doublet(aaQ, aaQ, false, true);
    aaQ.reset();

    // d_tuvw I_tuvw -> dI_t, D_tu IF_tu -> DIF_t
    double** actMOp = actMO->pointer();
    double** TPDMp = matrices_["TPDM"]->pointer();
    int relact = 0;
    int nact3 = nact_*nact_*nact_;

    SharedVector dI(new Vector("dI", nirrep_, nactpi_));
    SharedVector DIF(new Vector("IF * OPDM", nirrep_, nactpi_));
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

            // Need to figure out the right rotations
            std::vector<std::pair<int, int>> iter;

            int offset_col = 0;
            int offset_row = 0;
            // Loop over spaces last space will have no rotations
            for (int nras = 0; nras < ras_spaces_.size()-1; nras++){
                int ras_size = ras_spaces_[nras][h];
                offset_col += ras_size;

                // Loop over pairs
                for (int i=offset_row; i<(offset_row + ras_size); i++){
                for (int a=offset_col; a<nactpi_[h]; a++){
                    // outfile->Printf("Running %d %d %d\n", nras, i, a);
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
                    int end = offset_act + nactpi_[h];
                    for (int u=offset_act; u<end; u++){
                        int pu = p * nact_ + u;
                        int qu = q * nact_ + u;
                        for (int v=offset_act; v<end; v++){
                            int pv = p*nact_ + v;
                            int qv = q*nact_ + v;
                            int uv = u*nact_ + v;
                            value += 2.0 * TPDMp[pu][pv] * actMOp[qu][qv];
                            value += 2.0 * TPDMp[qu][qv] * actMOp[pu][pv];
                            value += TPDMp[pp][uv] * actMOp[qq][uv];
                            value += TPDMp[qq][uv] * actMOp[pp][uv];
                            value -= 4.0 * TPDMp[pv][qu] * actMOp[pu][qv];
                            value -= 2.0 * TPDMp[pq][uv] * actMOp[pq][uv];
                        }
                    }
                    Hp[oi][a] += 2.0 * value;
                    // outfile->Printf("double loop thing %lf\n", value);
                    // outfile->Printf("%d %d | %lf\n", oi, a, Hp[oi][a]);

                    // Hp[oa][a] += 2.0 * value;

                } // End i loop
                } // End a loop
                offset_row += ras_size;
            } // End ras loop
            offset_act += nactpi_[h];
        } // End aa RASSCF block
    }
    return H;
}
SharedMatrix SOMCSCF::Hk(SharedMatrix x)
{

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
    // U->print();
    // Uact->print();

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
    K[0]->gemm(false, false, 1.0, J[0], matrices_["C"], 0.0);
    IFk->gemm(true, false, 1.0, matrices_["C"], K[0], 1.0);
    // outfile->Printf("The rotated inactive Fock matrix\n");
    // IFk->print();

    // Rotated active fock
    SharedMatrix ret = Matrix::doublet(matrices_["AFock"], U, false, true);
    ret->gemm(false, false, 1.0, U, matrices_["AFock"], 1.0);

    J[1]->scale(2.0);
    K[1]->scale(0.5);
    J[1]->subtract(K[1]);
    J[1]->subtract(K[1]->transpose());
    K[1]->gemm(false, false, 1.0, J[1], matrices_["C"], 0.0);
    ret->gemm(true, false, 1.0, matrices_["C"], K[1], 1.0);
    // outfile->Printf("The rotated active Fock matrix\n");
    // ret->print();

    ret->add(IFk);
    ret->scale(2.0);

    /// Build Qk

    // Remove U symmetry
    SharedMatrix dUact;
    if (nirrep_ == 1){
        dUact = Uact;
    }
    else{
        dUact = SharedMatrix(new Matrix("Dense Uact", nact_, nmo_));
        double** dUactp = dUact->pointer();
        int offset_act = 0;
        int offset_nmo = 0;
        for (int h=0; h<nirrep_; h++){
            if (!nactpi_[h]) continue;
            double** Uactp = Uact->pointer(h);
            for (int a=0; a<nactpi_[h]; a++) {
                C_DCOPY(nmopi_[h], Uactp[a], 1, dUactp[offset_act+a]+offset_nmo, 1);
            }
            offset_act += nactpi_[h];
            offset_nmo += nmopi_[h];
        }
    }
    // dUact->print();
    // Uact->print();

    int nQ = dferi_->size_Q();
    int nact2 = nact_*nact_;
    int nact3 = nact2*nact_;
    double* TPDMp = matrices_["TPDM"]->pointer()[0];
    std::map<std::string, boost::shared_ptr<Tensor> >& dfints = dferi_->ints();

    // Read NaQ
    boost::shared_ptr<Tensor> NaQT = dfints["NaQ"];
    SharedMatrix NaQ(new Matrix("NaQ", nmo_ * nact_, nQ));
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
    boost::shared_ptr<Tensor> NNQT = dfints["NNQ"];

    SharedMatrix wnQ(new Matrix("nwQ", nact_*nmo_, nQ));
    double** wnQp = wnQ->pointer();
    SharedMatrix NNQ(new Matrix("NNQ", chunk_size * nmo_, nQ));
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
    boost::shared_ptr<Tensor> aaQT = dfints["aaQ"];
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

    SharedMatrix Qk = Matrix::doublet(matrices_["Q"], U, false, true);

    int offset_act = 0;
    int offset_nmo = 0;
    for (int h=0; h<nirrep_; h++){
        if (!nactpi_[h]) continue;
        double** Qkp = Qk->pointer(h);
        for (int a=0; a<nactpi_[h]; a++){
            C_DAXPY(nmopi_[h], 1.0, dQkp[offset_act+a]+offset_nmo, 1, Qkp[a], 1);
        }
        offset_act += nactpi_[h];
        offset_nmo += nmopi_[h];
    }
    // outfile->Printf("Final Qk\n");
    // Qk->print();

    // Add in Q and zero out virtual
    for (int h=0; h<nirrep_; h++){
        if (nactpi_[h]) {
            double** Fkp = ret->pointer(h);
            double** Qkp = Qk->pointer(h);

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
    // ret->print();


    // => Orbtial Gradient <= //
    SharedMatrix hessx(new Matrix("Hessian x", nirrep_, noapi_, navpi_));

    for (int h=0; h<nirrep_; h++){
        if (!noapi_[h] || !navpi_[h]) continue;

        double** Hxp = hessx->pointer(h);
        double** Fkp = ret->pointer(h);

        for (int i=0; i<noapi_[h]; i++){
            for (int j=0; j<navpi_[h]; j++){
                int nj = noccpi_[h] + j;
                Hxp[i][j] = 2.0 * (Fkp[i][nj] - Fkp[nj][i]);

                // Ensure the hessian is zero for the active-active diagonal block
                if (nj == i){
                    Hxp[i][j] = 0.0;
                }
            }
        }
    }

    zero_redundant(hessx);
    // hessx->scale(4.0);
    // hessx->print();

    return hessx;

}
SharedMatrix SOMCSCF::approx_solve()
{
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
    // x->print();

    // Calc hessian vector product, find residual and conditioned residual
    SharedMatrix r = matrices_["Gradient"]->clone();

    SharedMatrix Ap = Hk(x);
    r->subtract(Ap);

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

        p->scale(beta);
        p->add(z);
    }// End iterations

    if (print){
        outfile->Printf("    %-4s %11.3E %10s\n", "Best", best_conv, "--");
        outfile->Printf("    ---------------------------------------\n");
        outfile->Printf("\n");
    }
    return best;
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

} // Namespace psi
