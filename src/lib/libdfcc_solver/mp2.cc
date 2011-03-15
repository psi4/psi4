#include "mp2.h"
#include <libmints/mints.h>
#include <lib3index/3index.h>
#include <libqt/qt.h>
#include <psiconfig.h>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef HAVE_MKL
#include <mkl.h>
#endif

using namespace boost;
using namespace psi;

namespace psi { namespace dfcc {

MP2::MP2(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
  : CC(options, psio, chkpt)
{
  print_header();
  CC::print_header();

  common_init();
}
MP2::~MP2()
{
}
void MP2::common_init()
{
    mp2_algorithm_ = options_.get_str("MP2_ALGORITHM");
}
double MP2::compute_energy()
{
  energies_["MP2J Energy"] = 0.0;
  energies_["MP2K Energy"] = 0.0;

  if (mp2_algorithm_ == "DF") {
    compute_DF_MP2();
  } else if (mp2_algorithm_ == "SOS") {
    compute_OS_MP2();
  } else if (mp2_algorithm_ == "MOS") {
    compute_OS_MP2();
  } else if (mp2_algorithm_ == "PS") {
    compute_PS_MP2();
  } else if (mp2_algorithm_ == "PS2") {
    compute_OS_MP2();
    compute_PS2_MP2();
  } else if (mp2_algorithm_ == "PS3") {
    compute_PS3_MP2();
  }

  energies_["Opposite-Spin Energy"] = 0.5*energies_["MP2J Energy"];
  energies_["Same-Spin Energy"] = 0.5*energies_["MP2J Energy"] +  energies_["MP2K Energy"];
  energies_["Correlation Energy"] = energies_["MP2J Energy"] + energies_["MP2K Energy"];
  energies_["Total Energy"] = energies_["Reference Energy"] + energies_["Correlation Energy"];

  energies_["SCS Opposite-Spin Energy"] = 0.5*oss_*energies_["MP2J Energy"];
  energies_["SCS Same-Spin Energy"] = 0.5*sss_*energies_["MP2J Energy"] +  sss_*energies_["MP2K Energy"];
  energies_["SCS Correlation Energy"] = energies_["SCS Opposite-Spin Energy"] + energies_["SCS Same-Spin Energy"];
  energies_["SCS Total Energy"] = energies_["Reference Energy"] + energies_["SCS Correlation Energy"];

  print_energies();
  return energies_["Total Energy"];
}
void MP2::compute_DF_MP2()
{
    shared_ptr<DFTensor> df(new DFTensor(psio_, basisset_, ribasis_));
    df->form_OV_integrals((ULI)(0.9*(double)doubles_), C_aocc_, C_avir_, false, fitting_algorithm_, fitting_condition_, schwarz_cutoff_);
    int nocc = df->nocc();
    int nvir = df->nvir();
    int naux = df->naux();

    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif
    int rank = 0;

    #ifdef HAVE_MKL
        int mkl_nthread = mkl_get_max_threads();
        mkl_set_num_threads(1);
    #endif

    double*** Iab = new double**[nthread];
    for (int thread = 0; thread < nthread; thread++)
        Iab[thread] = block_matrix(nvir, nvir);

    shared_ptr<TensorChunk> chunk1;
    shared_ptr<TensorChunk> chunk2;

    unsigned long int scratch = nthread*nvir*(ULI)nvir;
    if (nocc*nvir*(ULI)naux < doubles_ - scratch) {
        // 1 block
        chunk1 = df->get_ov_iterator(doubles_ - scratch);
        chunk2 = chunk1;
    } else {
        // 2 blocks
        chunk1 = df->get_ov_iterator((doubles_ - scratch) / 2L);
        chunk2 = df->get_ov_iterator((doubles_ - scratch) / 2L);
    }
    int nblock = chunk1->nblock();
    shared_ptr<Matrix> Qia1 = chunk1->chunk();
    shared_ptr<Matrix> Qia2 = chunk2->chunk();
    double** Qia1p = Qia1->pointer();
    double** Qia2p = Qia2->pointer();

    // ==> Energy Evaluation <== //
    double E_MP2J = 0.0;
    double E_MP2K = 0.0;

    // Outer loop is blocks
    for (int block1 = 0; block1 < nblock; block1++) {
        chunk1->read_block(block1);
        Qia1->print();
        int istart = chunk1->block_start() / nvir;
        int isize  = chunk1->block_size()  / nvir;

        // First do B = B blocks (diagonal)
        // TODO unroll and thread
        for (int ilocal = 0; ilocal < isize; ilocal++) {
            rank = 0;

            for (int jlocal = 0; jlocal <= ilocal; jlocal++) {
                int i = ilocal + istart;
                int j = jlocal + istart;

                C_DGEMM('N','T', nvir, nvir, naux, 1.0, &Qia1p[ilocal*nvir][0], naux,
                    &Qia1p[jlocal*nvir][0], naux, 0.0, Iab[rank][0], nvir);

                for (int a = 0; a < nvir; a++) {
                    for (int b = 0; b < nvir; b++) {
                        double iajb = Iab[rank][a][b];
                        double ibja = Iab[rank][b][a];
                        double denom = 1.0 / (evals_avirp_[a] + evals_avirp_[b]
                                             -evals_aoccp_[i] - evals_aoccp_[j]);
                        double perm_scale = ((i == j) ? 1.0 : 2.0);

                        E_MP2J += -2.0 * perm_scale*iajb*iajb*denom;
                        E_MP2K +=  1.0 * perm_scale*iajb*ibja*denom;
                    }
                }
            }
        }

        // Now do B != C blocks (lower triangle)
        // TODO unroll and thread
        for (int block2 = 0; block2 < block1; block2++) {
            chunk2->read_block(block2);
            int jstart = chunk2->block_start() / nvir;
            int jsize  = chunk2->block_size()  / nvir;

            for (int ilocal = 0; ilocal < isize; ilocal++) {
                rank = 0;

                for (int jlocal = 0; jlocal < jsize; jlocal++) {
                    int i = ilocal + istart;
                    int j = jlocal + jstart;

                    C_DGEMM('N','T', nvir, nvir, naux, 1.0, &Qia1p[ilocal*nvir][0], naux,
                        &Qia2p[jlocal*nvir][0], naux, 0.0, Iab[rank][0], nvir);

                    for (int a = 0; a < nvir; a++) {
                        for (int b = 0; b < nvir; b++) {
                            double iajb = Iab[rank][a][b];
                            double ibja = Iab[rank][b][a];
                            double denom = 1.0 / (evals_avirp_[a] + evals_avirp_[b]
                                                 -evals_aoccp_[i] - evals_aoccp_[j]);

                            E_MP2J += -4.0 * iajb*iajb*denom;
                            E_MP2K +=  2.0 * iajb*ibja*denom;
                        }
                    }
                }
            }
        }
    }

    #ifdef HAVE_MKL
        mkl_set_num_threads(mkl_nthread);
    #endif

    // Free scratch
    for (int thread = 0; thread < nthread; thread++)
        free_block(Iab[thread]);
    delete[] Iab;

    energies_["MP2J Energy"] = E_MP2J;
    energies_["MP2K Energy"] = E_MP2K;
}
void MP2::compute_OS_MP2()
{
    shared_ptr<DFTensor> df(new DFTensor(psio_, basisset_, ribasis_));
    df->form_OV_integrals((ULI)(0.9*(double)doubles_), C_aocc_, C_avir_, false, fitting_algorithm_, fitting_condition_, schwarz_cutoff_);
    int nocc = df->nocc();
    int nvir = df->nvir();
    int naux = df->naux();

    shared_ptr<LaplaceDenominator> denom(new LaplaceDenominator(evals_aocc_, evals_avir_, denominator_delta_));
    shared_ptr<Matrix> Tau = denom->denominator();
    denom->debug();
    int nvector = denom->nvector();

    unsigned long int scratch = naux*(ULI)naux;
    shared_ptr<Matrix> Z(new Matrix("Z_QQ'^w", naux, naux));

    shared_ptr<TensorChunk> chunk = df->get_ov_iterator((doubles_ - scratch)/2);
    int max_rows = chunk->max_rows();
    int nblock = chunk->nblock();
    shared_ptr<Matrix> Qia = chunk->chunk();
    shared_ptr<Matrix> Qiaw(new Matrix("(Q|ia)^w", max_rows, naux));

    double** Qiap = Qia->pointer();
    double** Qiawp = Qiaw->pointer();
    double** Zp = Z->pointer();
    double** Taup = Tau->pointer();

    // ==> Energy Evaluation <== //
    double E_MP2J = 0.0;
    for (int block = 0; block < nblock; block++) {

        chunk->read_block(block);
        int start = chunk->block_start();
        int size = chunk->block_size();

        for (int w = 0; w < nvector; w++) {
            memcpy(static_cast<void*>(Qiawp[0]), static_cast<void*>(Qiap[0]), size*(ULI)naux*sizeof(double));

            for (int ia = 0; ia < size; ia++)
                C_DSCAL(naux, sqrt(Taup[w][ia + start]), Qiawp[ia], 1);

            C_DGEMM('T','N', naux, naux, size, 1.0, Qiawp[0], naux, Qiawp[0], naux, 0.0, Zp[0], naux);

            E_MP2J -= 2.0 * C_DDOT(naux*(ULI)naux, Zp[0], 1, Zp[0], 1);
        }
    }

    energies_["MP2J Energy"] = E_MP2J;
}
void MP2::compute_PS_MP2()
{
    throw FeatureNotImplemented("libdfcc_solver", "psi::dfcc::MP2::compute_PS_MP2", __FILE__, __LINE__);
}
void MP2::compute_PS2_MP2()
{
    shared_ptr<PSTensor> ps(new PSTensor(psio_, basisset_, dealias_, grid_));
    ps->form_OV_integrals((ULI)(0.9*(double)doubles_), C_aocc_, C_avir_, false, schwarz_cutoff_);

    int naux = ps->naux();
    int nocc = ps->nocc();
    int nvir = ps->nvir();

    shared_ptr<LaplaceDenominator> denom(new LaplaceDenominator(evals_aocc_, evals_avir_, denominator_delta_));
    shared_ptr<Matrix> Tau = denom->denominator();
    shared_ptr<Matrix> Tau_i = denom->denominator_occ();
    shared_ptr<Matrix> Tau_a = denom->denominator_vir();
    int nvector = denom->nvector();

    double** Taup = Tau->pointer();
    double** Tau_ip = Tau_i->pointer();
    double** Tau_ap = Tau_a->pointer();

    ULI scratch = 2L * (naux*nocc + naux*nvir) + 2L * (naux*naux);
    shared_ptr<TensorChunk> chunk = ps->get_ov_iterator((doubles_ - scratch)/2);
    int max_rows = chunk->max_rows();
    int nblock = chunk->nblock();

    shared_ptr<Matrix> Q = ps->form_Q(C_aocc_);
    shared_ptr<Matrix> X = ps->form_X(C_avir_);
    shared_ptr<Matrix> A = chunk->chunk();

    // Copies of Q and X for Tau scaling
    shared_ptr<Matrix> Qw(new Matrix("Q_i^Pw (Scaled)", nocc, naux));
    shared_ptr<Matrix> Xw(new Matrix("X_a^Pw (Scaled)", nvir, naux));

    // Copy of A_jb^P for Tau scaling
    shared_ptr<Matrix> Aw(new Matrix("A_jb^Pw (Scaled)", max_rows, naux));

    double** Qp = Q->pointer();
    double** Qwp = Qw->pointer();
    double** Xp = X->pointer();
    double** Xwp = Xw->pointer();
    double** Ap = A->pointer();
    double** Awp = Aw->pointer();

    // The Z, Lambda, and Gamma intermediates
    shared_ptr<Matrix> Z(new Matrix("Z_j^PQw", naux, naux));
    shared_ptr<Matrix> L(new Matrix("L^PQw", naux, naux));

    double** Zp = Z->pointer();
    double** Lp = L->pointer();

    double E_MP2K = 0.0;

    for (int block = 0; block < nblock; block++) {

        chunk->read_block(block);
        int start = chunk->block_start();
        int size = chunk->block_size();
        int istart = start / nvir;
        int isize = size / nvir;

        for (int w = 0; w < nvector; w++) {

            memcpy(static_cast<void*>(Awp[0]), static_cast<void*>(Ap[0]), size*(ULI)naux*sizeof(double));
            for (int ia = 0; ia < size; ia++)
                C_DSCAL(naux, sqrt(Taup[w][ia + start]), Awp[ia], 1);

            memcpy(static_cast<void*>(Xwp[0]), static_cast<void*>(Xp[0]), nvir*(ULI)naux*sizeof(double));
            for (int a = 0; a < nvir; a++)
                C_DSCAL(naux, sqrt(Tau_ap[w][a]), Xwp[a], 1);

            memcpy(static_cast<void*>(Qwp[0]), static_cast<void*>(Qp[0]), nocc*(ULI)naux*sizeof(double));
            for (int i = 0; i < nocc; i++)
                C_DSCAL(naux, sqrt(Tau_ip[w][i]), Qwp[i], 1);

            C_DGEMM('T','N',naux,naux,nocc,1.0,Qwp[0],naux,Qwp[0],naux,0.0,Lp[0],naux);
            L->print();

            for (int j = 0; j < isize; j++) {
                C_DGEMM('T','N',naux,naux,nvir,1.0,Xwp[0],naux,Awp[j*nvir],naux,0.0,Zp[0],naux);
                Z->print();

                for (ULI QP = 0; QP < naux*(ULI)naux; QP++)
                    Zp[0][QP] *= Zp[0][QP];
                Z->print();

                E_MP2K += C_DDOT(naux*(ULI)naux, Zp[0], 1, Lp[0], 1);
            }
        }
    }

    energies_["MP2K Energy"] = E_MP2K;
}
void MP2::compute_PS3_MP2()
{
    shared_ptr<PSTensor> ps(new PSTensor(psio_, basisset_, dealias_, grid_));
    ps->form_OV_integrals((ULI)(0.9*(double)doubles_), C_aocc_, C_avir_, false, schwarz_cutoff_);

    int naux = ps->naux();
    int nocc = ps->nocc();
    int nvir = ps->nvir();

    shared_ptr<LaplaceDenominator> denom(new LaplaceDenominator(evals_aocc_, evals_avir_, denominator_delta_));
    shared_ptr<Matrix> Tau = denom->denominator();
    shared_ptr<Matrix> Tau_i = denom->denominator_occ();
    shared_ptr<Matrix> Tau_a = denom->denominator_vir();
    int nvector = denom->nvector();

    double** Taup = Tau->pointer();
    double** Tau_ip = Tau_i->pointer();
    double** Tau_ap = Tau_a->pointer();

    ULI scratch = 2L * (naux*nocc + naux*nvir) + 3L * (naux*naux);
    shared_ptr<TensorChunk> chunk = ps->get_ov_iterator((doubles_ - scratch)/2);
    int max_rows = chunk->max_rows();
    int nblock = chunk->nblock();

    shared_ptr<Matrix> Q = ps->form_Q(C_aocc_);
    shared_ptr<Matrix> X = ps->form_X(C_avir_);
    shared_ptr<Matrix> A = chunk->chunk();

    // Copies of Q and X for Tau scaling
    shared_ptr<Matrix> Qw(new Matrix("Q_i^Pw (Scaled)", nocc, naux));
    shared_ptr<Matrix> Xw(new Matrix("X_a^Pw (Scaled)", nvir, naux));

    // Copy of A_jb^P for Tau scaling
    shared_ptr<Matrix> Aw(new Matrix("A_jb^Pw (Scaled)", max_rows, naux));

    double** Qp = Q->pointer();
    double** Qwp = Qw->pointer();
    double** Xp = X->pointer();
    double** Xwp = Xw->pointer();
    double** Ap = A->pointer();
    double** Awp = Aw->pointer();

    // The Z, Lambda, and Gamma intermediates
    shared_ptr<Matrix> Z(new Matrix("Z_j^PQw", naux, naux));
    shared_ptr<Matrix> L(new Matrix("L^PQw", naux, naux));
    shared_ptr<Matrix> G(new Matrix("G^PQw", naux, naux));

    double** Zp = Z->pointer();
    double** Lp = L->pointer();
    double** Gp = G->pointer();

    double E_MP2J = 0.0;
    double E_MP2K = 0.0;

    for (int block = 0; block < nblock; block++) {

        chunk->read_block(block);
        int start = chunk->block_start();
        int size = chunk->block_size();
        int istart = start / nvir;
        int isize = size / nvir;

        for (int w = 0; w < nvector; w++) {

            memcpy(static_cast<void*>(Awp[0]), static_cast<void*>(Ap[0]), size*(ULI)naux*sizeof(double));
            for (int ia = 0; ia < size; ia++)
                C_DSCAL(naux, sqrt(Taup[w][ia + start]), Awp[ia], 1);

            memcpy(static_cast<void*>(Xwp[0]), static_cast<void*>(Xp[0]), nvir*(ULI)naux*sizeof(double));
            for (int a = 0; a < nvir; a++)
                C_DSCAL(naux, sqrt(Tau_ap[w][a]), Xwp[a], 1);

            memcpy(static_cast<void*>(Qwp[0]), static_cast<void*>(Qp[0]), nocc*(ULI)naux*sizeof(double));
            for (int i = 0; i < nocc; i++)
                C_DSCAL(naux, sqrt(Tau_ip[w][i]), Qwp[i], 1);

            // MP2J
            C_DGEMM('T','N',naux,naux,nocc,1.0,Qwp[0],naux,Qwp[0],naux,0.0,Lp[0],naux);
            C_DGEMM('T','N',naux,naux,nvir,1.0,Xwp[0],naux,Xwp[0],naux,0.0,Gp[0],naux);

            for (ULI PQ = 0L; PQ < naux*(ULI)naux; PQ++)
                Gp[0][PQ] *= Lp[0][PQ];

            C_DGEMM('T','N', naux, naux, size, 1.0, Awp[0], naux, Awp[0], naux, 0.0, Zp[0], naux);

            E_MP2J -= 2.0 * C_DDOT(naux*(ULI)naux, Zp[0], 1, Gp[0], 1);

            // MP2K
            for (int j = 0; j < isize; j++) {
                C_DGEMM('T','N',naux,naux,nvir,1.0,Xwp[0],naux,Awp[j*nvir],naux,0.0,Zp[0],naux);

                for (ULI QP = 0; QP < naux*(ULI)naux; QP++)
                    Zp[0][QP] *= Zp[0][QP];

                E_MP2K += C_DDOT(naux*(ULI)naux, Zp[0], 1, Lp[0], 1);
            }
        }
    }

    energies_["MP2J Energy"] = E_MP2J;
    energies_["MP2K Energy"] = E_MP2K;
}

void MP2::print_header()
{
    fprintf(outfile, "\t ********************************************************\n");
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t *                        DF-MP2                        *\n");
    fprintf(outfile, "\t *    2nd-Order Density-Fitted Moller-Plesset Theory    *\n");
    fprintf(outfile, "\t *        with Laplace and Pseudospectral Grids         *\n");
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t *            Rob Parrish and Ed Hohenstein             *\n");
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t ********************************************************\n");
    fprintf(outfile, "\n");
}
void MP2::print_energies()
{
    fprintf(outfile, "\n");
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t ====================> MP2 Energies <==================== \n");
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Reference Energy",         energies_["Reference Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "MP2J Energy",              energies_["MP2J Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "MP2K Energy",              energies_["MP2K Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Same-Spin Energy",         energies_["Same-Spin Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Opposite-Spin Energy",     energies_["Opposite-Spin Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Correlation Energy",       energies_["Correlation Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Total Energy",             energies_["Total Energy"]);
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t ==================> SCS-MP2 Energies <================== \n");
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t %-25s = %24.16f [-]\n", "SCS Same-Spin Scale",      sss_);
    fprintf(outfile, "\t %-25s = %24.16f [-]\n", "SCS Opposite-Spin Scale",  oss_);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "SCS Same-Spin Energy",     energies_["SCS Same-Spin Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "SCS Opposite-Spin Energy", energies_["SCS Opposite-Spin Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "SCS Correlation Energy",   energies_["SCS Correlation Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "SCS Total Energy",         energies_["SCS Total Energy"]);
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\n");
    fflush(outfile);
}

}}
