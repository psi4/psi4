/*
 *  SAPT_DFT.CC
 *
 */

#ifdef HAVE_MKL
#include <mkl.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>
#include <time.h>

#include <psifiles.h>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <libutil/quad.h>
#include <psifiles.h>

#include <libmints/mints.h>
#include <libscf_solver/integrator.h>

#include "sapt_dft.h"

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace sapt {

SAPT_DFT::SAPT_DFT(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
    : SAPT0(options, psio, chkpt)
{
    E_UCHF_disp_ = 0.0;
    E_TDDFT_disp_ = 0.0;
    E_MP2C_delta_ = E_TDDFT_disp_ - E_UCHF_disp_;
    E_MP2C_int_ = E_MP2_int_ + E_MP2C_delta_;
    get_quad_basis();
}

SAPT_DFT::~SAPT_DFT()
{
}

double SAPT_DFT::compute_energy()
{
    // THe obligatory
    print_header();
    // Allocate space
    allocate_arrays();
    // DF integrals for use with X
    fprintf(outfile, "  Computing density-fitting MO integrals.\n");
    compute_integrals();

    // This needs to go away
    compute_amplitudes();
    elst10();
    exch10();
    exch_disp20();
    disp20();
    cphf_induction();
    ind20();
    exch_ind20();
    //exit(0);
    // Yeah, all of it, sorry Ed

    // The density matrices (AO)
    // The J matrix
    fprintf(outfile, "  Computing auxiliary fitting matrix.\n");
    compute_J(); // Confirmed
    // The bloody W matrix of Dyson's Equation
    fprintf(outfile, "  Computing interelectronic interaction matrix.\n");
    compute_W();


    // Get a quadrature object for imaginary frequencies
    shared_ptr<ChebyshevIIQuadrature> quad = shared_ptr<ChebyshevIIQuadrature> ( \
        new ChebyshevIIQuadrature(options_.get_int("N_OMEGA")));

    // Perform quadrature over imaginary frequencies
    fprintf(outfile, "\n");
    fprintf(outfile, "  ---------------------------------------------------------------------------------------\n");
    fprintf(outfile, "   ===========================> CASIMIR-POLDER INTEGRATION <============================\n");
    fprintf(outfile, "  ---------------------------------------------------------------------------------------\n");
    fprintf(outfile, "   Point      Omega            Weight             E_UCHF [mH]            E_TDDFT [mH]\n");
    fprintf(outfile, "  ---------------------------------------------------------------------------------------\n");
    fflush(outfile);
    int n = 0;
    for (quad->reset(); !quad->isDone(); quad->nextPoint() ) {

        double omega = quad->getPoint();
        double weight = quad->getWeight();

        // Uncoupled response functions
        compute_X_0(omega);
        // UCHF Dispersion
        double UCHF = compute_UCHF_disp();
        // Coupled response functions
        compute_X_coup(omega);
        // TDDFT Dispersion
        double TDDFT = compute_TDDFT_disp();

        n++;
        fprintf(outfile, "   %3d   %12.8E   %12.8E   %20.12f   %20.12f\n", \
           n, omega, weight, UCHF*1000.0,TDDFT*1000.0);
        E_UCHF_disp_ += weight*UCHF;
        E_TDDFT_disp_ += weight*TDDFT;
    }
    fprintf(outfile, "  ---------------------------------------------------------------------------------------\n");
    fprintf(outfile, "    @ UCHF Dispersion Energy:  %18.12f [mH] %18.12f [kcal mol^-1]\n", E_UCHF_disp_*1000.0, \
        E_UCHF_disp_*627.509);
    fprintf(outfile, "    @ TDDFT Dispersion Energy: %18.12f [mH] %18.12f [kcal mol^-1]\n", E_TDDFT_disp_*1000.0, \
        E_TDDFT_disp_*627.509);
    fprintf(outfile, "  ---------------------------------------------------------------------------------------\n");
    fflush(outfile);


    // Free arrays
    free_arrays();

    // More fun
    print_results();

    return E_MP2C_int_;
}

void SAPT_DFT::print_header()
{
     fprintf(outfile,"                 SAPT DFT  \n");
     fprintf(outfile,"       Rob Parrish and Ed Hohenstein\n") ;
     fprintf(outfile,"             9 November 2009\n") ;
     fprintf(outfile,"\n");
     fprintf(outfile,"          Orbital Information\n");
     fprintf(outfile,"        -----------------------\n");
     fprintf(outfile,"          NSO     = %9d\n",calc_info_.nso);
     fprintf(outfile,"          NMO     = %9d\n",calc_info_.nmo);
     fprintf(outfile,"          NRI     = %9d\n",calc_info_.nri);
     fprintf(outfile,"          NOCC_A  = %9d\n",calc_info_.noccA);
     fprintf(outfile,"          NOCC_B  = %9d\n",calc_info_.noccB);
     fprintf(outfile,"          NVIR_A  = %9d\n",calc_info_.nvirA);
     fprintf(outfile,"          NVIR_B  = %9d\n\n",calc_info_.nvirB);

     #ifdef _OPENMP
     fprintf(outfile,"  Running SAPT_DFT with %d OMP threads\n\n",
       omp_get_max_threads());
     #endif

     /**
     int m = 6;
     int n = 4;

     double alpha = 1.0;
     double beta = 0.0;

     double** A = block_matrix(m,m);
     double** Atemp = block_matrix(m,m);
     double** B1 = block_matrix(m,n);
     double** B2 = block_matrix(n,m);
     double** C1 = block_matrix(m,n);
     double** C2 = block_matrix(n,m);

     for (int a = 0; a < m; a++)
        for (int b = 0; b < n; b++) {
            B1[a][b] = random() / (double)RAND_MAX;
            B2[b][a] = random() / (double)RAND_MAX;
        }

     for (int a = 0; a < m; a++)
        for (int b = 0; b < m; b++) {
            A[a][b] = random() / (double)RAND_MAX;
        }

     for (int a = 0; a < m; a++)
        for (int b = 0; b < a; b++) {
            A[a][b] = (A[a][b] + A[b][a]) / 2.0;
            A[b][a] = A[a][b];
        }

     C_DCOPY(m*m,A[0],1,Atemp[0],1);

     fprintf(outfile," A: \n");
     print_mat(A,m,m,outfile);
     fprintf(outfile," B1: \n");
     print_mat(B1,m,n,outfile);
     fprintf(outfile," B2: \n");
     print_mat(B2,n,m,outfile);

     C_DSYMM('L','U',m,n,alpha, A[0], m, B1[0], n, beta, C1[0], n);

     fprintf(outfile," A*B: \n");
     print_mat(C1,m,n,outfile);

     C_DSYMM('R','U',n,m,alpha, A[0], m, B2[0], m, beta, C2[0], m);

     fprintf(outfile," B*A: \n");
     print_mat(C2,n,m,outfile);

     C_DCOPY(m*m,Atemp[0],1,A[0],1);
     for (int a = 0; a < m; a++)
        for (int b = 0; b < a; b++)
            A[b][a] = 0.0;

     fprintf(outfile, " AL:\n");
     print_mat(A,m,m,outfile);

     C_DSYMM('L','L',m,n,alpha, A[0], m, B1[0], n, beta, C1[0], n);

     fprintf(outfile," A*B: \n");
     print_mat(C1,m,n,outfile);

     C_DSYMM('R','L',n,m,alpha, A[0], m, B2[0], m, beta, C2[0], m);

     fprintf(outfile," B*A: \n");
     print_mat(C2,n,m,outfile);

     C_DCOPY(m*m,Atemp[0],1,A[0],1);
     for (int a = 0; a < m; a++)
        for (int b = 0; b < a; b++)
            A[a][b] = 0.0;

     fprintf(outfile, " AU:\n");
     print_mat(A,m,m,outfile);

     C_DSYMM('L','U',m,n,alpha, A[0], m, B1[0], n, beta, C1[0], n);

     fprintf(outfile," A*B: \n");
     print_mat(C1,m,n,outfile);

     C_DSYMM('R','U',n,m,alpha, A[0], m, B2[0], m, beta, C2[0], m);

     fprintf(outfile," B*A: \n");
     print_mat(C2,n,m,outfile);
     fflush(outfile);
     **/
}
double SAPT_DFT::print_results() {
    return 0.0;
}
void SAPT_DFT::allocate_arrays() {

    int naux = calc_info_.nri;
    // Allocation of UCHF response params
    X0_A_ = block_matrix(naux,naux);
    X0_B_ = block_matrix(naux,naux);
    // Allocation of TDDFT response params
    XC_A_ = block_matrix(naux,naux);
    XC_B_ = block_matrix(naux,naux);
}
void SAPT_DFT::free_arrays() {
    free_block(X0_A_);
    free_block(X0_B_);
    free_block(XC_A_);
    free_block(XC_B_);
}
void SAPT_DFT::compute_J() {

    int naux = calc_info_.nri;
    IntegralFactory rifactory_J(ribasis_, zero_, ribasis_, zero_);

    TwoBodyAOInt* Jint = rifactory_J.eri();
    //double **calc_info_.J = block_matrix(ribasis_->nbf(), ribasis_->nbf());
    J_ = block_matrix(ribasis_->nbf(), ribasis_->nbf());
    const double *Jbuffer = Jint->buffer();

    int index = 0;

     for (int MU=0; MU < ribasis_->nshell(); ++MU) {
        int nummu = ribasis_->shell(MU)->nfunction();

        for (int NU=0; NU <= MU; ++NU) {
            int numnu = ribasis_->shell(NU)->nfunction();

            Jint->compute_shell(MU, 0, NU, 0);

            index = 0;
            for (int mu=0; mu < nummu; ++mu) {
                int omu = ribasis_->shell(MU)->function_index() + mu;

                for (int nu=0; nu < numnu; ++nu, ++index) {
                    int onu = ribasis_->shell(NU)->function_index() + nu;


                    J_[omu][onu] = Jbuffer[index];
                    J_[onu][omu] = Jbuffer[index];
                }
            }
        }
    }



    Jinv_ = block_matrix(naux,naux);
    C_DCOPY(naux*(ULI)naux, J_[0], 1, Jinv_[0], 1);

    // Cholesky factorization
    int stat = C_DPOTRF('U',naux,Jinv_[0],naux);
    // Inverse
    stat = C_DPOTRI('U',naux,Jinv_[0],naux);

    // Rexpand J^-1 (Upper triangle)
    for (int P = 0; P < naux; P++)
        for (int Q = P+1; Q < naux; Q++)
            Jinv_[P][Q] = Jinv_[Q][P];

    if (debug_) {
        fprintf(outfile,"J:\n");
        print_mat(J_,naux,naux,outfile);
        fprintf(outfile,"J^-1:\n");
        print_mat(Jinv_,naux,naux,outfile);
    }
}
void SAPT_DFT::compute_Q() {

    int naux_quad = quadbasis_->nbf();
    IntegralFactory rifactory_J(quadbasis_, zero_, quadbasis_, zero_);

    TwoBodyAOInt* Jint = rifactory_J.eri();
    //double **calc_info_.J = block_matrix(ribasis_->nbf(), ribasis_->nbf());
    Q_ = block_matrix(quadbasis_->nbf(), quadbasis_->nbf());
    const double *Jbuffer = Jint->buffer();

    int index = 0;

     for (int MU=0; MU < quadbasis_->nshell(); ++MU) {
        int nummu = quadbasis_->shell(MU)->nfunction();

        for (int NU=0; NU <= MU; ++NU) {
            int numnu = quadbasis_->shell(NU)->nfunction();

            Jint->compute_shell(MU, 0, NU, 0);

            index = 0;
            for (int mu=0; mu < nummu; ++mu) {
                int omu = quadbasis_->shell(MU)->function_index() + mu;

                for (int nu=0; nu < numnu; ++nu, ++index) {
                    int onu = quadbasis_->shell(NU)->function_index() + nu;


                    Q_[omu][onu] = Jbuffer[index];
                    Q_[onu][omu] = Jbuffer[index];
                }
            }
        }
    }



    Qinv_ = block_matrix(naux_quad,naux_quad);
    C_DCOPY(naux_quad*(ULI)naux_quad, Q_[0], 1, Qinv_[0], 1);

    // Cholesky factorization
    int stat = C_DPOTRF('U',naux_quad,Qinv_[0],naux_quad);
    // Inverse
    stat = C_DPOTRI('U',naux_quad,Qinv_[0],naux_quad);

    // Rexpand J^-1 (Upper triangle)
    for (int P = 0; P < naux_quad; P++)
        for (int Q = P+1; Q < naux_quad; Q++)
            Qinv_[P][Q] = Qinv_[Q][P];

    if (debug_) {
        fprintf(outfile,"Q:\n");
        print_mat(Q_,naux_quad,naux_quad,outfile);
        fprintf(outfile,"Q^-1:\n");
        print_mat(Qinv_,naux_quad,naux_quad,outfile);
    }
}
void SAPT_DFT::get_quad_basis()
{
    if (options_.get_str("RI_BASIS_QUAD") == "") {
        quadbasis_ = ribasis_;
    } else {
        shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(
          options_.get_str("BASIS_PATH")));
        quadbasis_ = BasisSet::construct(parser, molecule_, options_.get_str(
          "RI_BASIS_QUAD"));
    }
}
void SAPT_DFT::compute_d()
{
    compute_Q();

    // Get my bearings
    int naux = calc_info_.nri;
    int naux_quad = quadbasis_->nbf();
    int norbs = calc_info_.nso;
    int nmo = calc_info_.nso;
    int noccA = calc_info_.noccA;
    int noccB = calc_info_.noccB;
    int nocc = noccA + noccB;
    int nvirA = calc_info_.nvirA;
    int nvirB = calc_info_.nvirB;
    int nmoA = noccA+nvirA;
    int nmoB = noccB+nvirB;
    double **CA = calc_info_.CA;
    double **CB = calc_info_.CB;

    double** DA = block_matrix(norbs,norbs);
    C_DGEMM('N','T',norbs,norbs,noccA,1.0,&CA[0][0],nmo,&CA[0][0],nmo,0.0,&DA[0][0],norbs);
    double** DB = block_matrix(norbs,norbs);
    C_DGEMM('N','T',norbs,norbs,noccB,1.0,&CB[0][0],nmo,&CB[0][0],nmo,0.0,&DB[0][0],norbs);

    if (debug_) {
        fprintf(outfile,"CA:\n");
        print_mat(CA,norbs,nmoA,outfile);
        fprintf(outfile,"CB:\n");
        print_mat(CB,norbs,nmoA,outfile);
        fprintf(outfile,"DA:\n");
        print_mat(DA,norbs,norbs,outfile);
        fprintf(outfile,"DB:\n");
        print_mat(DB,norbs,norbs,outfile);
    }

    WA_ = block_matrix(naux,naux);
    WB_ = block_matrix(naux,naux);

    double* c_A = init_array(naux_quad);
    double* c_B = init_array(naux_quad);
    d_A_ = init_array(naux_quad);
    d_B_ = init_array(naux_quad);

    // =============== Density Coefficients ================//

    IntegralFactory rifactory(basisset_, basisset_, quadbasis_, zero_);

    const double *buffer;
    shared_ptr<TwoBodyAOInt> eri = shared_ptr<TwoBodyAOInt>(rifactory.eri());
    buffer = eri->buffer();

    // Indexing for three-index AO tensors
    int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
    int index;

    int max_P_shell = 0;
    for (Pshell=0; Pshell < quadbasis_->nshell(); ++Pshell)
        if (quadbasis_->shell(Pshell)->nfunction() > max_P_shell)
            max_P_shell = quadbasis_->shell(Pshell)->nfunction();

    double** Amn = block_matrix(max_P_shell,norbs*(ULI)norbs);

    // A bit naive at the moment (no sieves or threading)
    for (Pshell=0; Pshell < quadbasis_->nshell(); ++Pshell) {
        numP = quadbasis_->shell(Pshell)->nfunction();
        for (MU=0; MU < basisset_->nshell(); ++MU) {
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU < basisset_->nshell(); ++NU) {
                numnu = basisset_->shell(NU)->nfunction();
                eri->compute_shell(MU, NU, Pshell, 0);
                for (mu=0 ; mu < nummu; ++mu) {
                    omu = basisset_->shell(MU)->function_index() + mu;
                    for (nu=0; nu < numnu; ++nu) {
                        onu = basisset_->shell(NU)->function_index() + nu;
                        for (P=0; P < numP; ++P) {
                            Amn[P][omu*norbs+onu] = buffer[mu*numnu*numP+nu*numP+P];
                        }
                    }
                }
            }
        }
        C_DGEMV('n', numP, norbs*(ULI)norbs, 1.0, Amn[0], norbs*(ULI)norbs, &DA[0][0], 1, 0.0, \
            &c_A[quadbasis_->shell(Pshell)->function_index()], 1);
        C_DGEMV('n', numP, norbs*(ULI)norbs, 1.0, Amn[0], norbs*(ULI)norbs, &DB[0][0], 1, 0.0, \
            &c_B[quadbasis_->shell(Pshell)->function_index()], 1);
    }

    // d = J^-1 c
    C_DGEMV('n', naux_quad, naux_quad, 1.0, Qinv_[0], naux_quad, c_A, 1, 0.0, d_A_, 1);
    C_DGEMV('n', naux_quad, naux_quad, 1.0, Qinv_[0], naux_quad, c_B, 1, 0.0, d_B_, 1);

    if (debug_) {
        fprintf(outfile, "  c_A Coefficients:\n");
        for (P = 0; P < naux_quad; P++)
            fprintf(outfile, "    P = %d: c = %20.14f\n", \
                P, c_A[P]);
        fprintf(outfile, "  c_B Coefficients:\n");
        for (P = 0; P < naux_quad; P++)
            fprintf(outfile, "    P = %d: c = %20.14f\n", \
                P, c_B[P]);

        fprintf(outfile, "  d_A Coefficients:\n");
        for (P = 0; P < naux_quad; P++)
            fprintf(outfile, "    P = %d: d = %20.14f\n", \
                P, d_A_[P]);
        fprintf(outfile, "  d_B Coefficients:\n");
        for (P = 0; P < naux_quad; P++)
            fprintf(outfile, "    P = %d: d = %20.14f\n", \
                P, d_B_[P]);
    }

    free(c_A);
    free(c_B);

    // Test by forming Coulomb matrix
    if (debug_) {
        double **JA = block_matrix(norbs,norbs);
        double **JB = block_matrix(norbs,norbs);
        for (Pshell=0; Pshell < quadbasis_->nshell(); ++Pshell) {
            numP = quadbasis_->shell(Pshell)->nfunction();
            for (MU=0; MU < basisset_->nshell(); ++MU) {
                nummu = basisset_->shell(MU)->nfunction();
                for (NU=0; NU < basisset_->nshell(); ++NU) {
                    numnu = basisset_->shell(NU)->nfunction();
                    eri->compute_shell(MU, NU, Pshell, 0);
                    for (mu=0 ; mu < nummu; ++mu) {
                        omu = basisset_->shell(MU)->function_index() + mu;
                        for (nu=0; nu < numnu; ++nu) {
                            onu = basisset_->shell(NU)->function_index() + nu;
                            for (P=0; P < numP; ++P) {
                                Amn[P][omu*norbs+onu] = buffer[mu*numnu*numP+nu*numP+P];
                            }
                        }
                    }
                }
            }
            C_DGEMV('t', numP, norbs*(ULI)norbs, 1.0, Amn[0], norbs*(ULI)norbs,
                &d_A_[quadbasis_->shell(Pshell)->function_index()], 1, 1.0, JA[0], 1);
            C_DGEMV('t', numP, norbs*(ULI)norbs, 1.0, Amn[0], norbs*(ULI)norbs,
                &d_B_[quadbasis_->shell(Pshell)->function_index()], 1, 1.0, JB[0], 1);
        }
        fprintf(outfile, "  Trial Coulomb (A):\n");
        print_mat(JA, norbs, norbs, outfile);
        fprintf(outfile, "  Trial Coulomb (B):\n");
        print_mat(JB, norbs, norbs, outfile);
        free_block(JA);
        free_block(JB);
    }
    // End Test

    free_block(Amn);
    free_block(Q_);
    free_block(Qinv_);

    // =============== Numerical Test ===============//
    if (debug_) {

        double** WtestA = block_matrix(naux, naux);
        double** WptestA = block_matrix(naux, naux);
        double** WtestB = block_matrix(naux, naux);
        double** WptestB = block_matrix(naux, naux);

        shared_ptr<Matrix> DmatA(new Matrix("DA", 1, &norbs, &norbs));
        shared_ptr<Matrix> DmatB(new Matrix("DB", 1, &norbs, &norbs));
        for (int mu = 0; mu < norbs; mu ++)
            for (int nu = 0; nu < norbs; nu ++) {
                DmatA->set(0, mu, nu, DA[mu][nu]);
                DmatB->set(0, mu, nu, DB[mu][nu]);
            }

        DmatA->print(outfile);
        DmatB->print(outfile);

        shared_ptr<Properties> props(new Properties(basisset_, options_.get_int("N_BLOCK")));
        shared_ptr<BasisPoints> points(new BasisPoints(ribasis_, options_.get_int("N_BLOCK")));
        shared_ptr<psi::scf::Integrator> integrator(new psi::scf::Integrator(basisset_->molecule(), options_));

        const double* rhog = props->getDensity();
        double **basisg = points->getPoints();

        double rho_A = 0.0;
        double rho_B = 0.0;

        for (integrator->reset(); !integrator->isDone(); ) {
            shared_ptr<GridBlock> q = integrator->getNextBlock();
            int ntrue = q->getTruePoints();
            double* wg = q->getWeights();

            props->computeProperties(q, DmatA);
            points->computePoints(q);

            for (int grid_index = 0; grid_index<ntrue; grid_index++) {

               rho_A += wg[grid_index] * rhog[grid_index];
               for (int P = 0; P < naux; P++)
                   for (int Q = 0; Q < naux; Q++) {
                       WtestA[P][Q] += wg[grid_index] * rhog[grid_index]* \
                           basisg[grid_index][P] * basisg[grid_index][Q];
                   }
                if (rhog[grid_index] > 1E-20)
                    for (int P = 0; P < naux; P++)
                        for (int Q = 0; Q < naux; Q++) {
                            WptestA[P][Q] += wg[grid_index] * -C_x_*(8.0/9.0)*pow(rhog[grid_index],-2.0/3.0)* \
                                basisg[grid_index][P] * basisg[grid_index][Q];
                        }
            }
        }
        for (integrator->reset(); !integrator->isDone(); ) {
            shared_ptr<GridBlock> q = integrator->getNextBlock();
            int ntrue = q->getTruePoints();
            double* wg = q->getWeights();

            props->computeProperties(q, DmatB);
            points->computePoints(q);

            for (int grid_index = 0; grid_index<ntrue; grid_index++) {

               rho_B += wg[grid_index] * rhog[grid_index];
               for (int P = 0; P < naux; P++)
                   for (int Q = 0; Q < naux; Q++) {
                       WtestB[P][Q] += wg[grid_index] * rhog[grid_index]* \
                           basisg[grid_index][P] * basisg[grid_index][Q];
                   }
                if (rhog[grid_index] > 1E-20)
                    for (int P = 0; P < naux; P++)
                        for (int Q = 0; Q < naux; Q++) {
                            WptestB[P][Q] += wg[grid_index] * -C_x_*(8.0/9.0)*pow(rhog[grid_index],-2.0/3.0)* \
                                basisg[grid_index][P] * basisg[grid_index][Q];
                        }
            }
        }

        fprintf(outfile, "  Numerical Densities: A = %14.10f, B = %14.10f.\n", rho_A, rho_B);
        fprintf(outfile, "  Numerical W_A (P|\\rho|Q)\n");
        print_mat(WtestA, naux, naux, outfile);
        fprintf(outfile, "  Numerical W_B (P|\\rho|Q)\n");
        print_mat(WtestB, naux, naux, outfile);
        fprintf(outfile, "  Numerical W_A (P|f_xc(\\rho)|Q)\n");
        print_mat(WptestA, naux, naux, outfile);
        fprintf(outfile, "  Numerical W_B (P|f_xc(\\rho)|Q)\n");
        print_mat(WptestB, naux, naux, outfile);

        //C_DCOPY(naux*(ULI)naux, WtestA[0], 1, WA_[0], 1);
        //C_DCOPY(naux*(ULI)naux, WtestB[0], 1, WB_[0], 1);
        free_block(WtestA);
        free_block(WptestA);
        free_block(WtestB);
        free_block(WptestB);
    }
    free_block(DA);
    free_block(DB);
}
void SAPT_DFT::compute_W()
{
    int naux = calc_info_.nri;
    int naux_quad = quadbasis_->nbf();

    C_x_ = 3.0/8.0*pow(3.0,1.0/3.0)*pow(4.0,2.0/3.0)*pow(M_PI,-1.0/3.0);
    compute_d();

    // =============== Heavy Three-Index Tensor ================//

    IntegralFactory PQRfactory(ribasis_, ribasis_, quadbasis_, zero_);

    int max_P_shell = 0;
    for (int Pshell=0; Pshell < ribasis_->nshell(); ++Pshell)
        if (ribasis_->shell(Pshell)->nfunction() > max_P_shell)
            max_P_shell = ribasis_->shell(Pshell)->nfunction();

    double** PQR = block_matrix(max_P_shell*max_P_shell,naux_quad);
    double* TempA = init_array(max_P_shell*max_P_shell);
    double* TempB = init_array(max_P_shell*max_P_shell);

    shared_ptr<ThreeCenterOverlapInt> o3(PQRfactory.overlap_3c());
    const double* buffer = o3->buffer();

    int P, numP, Q, numQ, R, numR, q, p, r, op, oq, oR;
    // A bit naive at the moment (no sieves or threading)
    for (P=0; P < ribasis_->nshell(); ++P) {
        numP = ribasis_->shell(P)->nfunction();
        for (Q=0; Q< ribasis_->nshell(); ++Q) {
            numQ = ribasis_->shell(Q)->nfunction();
            for (R=0; R < quadbasis_->nshell();  ++R) {
                numR = quadbasis_->shell(R)->nfunction();
                o3->compute_shell(P, Q, R);
                for (p=0 ; p < numP; ++p) {
                    op = ribasis_->shell(P)->function_index() + p;
                    for (q=0; q < numQ; ++q) {
                        oq = ribasis_->shell(Q)->function_index() + q;
                        for (r=0; r < numR; ++r) {
                            oR = quadbasis_->shell(R)->function_index() + r;
                            if (debug_)
                                fprintf(outfile, "(%4d %4d %4d) : %14.10f\n", op, oq, oR, buffer[p*numQ*numR + q*numR + r]);
                            PQR[p*numQ + q][r + quadbasis_->shell(R)->function_index()] = \
                                buffer[p*numQ*numR+q*numR+r];
                        }
                    }
                }
            }
            // Multiply
            //if (debug_) {
            //    fprintf(outfile, "  PQ Shell: (%3d %3d | R):\n", P, Q);
            //    print_mat(PQR, numP*numQ, naux, outfile);
            //}
            C_DGEMV('n', numP*numQ, naux_quad, 1.0, PQR[0], naux_quad, d_A_, 1, 0.0, TempA, 1);
            C_DGEMV('n', numP*numQ, naux_quad, 1.0, PQR[0], naux_quad, d_B_, 1, 0.0, TempB, 1);
            //if (debug_) {
            //    fprintf(outfile, "  Bits of W for  shells (%3d|\\rho |%3d): \n", P, Q);
            //}
            for (p=0 ; p < numP; ++p) {
                op = ribasis_->shell(P)->function_index() + p;
                for (q=0; q < numQ; ++q) {
                    oq = ribasis_->shell(Q)->function_index() + q;
                    //if (debug_) {
                    //    fprintf(outfile, "    (%d | \\rho |%d)  = %14.10f, %14.10f\n", op, oq, TempA[p*numQ + q], TempB[p*numQ + q]);
                    //}
                    WA_[op][oq] = TempA[p*numQ + q];
                    WB_[op][oq] = TempB[p*numQ + q];
                }
            }
        }
    }

    if (debug_) {
        fprintf(outfile, "  W_A (P|\\rho|Q):\n");
        print_mat(WA_, naux, naux, outfile);
        fprintf(outfile, "  W_B (P|\\rho|Q):\n");
        print_mat(WB_, naux, naux, outfile);
    }
    free_block(PQR);
    free(TempA);
    free(TempB);
    free(d_A_);
    free(d_B_);


    // =============== M -> Mprime ================//

    IntegralFactory rifactory_JS(ribasis_, ribasis_, zero_,  zero_);
    shared_ptr<OneBodySOInt> Sint = shared_ptr<OneBodySOInt>(rifactory_JS.so_overlap());
    MatrixFactory matJS;
    matJS.init_with(1,&naux,&naux);

    //Put the integrals in a good old double**
    SharedMatrix S_J = shared_ptr<Matrix>(matJS.create_matrix("S_J"));
    //Compute those integrals
    Sint->compute(S_J);
    if (debug_)
        S_J->print(outfile);
    double** S = S_J->to_block_matrix();

    double **V = block_matrix(naux, naux);
    double **Vtemp = block_matrix(naux, naux);
    double **Vcopy = block_matrix(naux, naux);
    C_DCOPY(naux*(ULI)naux,S[0],1,Vtemp[0],1);
    // Form V : V'SV = 1
    // First, diagonalize V
    // the C_DSYEV call replaces the original matrix J with its eigenvectors
    double* eigval = init_array(naux);
    int lwork = naux * 3;
    double* work = init_array(lwork);
    int stat = C_DSYEV('v','u',naux,Vtemp[0],naux,eigval, work,lwork);
    if (stat != 0) {
        fprintf(outfile, "C_DSYEV failed\n");
        exit(PSI_RETURN_FAILURE);
    }

    C_DCOPY(naux*(ULI)naux,Vtemp[0],1,Vcopy[0],1);

    for (int ind=0; ind<naux; ind++) {
        if (eigval[ind] < 1.0E-10)
            eigval[ind] = 0.0;
        else {
            eigval[ind] = - 1.0 / sqrt(eigval[ind]);
        }
        // scale one set of eigenvectors by the diagonal elements j^{-1/2}
        C_DSCAL(naux, eigval[ind], Vcopy[ind], 1);
    }
    free(eigval);

    // S'^-1/2 = Sp_copy(T) * Sp
    C_DGEMM('t','n',naux,naux,naux,1.0,
            Vcopy[0],naux,Vtemp[0],naux,0.0,V[0],naux);

    free_block(Vcopy);
    free_block(Vtemp);

    if (debug_) {
        fprintf(outfile, "  V (S^-1/2):\n");
        print_mat(V, naux, naux, outfile);
    }

    double** Temp1 = block_matrix(naux,naux);

    C_DGEMM('T','N', naux, naux, naux, 1.0, V[0], naux, WA_[0], naux, \
        0.0, Temp1[0], naux);
    C_DGEMM('N','N', naux, naux, naux, 1.0, Temp1[0], naux, V[0], naux, \
        0.0, WA_[0], naux);
    C_DGEMM('T','N', naux, naux, naux, 1.0, V[0], naux, WB_[0], naux, \
        0.0, Temp1[0], naux);
    C_DGEMM('N','N', naux, naux, naux, 1.0, Temp1[0], naux, V[0], naux, \
        0.0, WB_[0], naux);

    if (debug_) {
        fprintf(outfile, "  Wp A (V'*W*V):\n");
        print_mat(WA_, naux, naux, outfile);
        fprintf(outfile, "  Wp B (V'*W*V):\n");
        print_mat(WB_, naux, naux, outfile);
    }

    // =============== Apply functional kernel ================//

    double** U = block_matrix(naux,naux);

    // Monomer A
    C_DCOPY(naux*(ULI)naux,WA_[0],1,U[0],1);
    // Form V : V'SV = 1
    // First, diagonalize V
    // the C_DSYEV call replaces the original matrix J with its eigenvectors
    double* lambda = init_array(naux);
    stat = C_DSYEV('v','u',naux,U[0],naux,lambda, work,lwork);
    if (stat != 0) {
        fprintf(outfile, "C_DSYEV failed\n");
        exit(PSI_RETURN_FAILURE);
    }

    if (debug_) {
        fprintf(outfile, "  U_A:\n");
        print_mat(U, naux, naux, outfile);
        fprintf(outfile, "  rho_A Coefficients:\n");
        for (P = 0; P < naux; P++)
            fprintf(outfile, "    P = %d: rho = %20.14f\n", \
                P, lambda[P]);
    }

    // Here's the functional
    C_DCOPY(naux*(ULI)naux,U[0],1,Temp1[0],1);
    bool warn = false;
    for (int k = 0; k < naux; k++) {
        if (lambda[k] < 1.0E-10) {
            warn = true;
            lambda[k] = 0.0;
        } else {
            lambda[k] = -8.0/9.0 * C_x_ * pow(lambda[k], -2.0/3.0);
        }
        C_DSCAL(naux, lambda[k], U[k], 1);
    }
    if (warn)
        fprintf(outfile, "  WARNING: Small/negative eigenvalue detected in (P|\\rho|Q)\n");
    if (debug_) {
        fprintf(outfile, "  lambda_A Coefficients:\n");
        for (P = 0; P < naux; P++)
            fprintf(outfile, "    P = %d: lambda = %20.14f\n", \
                P, lambda[P]);
    }

    // Now revert to M
    C_DGEMM('T','N',naux, naux, naux, 1.0, U[0], naux, Temp1[0], naux, \
        0.0, WA_[0], naux);
    C_DGEMM('N','N',naux, naux, naux, 1.0, V[0], naux, WA_[0], naux, \
        0.0, Temp1[0], naux);
    C_DGEMM('N','N',naux, naux, naux, 1.0, S[0], naux, Temp1[0], naux, \
        0.0, WA_[0], naux);
    C_DGEMM('N','T',naux, naux, naux, 1.0, WA_[0], naux, V[0], naux, \
        0.0, Temp1[0], naux);
    C_DGEMM('N','T',naux, naux, naux, 1.0, Temp1[0], naux, S[0], naux, \
        0.0, WA_[0], naux);

    // Monomer B
    C_DCOPY(naux*(ULI)naux,WB_[0],1,U[0],1);
    // Form V : V'SV = 1
    // First, diagonalize V
    // the C_DSYEV call replaces the original matrix J with its eigenvectors
    stat = C_DSYEV('v','u',naux,U[0],naux,lambda, work,lwork);
    if (stat != 0) {
        fprintf(outfile, "C_DSYEV failed\n");
        exit(PSI_RETURN_FAILURE);
    }

    if (debug_) {
        fprintf(outfile, "  U:\n");
        print_mat(U, naux, naux, outfile);
        fprintf(outfile, "  rho_B Coefficients:\n");
        for (P = 0; P < naux; P++)
            fprintf(outfile, "    P = %d: rho = %20.14f\n", \
                P, lambda[P]);
    }

    // Here's the functional
    C_DCOPY(naux*(ULI)naux,U[0],1,Temp1[0],1);
    warn = false;
    for (int k = 0; k < naux; k++) {
        if (lambda[k] < 1.0E-10) {
            warn = true;
            lambda[k] = 0.0;
        } else {
            lambda[k] = -8.0/9.0 * C_x_ * pow(lambda[k], -2.0/3.0);
        }
        C_DSCAL(naux, lambda[k], U[k], 1);
    }
    if (warn)
        fprintf(outfile, "  WARNING: Small/negative eigenvalue detected in (P|\\rho|Q)\n");
    if (debug_) {
        fprintf(outfile, "  lambda_B Coefficients:\n");
        for (P = 0; P < naux; P++)
            fprintf(outfile, "    P = %d: lambda = %20.14f\n", \
                P, lambda[P]);
    }

    // Now revert to M
    C_DGEMM('T','N',naux, naux, naux, 1.0, U[0], naux, Temp1[0], naux, \
        0.0, WB_[0], naux);
    C_DGEMM('N','N',naux, naux, naux, 1.0, V[0], naux, WB_[0], naux, \
        0.0, Temp1[0], naux);
    C_DGEMM('N','N',naux, naux, naux, 1.0, S[0], naux, Temp1[0], naux, \
        0.0, WB_[0], naux);
    C_DGEMM('N','T',naux, naux, naux, 1.0, WB_[0], naux, V[0], naux, \
        0.0, Temp1[0], naux);
    C_DGEMM('N','T',naux, naux, naux, 1.0, Temp1[0], naux, S[0], naux, \
        0.0, WB_[0], naux);

    // =============== Release Memory ================//

    free(lambda);
    free(work);

    free_block(Temp1);
    free_block(U);
    free_block(V);

    if (debug_) {
        fprintf(outfile, "  WA (P|f_xc(\\rho)|Q):\n");
        print_mat(WA_, naux, naux, outfile);
        fprintf(outfile, "  WB (P|f_xc(\\rho)|Q):\n");
        print_mat(WB_, naux, naux, outfile);
    }
    for (P = 0; P < naux; P++) {
        for (Q = 0; Q < naux; Q++) {
            WA_[P][Q] += J_[P][Q];
            WB_[P][Q] += J_[P][Q];
        }
    }
    if (debug_) {
        fprintf(outfile, "  WA (P|1/r + f_xc(\\rho)|Q):\n");
        print_mat(WA_, naux, naux, outfile);
        fprintf(outfile, "  WB (P|1/r + f_xc(\\rho)|Q):\n");
        print_mat(WB_, naux, naux, outfile);
    }
}
void SAPT_DFT::compute_X_0(double omega) {

    // Get my bearings
    int naux = calc_info_.nri;
    int norbs = calc_info_.nso;
    int noccA = calc_info_.noccA;
    int noccB = calc_info_.noccB;
    int nvirA = calc_info_.nvirA;
    int nvirB = calc_info_.nvirB;

    double lambda, eps_ia;
    double omega2 = omega*omega;
    double *eps_i, *eps_a;

    // =============== MONOMER A ================//
    if (debug_)
        fprintf(outfile, "  Monomer A X_0:\n");
    double** A_ints = block_matrix(naux, nvirA*(ULI)noccA);
    eps_i = calc_info_.evalsA;
    eps_a = calc_info_.evalsA + noccA;

    if (debug_) {
        fprintf(outfile,"  Epsilon occ:\n");
        for (int i = 0; i < noccA; i++)
            fprintf(outfile, "   i = %d: eps = %20.14f\n", i, eps_i[i]);
        fprintf(outfile,"  Epsilon virt:\n");
        for (int i = 0; i < nvirA; i++)
            fprintf(outfile, "   a = %d: eps = %20.14f\n", i, eps_a[i]);
    }

    // Read the AR DF integrals
    psio_address next_DF_AA = PSIO_ZERO;
    psio_->read(PSIF_SAPT_AA_DF_INTS, "AR Bare RI Integrals", (char*) \
        &A_ints[0][0], naux*nvirA*noccA*(ULI)sizeof(double), next_DF_AA, \
        &next_DF_AA);

    if (debug_) {
        fprintf(outfile, "  AR Bare RI Integrals (A|ia)\n");
        print_mat(A_ints, naux, nvirA*(ULI)noccA, outfile);
    }

    // Scale the products ia by \sqrt(lambda_ia)
    if (debug_)
        fprintf(outfile, "  Applying response, omega = %20.14f\n",omega);
    for (int i = 0; i < noccA; i++) {
        for (int a = 0; a < nvirA; a++) {
            eps_ia = eps_a[a] - eps_i[i];
            lambda = 4.0 * eps_ia / (eps_ia*eps_ia + omega2);
            if (debug_)
                fprintf(outfile, "   i = %d, a = %d, lambda = %20.14f\n", i, a, lambda);
            C_DSCAL(naux, sqrt(lambda), &A_ints[0][i*nvirA + a], nvirA*(ULI)noccA);
        }
    }

    if (debug_) {
        fprintf(outfile, "  Symmetric Response AR RI Integrals (A|ia)\n");
        print_mat(A_ints, naux, nvirA*(ULI)noccA, outfile);
    }

    // The only O(N^4) DGEMM in the whole thing
    C_DGEMM('N','T', naux, naux, noccA*(ULI)nvirA, 1.0, A_ints[0], noccA*(ULI)nvirA, \
        A_ints[0], noccA*(ULI)nvirA, 0.0, X0_A_[0], naux);

    if (debug_) {
        fprintf(outfile, "X_0^A:\n");
        print_mat(X0_A_, naux, naux, outfile);
    }

    free_block(A_ints);

    // =============== MONOMER B ================//
    //fprintf(outfile, "  Monomer B X_0:\n");
    double** B_ints = block_matrix(naux, nvirB*(ULI)noccB);
    eps_i = calc_info_.evalsB;
    eps_a = calc_info_.evalsB + noccB;

    if (debug_) {
        fprintf(outfile,"  Epsilon occ:\n");
        for (int i = 0; i < noccB; i++)
            fprintf(outfile, "   i = %d: eps = %20.14f\n", i, eps_i[i]);
        fprintf(outfile,"  Epsilon virt:\n");
        for (int i = 0; i < nvirB; i++)
            fprintf(outfile, "   a = %d: eps = %20.14f\n", i, eps_a[i]);
    }

    // Read the BS DF integrals
    psio_address next_DF_BB = PSIO_ZERO;
    psio_->read(PSIF_SAPT_BB_DF_INTS, "BS Bare RI Integrals", (char*) \
        &B_ints[0][0], naux*nvirB*noccB*(ULI)sizeof(double), next_DF_BB, \
        &next_DF_BB);

    if (debug_) {
        fprintf(outfile, "  BS Bare RI Integrals (A|ia)\n");
        print_mat(A_ints, naux, nvirB*(ULI)noccB, outfile);
    }

    // Scale the products ia by \sqrt(lambda_ia)
    if (debug_)
        fprintf(outfile, "  Applying response, omega = %20.14f\n",omega);
    for (int i = 0; i < noccB; i++) {
        for (int a = 0; a < nvirB; a++) {
            eps_ia = eps_a[a] - eps_i[i];
            lambda = 4.0 * eps_ia / (eps_ia*eps_ia + omega2);
            if (debug_)
                fprintf(outfile, "   i = %d, a = %d, lambda = %20.14f\n", i, a, lambda);
            C_DSCAL(naux, sqrt(lambda), &B_ints[0][i*nvirB + a], noccB*(ULI)nvirB);
        }
    }

    if (debug_) {
        fprintf(outfile, "  Symmetric Response BS RI Integrals (A|ia)\n");
        print_mat(B_ints, naux, nvirB*(ULI)noccB, outfile);
    }

    // OK, I lied
    C_DGEMM('N','T', naux, naux, noccB*(ULI)nvirB, 1.0, B_ints[0], noccB*(ULI)nvirB, \
        B_ints[0], noccB*(ULI)nvirB, 0.0, X0_B_[0], naux);

    if (debug_) {
        fprintf(outfile, "X_0^B:\n");
        print_mat(X0_B_, naux, naux, outfile);
    }

    free_block(B_ints);
}
void SAPT_DFT::compute_X_coup(double omega) {

    // Get my bearings
    int naux = calc_info_.nri;

    double ** Temp1 = block_matrix(naux,naux);
    double ** Temp2 = block_matrix(naux,naux);
    double ** Temp3 = block_matrix(naux,naux);

    // ========= XC_A ========== //
    // Form S - X_0 S^-1 W
    C_DGEMM('N','N', naux, naux, naux, 1.0, X0_A_[0], naux, Jinv_[0], naux, \
        0.0, Temp1[0], naux);
    C_DGEMM('N','N', naux, naux, naux, 1.0, Temp1[0], naux, WA_[0], naux, \
        0.0, Temp2[0], naux);
    for (ULI k = 0; k < naux*(ULI)naux; k++)
        Temp1[0][k] = J_[0][k] - Temp2[0][k];

    // Invert S - X_0 S^-1 W
    // LU factorization
    int* piv = init_int_array(naux);
    int lwork = 3 * naux;
    double* work = init_array(lwork);
    int stat = C_DGETRF(naux,naux,Temp1[0],naux, piv);
    // Inverse
    stat = stat | C_DGETRI(naux,Temp1[0],naux,piv,work,lwork);

    if (stat) {
        fprintf(outfile, "  Call Rob at CCMST and tell him MP2C's XC Dyson Equation Failed.\n");
        fprintf(outfile, "  (Buy him a beer when he fixes it).\n");
        exit(PSI_RETURN_FAILURE);
    }

    // Form X_0 S^1 W (S - X_0 S^-1 W) ^ -1
    C_DGEMM('N','N', naux, naux, naux, 1.0, Temp2[0], naux, Temp1[0], naux, \
        0.0, Temp3[0], naux);

    // Form XC = X_0 S^1 W (S - X_0 S^-1 W) ^ -1 X_0
    C_DGEMM('N','N', naux, naux, naux, 1.0, Temp3[0], naux, X0_A_[0], naux, \
        0.0, XC_A_[0], naux);

    if (debug_) {
        fprintf(outfile, "  XC_A:\n");
        print_mat(XC_A_, naux, naux, outfile);
    }

    // ========= XC_B ========== //
    // Form S - X_0 S^-1 W
    C_DGEMM('N','N', naux, naux, naux, 1.0, X0_B_[0], naux, Jinv_[0], naux, \
        0.0, Temp1[0], naux);
    C_DGEMM('N','N', naux, naux, naux, 1.0, Temp1[0], naux, WB_[0], naux, \
        0.0, Temp2[0], naux);
    for (ULI k = 0; k < naux*(ULI)naux; k++)
        Temp1[0][k] = J_[0][k] - Temp2[0][k];

    // Invert S - X_0 S^-1 W
    // LU factorization
    stat = C_DGETRF(naux,naux,Temp1[0],naux, piv);
    // Inverse
    stat = stat | C_DGETRI(naux,Temp1[0],naux,piv,work,lwork);
    free(work);
    free(piv);

    if (stat) {
        fprintf(outfile, "  Call Rob at CCMST and tell him MP2C's Dyson Equation Failed.\n");
        fprintf(outfile, "  (Buy him a beer when he fixes it).\n");
        exit(PSI_RETURN_FAILURE);
    }

    // Form X_0 S^1 W (S - X_0 S^-1 W) ^ -1
    C_DGEMM('N','N', naux, naux, naux, 1.0, Temp2[0], naux, Temp1[0], naux, \
        0.0, Temp3[0], naux);

    // Form XC = X_0 S^1 W (S - X_0 S^-1 W) ^ -1 X_0
    C_DGEMM('N','N', naux, naux, naux, 1.0, Temp3[0], naux, X0_B_[0], naux, \
        0.0, XC_B_[0], naux);

    if (debug_) {
        fprintf(outfile, "  XC_B:\n");
        print_mat(XC_B_, naux, naux, outfile);
    }

    free_block(Temp1);
    free_block(Temp2);
    free_block(Temp3);
}
double SAPT_DFT::compute_UCHF_disp()
{

    int naux = calc_info_.nri;

    double ** Temp1 = block_matrix(naux,naux);
    double ** Temp2 = block_matrix(naux,naux);

    // Form C_A
    C_DGEMM('T','T', naux, naux, naux, 1.0, X0_A_[0], naux, Jinv_[0], naux, \
        0.0, Temp1[0], naux);

    // Form C_B
    C_DGEMM('N','N', naux, naux, naux, 1.0, Jinv_[0], naux, X0_B_[0], naux, \
        0.0, Temp2[0], naux);

    double contribution = -1.0/(2.0*M_PI)*C_DDOT(naux*(ULI)naux, Temp1[0], 1, Temp2[0],1);

    if (debug_) {
        fprintf(outfile,"  Q_A:\n");
        print_mat(Temp1, naux, naux, outfile);

        fprintf(outfile,"  Q_B:\n");
        print_mat(Temp2, naux, naux, outfile);

        fprintf(outfile, "  Contribution: %20.14f\n", contribution);
    }

    free_block(Temp1);
    free_block(Temp2);

    return contribution;
}
double SAPT_DFT::compute_TDDFT_disp()
{

    int naux = calc_info_.nri;

    double ** Temp1 = block_matrix(naux,naux);
    double ** Temp2 = block_matrix(naux,naux);

    // Form C_A
    C_DGEMM('T','T', naux, naux, naux, 1.0, XC_A_[0], naux, Jinv_[0], naux, \
        0.0, Temp1[0], naux);

    // Form C_B
    C_DGEMM('N','N', naux, naux, naux, 1.0, Jinv_[0], naux, XC_B_[0], naux, \
        0.0, Temp2[0], naux);

    double contribution = -1.0/(2.0*M_PI)*C_DDOT(naux*(ULI)naux, Temp1[0], 1, Temp2[0],1);

    if (debug_) {
        fprintf(outfile,"  Q_A:\n");
        print_mat(Temp1, naux, naux, outfile);

        fprintf(outfile,"  Q_B:\n");
        print_mat(Temp2, naux, naux, outfile);

        fprintf(outfile, "  Contribution: %20.14f\n", contribution);
    }

    free_block(Temp1);
    free_block(Temp2);

    return contribution;
}

}}
