/* 
 *  SAPT.CC 
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

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include "sapt.h"
#include "structs.h"

#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/molecule.h>

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace sapt {

SAPT::SAPT(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
    : Wavefunction(options, psio, chkpt)
{
    setup_sapt();
}

SAPT::~SAPT()
{
    cleanup_calc_info();
}
void SAPT::setup_sapt()
{
    get_params();
    get_ribasis();
    get_calc_info();
    oetrans();
    print_header();
}
void SAPT::get_ribasis()
{
   ribasis_ = shared_ptr<BasisSet>(new BasisSet(chkpt_, "DF_BASIS_SAPT"));
   zero_ = BasisSet::zero_basis_set();
   calc_info_.nri = ribasis_->nbf();  
   calc_info_.nrio = ribasis_->nbf() + 3;  
}
void SAPT::get_params()
{
    //LR Ints?
    params_.lr_ints = options_.get_bool("LR_INTS");
    
    //Files restarts 
    params_.df_restart = options_.get_bool("DF_RESTART");
    params_.t2_restart = options_.get_bool("T2_RESTART");
    params_.logfile = options_.get_bool("LOGFILE");
    if (params_.logfile) 
      params_.logfilename = fopen("sapt.log","w");
 
    //CPHF convergence parameters
    params_.e_conv = pow(10.0,-options_.get_int("E_CONVERGE"));
    params_.d_conv = pow(10.0,-options_.get_int("D_CONVERGE"));
    params_.maxiter = options_.get_int("MAXITER");
    params_.diisvec = options_.get_int("DIISVECS");
    
    //Print 
    params_.print = options_.get_int("PRINT");

    //Schwarz cutoff
    params_.schwarz = options_.get_double("SCHWARZ_CUTOFF");
}
void SAPT::get_calc_info()
{
  params_.memory = (double) memory_;

  psio_open(PSIF_SAPT_DIMER,PSIO_OPEN_OLD);

  int errcod = 0;
  errcod = psio_read_entry(PSIF_SAPT_DIMER,"Dimer NSO",(char *) &calc_info_.nso, sizeof(int));
  errcod = psio_read_entry(PSIF_SAPT_DIMER,"Dimer NMO",(char *) &calc_info_.nmo, sizeof(int));
  errcod = psio_read_entry(PSIF_SAPT_DIMER,"Dimer HF Energy",(char *) &calc_info_.eHF_D, sizeof(double));
  errcod = psio_read_entry(PSIF_SAPT_DIMER,"Dimer Nuclear Repulsion Energy",(char *) &calc_info_.enuc_D,
    sizeof(double));

  calc_info_.nsotri = calc_info_.nso*(calc_info_.nso+1)/2;
  calc_info_.nmotri = calc_info_.nmo*(calc_info_.nmo+1)/2;
  calc_info_.ntei = calc_info_.nsotri*(calc_info_.nsotri+1)/2;

  /* Store overlap integrals */
  calc_info_.S = init_array(calc_info_.nsotri);
  errcod = psio_read_entry(PSIF_SAPT_DIMER,"Dimer Overlap Integrals",(char *) &calc_info_.S[0],
    sizeof(double)*calc_info_.nsotri);

  psio_close(PSIF_SAPT_DIMER,1);

  calc_info_.ioff = (int *) malloc (calc_info_.nsotri * sizeof(int));
  calc_info_.index2i = (int *) malloc (calc_info_.nsotri * sizeof(int));
  calc_info_.index2j = (int *) malloc (calc_info_.nsotri * sizeof(int));

  calc_info_.ioff[0] = 0; // Create ioff array

  for (int i=1; i < calc_info_.nsotri; i++) {
    calc_info_.ioff[i] = calc_info_.ioff[i-1] + i;
  }

  for (int i=0; i<calc_info_.nso; i++) {
    for (int j=0; j<=i; j++) {
      calc_info_.index2i[INDEX(i,j)] = i;
      calc_info_.index2j[INDEX(i,j)] = j;
    }}

  psio_open(PSIF_SAPT_MONOMERA,PSIO_OPEN_OLD);

  errcod = psio_read_entry(PSIF_SAPT_MONOMERA,"Monomer NSO",(char *) &calc_info_.nso, sizeof(int));
  errcod = psio_read_entry(PSIF_SAPT_MONOMERA,"Monomer NMO",(char *) &calc_info_.nmo, sizeof(int));
  errcod = psio_read_entry(PSIF_SAPT_MONOMERA,"Monomer NOCC",(char *) &calc_info_.noccA, sizeof(int));
  errcod = psio_read_entry(PSIF_SAPT_MONOMERA,"Monomer NVIR",(char *) &calc_info_.nvirA, sizeof(int));
  errcod = psio_read_entry(PSIF_SAPT_MONOMERA,"Monomer Number of Electrons",(char *) &calc_info_.NA,
    sizeof(int));
  errcod = psio_read_entry(PSIF_SAPT_MONOMERA,"Monomer HF Energy",(char *) &calc_info_.eHF_A,
    sizeof(double));
  errcod = psio_read_entry(PSIF_SAPT_MONOMERA,"Monomer Nuclear Repulsion Energy",
    (char *) &calc_info_.enuc_A, sizeof(double));


  calc_info_.evalsA = init_array(calc_info_.nmo);
  errcod = psio_read_entry(PSIF_SAPT_MONOMERA,"Monomer HF Eigenvalues",(char *) &(calc_info_.evalsA[0]),
    sizeof(double)*calc_info_.nmo);

  calc_info_.VA = init_array(calc_info_.nsotri);
  errcod = psio_read_entry(PSIF_SAPT_MONOMERA,"Monomer Nuclear Attraction Integrals",
    (char *) &(calc_info_.VA[0]), sizeof(double)*calc_info_.nsotri);

  calc_info_.CA = block_matrix(calc_info_.nso,calc_info_.nmo);
  errcod = psio_read_entry(PSIF_SAPT_MONOMERA,"Monomer HF Coefficients",(char *) &(calc_info_.CA[0][0]),
    sizeof(double)*calc_info_.nmo*calc_info_.nso);

  psio_close(PSIF_SAPT_MONOMERA,1);

  psio_open(PSIF_SAPT_MONOMERB,PSIO_OPEN_OLD);

  errcod = psio_read_entry(PSIF_SAPT_MONOMERB,"Monomer NOCC",(char *) &calc_info_.noccB, sizeof(int));
  errcod = psio_read_entry(PSIF_SAPT_MONOMERB,"Monomer NVIR",(char *) &calc_info_.nvirB, sizeof(int));
  errcod = psio_read_entry(PSIF_SAPT_MONOMERB,"Monomer Number of Electrons",(char *) &calc_info_.NB,
    sizeof(int));
  errcod = psio_read_entry(PSIF_SAPT_MONOMERB,"Monomer HF Energy",(char *) &calc_info_.eHF_B,
    sizeof(double));
  errcod = psio_read_entry(PSIF_SAPT_MONOMERB,"Monomer Nuclear Repulsion Energy",
    (char *) &calc_info_.enuc_B, sizeof(double));

  calc_info_.evalsB = init_array(calc_info_.nmo);
  errcod = psio_read_entry(PSIF_SAPT_MONOMERB,"Monomer HF Eigenvalues",(char *) &(calc_info_.evalsB[0]),
    sizeof(double)*calc_info_.nmo);

  calc_info_.VB = init_array(calc_info_.nsotri);
  errcod = psio_read_entry(PSIF_SAPT_MONOMERB,"Monomer Nuclear Attraction Integrals",
    (char *) &(calc_info_.VB[0]), sizeof(double)*calc_info_.nsotri);

  calc_info_.CB = block_matrix(calc_info_.nso,calc_info_.nmo);
  errcod = psio_read_entry(PSIF_SAPT_MONOMERB,"Monomer HF Coefficients",(char *) &(calc_info_.CB[0][0]),
    sizeof(double)*calc_info_.nmo*calc_info_.nso);

  psio_close(PSIF_SAPT_MONOMERB,1);

  results_.hf_int = calc_info_.eHF_D - calc_info_.eHF_A - calc_info_.eHF_B;

}
void SAPT::cleanup_calc_info()
{
  free_block(calc_info_.S_AB);
  free_block(calc_info_.VABB);
  free_block(calc_info_.VAAB);
  free_block(calc_info_.VBAA);
  free_block(calc_info_.VBAB);
  free_block(calc_info_.CHFA);
  free_block(calc_info_.CHFB);
  free_block(calc_info_.WABS);
  free_block(calc_info_.WBAR);

  free(calc_info_.diagAA);
  free(calc_info_.diagBB);
  free(calc_info_.evalsA);
  free(calc_info_.evalsB);
  free(calc_info_.ioff);
  free(calc_info_.index2i);
  free(calc_info_.index2j);
  if (params_.logfile) fclose(params_.logfilename);
}
void SAPT::oetrans()
{
    strans();
    vtrans();
}
void SAPT::strans()
{
  double **Sij, **SiB;

  Sij = block_matrix(calc_info_.nso,calc_info_.nso);
  SiB = block_matrix(calc_info_.nso,calc_info_.nmo);

  for (int i=0,ij=0; i<calc_info_.nso; i++) {
    for (int j=0; j<=i; j++,ij++) {
      Sij[i][j] = calc_info_.S[ij];
      Sij[j][i] = calc_info_.S[ij];
      }}

  free(calc_info_.S);

  C_DGEMM('N','N',calc_info_.nso,calc_info_.nmo,calc_info_.nso,1.0,
          &(Sij[0][0]),calc_info_.nso,&(calc_info_.CB[0][0]),calc_info_.nmo,
          0.0,&(SiB[0][0]),calc_info_.nmo);

  free_block(Sij);
  calc_info_.S_AB = block_matrix(calc_info_.nmo,calc_info_.nmo);

  C_DGEMM('T','N',calc_info_.nmo,calc_info_.nmo,calc_info_.nso,1.0,
          &(calc_info_.CA[0][0]),calc_info_.nmo,&(SiB[0][0]),calc_info_.nmo,
          0.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo);

  free_block(SiB);
}
void SAPT::vtrans()
{
  double **VAij, **VBij;
  double **VAiB, **VBAj;

  VAij = block_matrix(calc_info_.nso,calc_info_.nso);
  VAiB = block_matrix(calc_info_.nso,calc_info_.nmo);

  for (int i=0,ij=0; i<calc_info_.nso; i++) {
    for (int j=0; j<=i; j++,ij++) {
      VAij[i][j] = calc_info_.VA[ij];
      VAij[j][i] = calc_info_.VA[ij];
      }}

  free(calc_info_.VA);

  C_DGEMM('N','N',calc_info_.nso,calc_info_.nmo,calc_info_.nso,1.0,
          &(VAij[0][0]),calc_info_.nso,&(calc_info_.CB[0][0]),calc_info_.nmo,
          0.0,&(VAiB[0][0]),calc_info_.nmo);

  free_block(VAij);

  calc_info_.VABB = block_matrix(calc_info_.nmo,calc_info_.nmo);

  C_DGEMM('T','N',calc_info_.nmo,calc_info_.nmo,calc_info_.nso,1.0,
          &(calc_info_.CB[0][0]),calc_info_.nmo,&(VAiB[0][0]),calc_info_.nmo,
          0.0,&(calc_info_.VABB[0][0]),calc_info_.nmo);

  calc_info_.VAAB = block_matrix(calc_info_.nmo,calc_info_.nmo);

  C_DGEMM('T','N',calc_info_.nmo,calc_info_.nmo,calc_info_.nso,1.0,
          &(calc_info_.CA[0][0]),calc_info_.nmo,&(VAiB[0][0]),calc_info_.nmo,
          0.0,&(calc_info_.VAAB[0][0]),calc_info_.nmo);

  free_block(VAiB);

  VBij = block_matrix(calc_info_.nso,calc_info_.nso);
  VBAj = block_matrix(calc_info_.nmo,calc_info_.nso);

  for (int i=0,ij=0; i<calc_info_.nso; i++) {
    for (int j=0; j<=i; j++,ij++) {
      VBij[i][j] = calc_info_.VB[ij];
      VBij[j][i] = calc_info_.VB[ij];
      }}

  free(calc_info_.VB);

  C_DGEMM('T','N',calc_info_.nmo,calc_info_.nso,calc_info_.nso,1.0,
          &(calc_info_.CA[0][0]),calc_info_.nmo,&(VBij[0][0]),calc_info_.nso,
          0.0,&(VBAj[0][0]),calc_info_.nso);

  free_block(VBij);

  calc_info_.VBAA = block_matrix(calc_info_.nmo,calc_info_.nmo);

  C_DGEMM('N','N',calc_info_.nmo,calc_info_.nmo,calc_info_.nso,1.0,
          &(VBAj[0][0]),calc_info_.nso,&(calc_info_.CA[0][0]),calc_info_.nmo,
          0.0,&(calc_info_.VBAA[0][0]),calc_info_.nmo);

  calc_info_.VBAB = block_matrix(calc_info_.nmo,calc_info_.nmo);

  C_DGEMM('N','N',calc_info_.nmo,calc_info_.nmo,calc_info_.nso,1.0,
          &(VBAj[0][0]),calc_info_.nso,&(calc_info_.CB[0][0]),calc_info_.nmo,
          0.0,&(calc_info_.VBAB[0][0]),calc_info_.nmo);

  free_block(VBAj);
}
void SAPT::print_header()
{
 fprintf(outfile,"       S A P T   \n");
 fprintf(outfile,"    Ed Hohenstein\n") ;
 fprintf(outfile,"     6 June 2009\n") ;
 fprintf(outfile,"\n");
 fprintf(outfile,"    Orbital Information\n");
 fprintf(outfile,"  -----------------------\n");
 fprintf(outfile,"    NSO     = %9d\n",calc_info_.nso);
 fprintf(outfile,"    NMO     = %9d\n",calc_info_.nmo);
 fprintf(outfile,"    NRI     = %9d\n",calc_info_.nri);
 fprintf(outfile,"    NOCC_A  = %9d\n",calc_info_.noccA);
 fprintf(outfile,"    NOCC_B  = %9d\n",calc_info_.noccB);
 fprintf(outfile,"    NVIR_A  = %9d\n",calc_info_.nvirA);
 fprintf(outfile,"    NVIR_B  = %9d\n\n",calc_info_.nvirB);

 #ifdef _OPENMP
 fprintf(outfile,"Running SAPT with %d OMP threads\n\n",omp_get_max_threads());
 #endif

 fflush(outfile);

}
double** SAPT::get_DF_ints(int filenum, char *label, int length)
{
  double **A = block_matrix(length,calc_info_.nrio);
  psio_read_entry(filenum,label,(char *) A[0],
                  sizeof(double)*length*(ULI) calc_info_.nrio);
  return(A);
}
void SAPT::zero_disk(int file, char *array, char *zero, int nri, int ijmax)
{
  psio_address next_PSIF = PSIO_ZERO;

  for (int ij=0; ij<ijmax; ij++) {
    psio_write(file,array,zero,sizeof(double)*(ULI) nri,next_PSIF,&next_PSIF);
  }
}
double SAPT::CHF(int dfnum, char *OO, char *OV, char *VV, double **W, double **CHF, 
           double *evals, int nocc, int nvir)
{
  time_t start = time(NULL);
  time_t stop;
  int iter=0;
  double conv, tval, denom, E_old, E;
  double **C, **C_old;
  double **resC, **oldC;

  resC = block_matrix(nocc*nvir,params_.diisvec);
  oldC = block_matrix(nocc*nvir,params_.diisvec);

  C = block_matrix(nocc,nvir);
  C_old = block_matrix(nocc,nvir);

  int a,r,ar;

  #pragma omp for private(a,r,denom) schedule(static)
  for (a=0; a < nocc; a++) {
    for (r=nocc; r < calc_info_.nmo; r++) {
      denom = evals[a] - evals[r];
      C_old[a][r-nocc] = W[a][r-nocc]/denom;
      }}

  E_old = 2.0*C_DDOT(nocc*nvir,&(C_old[0][0]),1,&(W[0][0]),1);

  conv = 1.0;

  if (params_.print)
    fprintf(outfile,"Iter      Energy (mH)         dE (mH)            RMS (mH)    Time (s)\n");

  do {

    if (iter > params_.diisvec-1) {
      diis_update(C_old, resC, oldC, nocc, nvir);
      }

  A_mat(dfnum, OO, OV, VV, C_old, C, nocc, nvir, iter);

  #pragma omp for private(a,r,denom) schedule(static)
  for (a=0; a<nocc; a++) {
    for (r=0; r<nvir; r++) {
      denom = evals[a] - evals[r+nocc];
      C[a][r] += W[a][r];
      C[a][r] /= denom;
  }}

    E = 2.0*C_DDOT(nocc*nvir,&(C[0][0]),1,&(W[0][0]),1);

    conv = 0.0;

    for (int a=0, i=0; a < nocc; a++) {
    for (int r=nocc; r < calc_info_.nmo; r++,i++) {
      resC[i][iter % params_.diisvec] = fabs(C[a][r-nocc] - 
                          C_old[a][r-nocc]);
      conv += pow(C[a][r-nocc] - C_old[a][r-nocc],2);
      }}

    conv = sqrt(conv/(nocc*nvir));

    C_DCOPY(nocc*nvir,&(C[0][0]),1,&(C_old[0][0]),1);
    C_DCOPY(nocc*nvir,&(C[0][0]),1,&(oldC[0][iter % params_.diisvec]),
            params_.diisvec);

    iter++;
    stop = time(NULL);
    if (params_.print) {
      fprintf(outfile,"%4d %16.8lf %17.9lf %17.9lf    %10ld\n",iter,E*1000.0,(E_old-E)*1000.0,conv*1000.0,stop-start);
      fflush(outfile);
    }

    E_old = E;
    }
  while(conv > params_.d_conv && iter < params_.maxiter);

  if (conv <= params_.d_conv) {
    if (params_.print)
      fprintf(outfile,"\nCHF Iterations converged\n\n");
    }
  else {
    if (params_.print)
      fprintf(outfile,"\nCHF Iterations did not converge\n\n");
    }

  C_DCOPY(nocc*nvir,&(C_old[0][0]),1,&(CHF[0][0]),1);

  free_block(resC);
  free_block(oldC);
  free_block(C);
  free_block(C_old);

  fflush(outfile);

  return(E);
}

void SAPT::A_mat(int dfnum, char *OO, char *OV, char *VV, double **C_old, 
  double **C_new, int nocc, int nvir, int iter)
{
  time_t start;
  time_t stop;

  double **B_p_AR = get_DF_ints(dfnum,OV,nocc*nvir);
  double *C_p = init_array(calc_info_.nrio);

  C_DGEMV('t',nocc*nvir,calc_info_.nrio,1.0,&(B_p_AR[0][0]),calc_info_.nrio,
          &(C_old[0][0]),1,0.0,C_p,1);
  C_DGEMV('n',nocc*nvir,calc_info_.nrio,4.0,&(B_p_AR[0][0]),calc_info_.nrio,
          C_p,1,0.0,&(C_new[0][0]),1);

  free(C_p);

  double **D_p_AA = block_matrix(nocc*nocc,calc_info_.nrio);

  for (int a=0; a<nocc; a++) {
    C_DGEMM('N','N',nocc,calc_info_.nrio,nvir,1.0,&(C_old[0][0]),nvir,
      &(B_p_AR[a*nvir][0]),calc_info_.nrio,0.0,&(D_p_AA[a][0]),
      nocc*calc_info_.nrio);
  }

  for (int a=0; a<nocc; a++) {
    C_DGEMM('N','T',nocc,nvir,calc_info_.nrio,-1.0,&(D_p_AA[a*nocc][0]),
      calc_info_.nrio,&(B_p_AR[a*nvir][0]),calc_info_.nrio,1.0,&(C_new[0][0]),
      nvir);
  }

  free_block(D_p_AA);
  memset(&(B_p_AR[0][0]),'\0',sizeof(double)*nocc*nvir*calc_info_.nrio);

  double avail_mem = params_.memory;
  avail_mem -= 8.0*((double) (nocc*nvir)*(double) calc_info_.nrio);

  int temp_size = (int) ((avail_mem) / (8.0*((double) (nvir*calc_info_.nrio))));
  
  if (temp_size > nvir)
    temp_size = nvir;

  int blocks = (nvir)/temp_size;
  if ((nvir)%temp_size) blocks++;
  
  if (temp_size < 1) {
    fprintf(outfile,"Not enough memory in A Matrix formation\n\n");
    exit(0);
  } 

  if (params_.logfile && !iter) {
    fprintf(params_.logfilename,"     Block      Start       Stop\n");
    fprintf(params_.logfilename,"    -------  ----------  ----------\n");
    for (int t_r=0; t_r<blocks; t_r++) {
      int r_start = temp_size*t_r;
      int r_stop = temp_size*(t_r+1);
      if (r_stop > nvir)
        r_stop = nvir;
      fprintf(params_.logfilename,"       %3d   %10d  %10d\n",t_r,r_start,r_stop);
    }
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }
  if (params_.logfile) {
    fprintf(params_.logfilename,"    Iteration %d\n\n",iter+1);
    fflush(params_.logfilename);
  }

  double **B_p_RR = block_matrix(temp_size*nvir,calc_info_.nrio);

  psio_address next_PSIF = PSIO_ZERO;
  for (int t_r=0; t_r<blocks; t_r++) {
    int r_start = temp_size*t_r;
    int r_stop = temp_size*(t_r+1);
    if (r_stop > nvir)
      r_stop = nvir;

    if (params_.logfile) {
      fprintf(params_.logfilename,"      Starting Block %3d ... ",t_r);
      start = time(NULL);
      fflush(params_.logfilename);
    }

    psio_read(dfnum,VV,(char *) &(B_p_RR[0][0]),sizeof(double)*
      (r_stop-r_start)*nvir*(ULI) calc_info_.nrio,next_PSIF,&next_PSIF);
    for (int r=r_start; r<r_stop; r++) {
      C_DGEMM('N','N',nocc,calc_info_.nrio,nvir,1.0,&(C_old[0][0]),nvir,
        &(B_p_RR[(r-r_start)*nvir][0]),calc_info_.nrio,1.0,&(B_p_AR[r][0]),
        nvir*calc_info_.nrio);
    }

    if (params_.logfile) {
      stop = time(NULL);
      fprintf(params_.logfilename,"finished %14ld seconds\n",stop-start);
      fflush(params_.logfilename);
    }
  }
  free_block(B_p_RR);

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }

  double **B_p_AA = get_DF_ints(dfnum,OO,nocc*nocc);

  for (int a=0; a<nocc; a++) {
    C_DGEMM('N','T',nocc,nvir,calc_info_.nrio,-1.0,&(B_p_AA[a*nocc][0]),
      calc_info_.nrio,&(B_p_AR[a*nvir][0]),calc_info_.nrio,1.0,&(C_new[0][0]),
      nvir);
  }

  free_block(B_p_AA);
  free_block(B_p_AR);
}

void SAPT::diis_update(double **C, double **resC, double **oldC, int nocc, int nvir)
{
  int *ipiv;
  double *Cvec;
  double **Bmat;

  ipiv = init_int_array(params_.diisvec+1);
  Bmat = block_matrix(params_.diisvec+1,params_.diisvec+1);
  Cvec = (double *) malloc((params_.diisvec+1)*sizeof(double));

  for (int i=0; i<params_.diisvec; i++) {
    for (int j=0; j<=i; j++) {
      Bmat[i][j] = Bmat[j][i] = C_DDOT(nocc*nvir,&(resC[0][i]),params_.diisvec,
                          &(resC[0][j]),params_.diisvec);
      }}

  for (int i=0; i<params_.diisvec; i++) {
    Bmat[params_.diisvec][i] = -1.0;
    Bmat[i][params_.diisvec] = -1.0;
    Cvec[i] = 0.0;
    }

  Bmat[params_.diisvec][params_.diisvec] = 0.0;
  Cvec[params_.diisvec] = -1.0;

  C_DGESV(params_.diisvec+1,1,&(Bmat[0][0]),params_.diisvec+1,
          &(ipiv[0]),&(Cvec[0]),params_.diisvec+1);

  for (int i=0,m=0; i<nocc; i++) {
    for (int j=0; j<nvir; j++,m++) {
      C[i][j] = C_DDOT(params_.diisvec,&(Cvec[0]),1,&(oldC[m][0]),1);
        }}

  free(ipiv);
  free(Cvec);
  free_block(Bmat);
}

double** SAPT::W_ints(int dfnum, char *OV, double *diagOO, double **Vints, int nocc, 
int nvir)
{
  double **W = block_matrix(nocc,nvir);

  for(int a=0; a<nocc; a++){
    C_DAXPY(nvir,1.0,&(Vints[a][nocc]),1,&(W[a][0]),1);
  }

  double **B_p_AR = get_DF_ints(dfnum,OV,nocc*nvir);
    
  C_DGEMV('n',nocc*nvir,calc_info_.nrio-3,2.0,&(B_p_AR[0][0]),calc_info_.nrio,
          diagOO,1,1.0,&(W[0][0]),1);

  free_block(B_p_AR);

  return(W);
}

}}
