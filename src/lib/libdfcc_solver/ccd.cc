#include "ccd.h"

#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libmints/mints.h>

using namespace std;
using namespace psi;

namespace psi { namespace dfcc {

CCD::CCD(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
  : CC(options, psio, chkpt)
{
  print_header();
  psio_->open(DFCC_INT_FILE,PSIO_OPEN_NEW);
  df_integrals();
  mo_integrals();
}

CCD::~CCD()
{
  psio_->close(DFCC_INT_FILE,1);
}

double CCD::compute_energy()
{
  tIAJB_ = init_array(naocc_*naocc_*navir_*(ULI) navir_);
  t2IAJB_ = init_array(naocc_*naocc_*navir_*(ULI) navir_);
  vIAJB_ = init_array(naocc_*naocc_*navir_*(ULI) navir_);
  xIAJB_ = init_array(naocc_*naocc_*navir_*(ULI) navir_);

  double emp2;

  term_1();
  apply_denom();
  emp2 = energy();

  fprintf(outfile,"  Reference Energy            %18.10lf\n",Eref_);
  fprintf(outfile,"  Correlation Energy          %18.10lf\n",emp2);
  fprintf(outfile,"  Total DF-MP2 Energy         %18.10lf\n\n",Eref_+emp2);
  fprintf(outfile,"  Iter       Energy (H)          dE (H)             RMS (H)\n");
  fflush(outfile);

  int iter = 1;
  double e_new;
  double e_old = emp2;

  do {

    term_1(); // t_ab^ij <- (ia|jb)
    term_2(); // t_ab^ij <- t_ac^ik [2(jb|kc)-(jk|bc)]
    term_3(); // t_ab^ij <- 2 t_ac^ik [2(kc|ld)-(kd|lc)] t_bd^jl
    term_4(); // t_ab^ij <- - t_ac^ij [2(kc|ld)-(kd|lc)] t_bd^kl
    term_5(); // t_ab^ij <- - t_ab^ik [2(kc|ld)-(kd|lc)] t_cd^jl

    iajb_ibja(tIAJB_);

    term_6(); // t_ab^ij <- - t_ac^ik [2(kc|ld)-(kd|lc)] t_bd^lj
    term_7(); // t_ab^ij <- - t_ac^ki (jb|kc)
              // t_ab^ij <- t_ac^ki (kc|ld) t_bd^lj

    iajb_ibja(vIAJB_);
    iajb_ibja(t2IAJB_);

    term_8(); // t_ab^ij <- t_cb^il (kc|ld) t_da^jk
    term_9(); // t_ab^ij <- - t_bc^ki (jk|ac)

    iajb_ibja(tIAJB_);
    iajb_ibja(t2IAJB_);

    iajb_ijab(tIAJB_);
    iajb_ijab(t2IAJB_);

    term_10(); // t_ab^ij <- t_ab^kl (ik|jl)
    term_11(); // t_ab^ij <- t_cd^ij (ac|bd)
    term_12(); // t_ab^ij <- t_cd^ij (kc|ld) t_ab^kl

    iajb_ijab(tIAJB_);
    iajb_ijab(t2IAJB_);

    symmetrize(); // t_ab^ij <- 0.5*( t_ab^ij + t_ba^ji )
    apply_denom();
    e_new = energy();

    fprintf(outfile,"  %4d %16.8lf %17.9lf %17.9lf\n",iter,e_new,e_old-e_new,0.0);
    fflush(outfile);

    e_old = e_new;

    iter++;
  }
  while(iter <= options_.get_int("MAXITER"));

  fprintf(outfile,"\n");
  fprintf(outfile,"  Reference Energy            %18.10lf\n",Eref_);
  fprintf(outfile,"  Correlation Energy          %18.10lf\n",e_old);
  fprintf(outfile,"  Total DF-CCD Energy         %18.10lf\n\n",Eref_+e_old);

  free(tIAJB_);
  free(t2IAJB_);
  free(vIAJB_);
  free(xIAJB_);

  return(0.0);
}

void CCD::print_header()
{
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*                       DF-CCD                         *\n");
    fprintf(outfile, "\t*               Coupled-Cluster Doubles                *\n");
    fprintf(outfile, "\t*                with all sorts of shit                *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*            Rob Parrish and Ed Hohenstein             *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\n");
    CC::print_header();

}

void CCD::df_integrals()
{
  IntegralFactory rifactory_J(ribasis_, zero_, ribasis_, zero_);
  TwoBodyAOInt* Jint = rifactory_J.eri();

  double **J = block_matrix(ndf_,ndf_);
  double **J_mhalf = block_matrix(ndf_,ndf_);
  const double *Jbuffer = Jint->buffer();

  for (int MU=0; MU < ribasis_->nshell(); ++MU) {
    int nummu = ribasis_->shell(MU)->nfunction();
    
    for (int NU=0; NU <= MU; ++NU) {
      int numnu = ribasis_->shell(NU)->nfunction();

      Jint->compute_shell(MU, 0, NU, 0);

      int index = 0;
      for (int mu=0; mu < nummu; ++mu) {
        int omu = ribasis_->shell(MU)->function_index() + mu;

        for (int nu=0; nu < numnu; ++nu, ++index) {
          int onu = ribasis_->shell(NU)->function_index() + nu;

          J[omu][onu] = Jbuffer[index];
        }
      }
    }
  }

  double* eigval = init_array(ndf_);
  int lwork = ndf_ * 3;
  double* work = init_array(lwork);
  int stat = C_DSYEV('v','u',ndf_,J[0],ndf_,eigval,work,lwork);
  if (stat != 0) {
    fprintf(outfile, "C_DSYEV failed\n");
    exit(PSI_RETURN_FAILURE);
  }
  free(work);

  double **J_copy = block_matrix(ndf_,ndf_);
  C_DCOPY(ndf_*ndf_,J[0],1,J_copy[0],1);

  for (int i=0; i<ndf_; i++) {
    if (eigval[i] < 1.0E-10)
      eigval[i] = 0.0;
    else {
      eigval[i] = 1.0 / sqrt(eigval[i]);
    }
    C_DSCAL(ndf_, eigval[i], J[i], 1);
  }
  free(eigval);

  C_DGEMM('t','n',ndf_,ndf_,ndf_,1.0,J_copy[0],ndf_,J[0],ndf_,0.0,J_mhalf[0],
    ndf_);

  free_block(J);
  free_block(J_copy);

  IntegralFactory rifactory(ribasis_, zero_, basisset_, basisset_);
  TwoBodyAOInt* eri = rifactory.eri();
  const double *buffer = eri->buffer();

  int maxPshell = 0;
  for (int Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
    int numPshell = ribasis_->shell(Pshell)->nfunction();
    if (numPshell > maxPshell) maxPshell = numPshell;
  }

  double** AO_RI = block_matrix(maxPshell,nso_*nso_);
  double* halftrans = init_array(nso_*namo_);
  double** MO_RI = block_matrix(ndf_,namo_*namo_);

  for (int Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
    int numPshell = ribasis_->shell(Pshell)->nfunction();
    for (int MU=0; MU < basisset_->nshell(); ++MU) {
      int nummu = basisset_->shell(MU)->nfunction();
      for (int NU=0; NU <= MU; ++NU) {
        int numnu = basisset_->shell(NU)->nfunction();

        eri->compute_shell(Pshell, 0, MU, NU);

        for (int P=0, index=0; P < numPshell; ++P) {

          for (int mu=0; mu < nummu; ++mu) {
            int omu = basisset_->shell(MU)->function_index() + mu;

            for (int nu=0; nu < numnu; ++nu, ++index) {
              int onu = basisset_->shell(NU)->function_index() + nu;

              AO_RI[P][omu*nso_+onu] = buffer[index];
              AO_RI[P][onu*nso_+omu] = buffer[index];
            }
          }
        }
      }
    }
    for (int P=0; P < numPshell; ++P) {
      int oP = ribasis_->shell(Pshell)->function_index() + P;
      C_DGEMM('T', 'N', namo_, nso_, nso_, 1.0, &(Cp_[0][nfocc_]), nmo_,
        AO_RI[P], nso_, 0.0, halftrans, nso_);
      C_DGEMM('N', 'N', namo_, namo_, nso_, 1.0, halftrans, nso_,
        &(Cp_[0][nfocc_]), nmo_, 0.0, MO_RI[oP], namo_);
    }
  }

  free(halftrans);
  free_block(AO_RI);

  double** MO_RI_J = block_matrix(namo_*namo_,ndf_);

  C_DGEMM('T','T',namo_*namo_,ndf_,ndf_,1.0,MO_RI[0],namo_*namo_,
    J_mhalf[0],ndf_,0.0,MO_RI_J[0],ndf_);

  free_block(J_mhalf);
  free_block(MO_RI);

  psio_address next_DF_OO = PSIO_ZERO;
  psio_address next_DF_OV = PSIO_ZERO;
  psio_address next_DF_VV = PSIO_ZERO;

  for (int i=0; i<naocc_; i++) {
    psio_->write(DFCC_INT_FILE,"OO DF Integrals",(char *)
      &(MO_RI_J[i*namo_][0]),naocc_*ndf_*sizeof(double),next_DF_OO,
      &next_DF_OO);
  }

  for (int i=0; i<naocc_; i++) {
    psio_->write(DFCC_INT_FILE,"OV DF Integrals",(char *)
      &(MO_RI_J[i*namo_+naocc_][0]),navir_*ndf_*sizeof(double),next_DF_OV,
      &next_DF_OV);
  }

  for (int i=naocc_; i<namo_; i++) {
    psio_->write(DFCC_INT_FILE,"VV DF Integrals",(char *)
      &(MO_RI_J[i*namo_+naocc_][0]),navir_*ndf_*sizeof(double),next_DF_VV,
      &next_DF_VV);
  }

  free_block(MO_RI_J);
}

void CCD::mo_integrals()
{
  double **B_p_OV = block_matrix(naocc_*navir_,ndf_);
  double **vOVOV = block_matrix(naocc_*navir_,naocc_*navir_);

  psio_->read_entry(DFCC_INT_FILE,"OV DF Integrals",(char *)
      &(B_p_OV[0][0]),naocc_*navir_*ndf_*sizeof(double));

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,ndf_,1.0,B_p_OV[0],ndf_,
    B_p_OV[0],ndf_,0.0,vOVOV[0],naocc_*navir_);

  psio_->write_entry(DFCC_INT_FILE,"OVOV Integrals",(char *)
    &(vOVOV[0][0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  free_block(B_p_OV);

  double **gOVOV = block_matrix(naocc_*navir_,naocc_*navir_);

  for (int i=0,ia=0; i<naocc_; i++) {
  for (int a=0; a<navir_; a++,ia++) {
    for (int j=0,jb=0; j<naocc_; j++) {
    for (int b=0; b<navir_; b++,jb++) {
      int ib = i*navir_+b;
      int ja = j*navir_+a;
      gOVOV[ia][jb] = 2.0*vOVOV[ia][jb] - vOVOV[ib][ja];
  }}}}

  psio_->write_entry(DFCC_INT_FILE,"G OVOV Integrals",(char *)
    &(gOVOV[0][0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  free_block(vOVOV);
  free_block(gOVOV);

  double **B_p_OO = block_matrix(naocc_*naocc_,ndf_);

  psio_->read_entry(DFCC_INT_FILE,"OO DF Integrals",(char *)
      &(B_p_OO[0][0]),naocc_*naocc_*ndf_*sizeof(double));

  double **vOOOO = block_matrix(naocc_*naocc_,naocc_*naocc_);

  for (int i=0,ij=0; i<naocc_; i++) {
    for (int j=0; j<naocc_; j++,ij++) {
      C_DGEMM('N','T',naocc_,naocc_,ndf_,1.0,B_p_OO[i*naocc_],ndf_,
        B_p_OO[j*naocc_],ndf_,0.0,vOOOO[ij],naocc_);
  }}

  psio_->write_entry(DFCC_INT_FILE,"OOOO Integrals",(char *)
    &(vOOOO[0][0]),naocc_*naocc_*naocc_*naocc_*sizeof(double));

  free_block(vOOOO);

  double **B_p_VV = block_matrix(navir_*navir_,ndf_);

  psio_->read_entry(DFCC_INT_FILE,"VV DF Integrals",(char *)
      &(B_p_VV[0][0]),navir_*navir_*ndf_*sizeof(double));

  double **vOVVO = block_matrix(naocc_*navir_,naocc_*navir_);

  for (int i=0,ia=0; i<naocc_; i++) {
    for (int a=0; a<navir_; a++,ia++) {
      C_DGEMM('N','T',naocc_,navir_,ndf_,1.0,B_p_OO[i*naocc_],ndf_,
        B_p_VV[a*navir_],ndf_,0.0,vOVVO[ia],navir_);
  }}

  psio_->write_entry(DFCC_INT_FILE,"OOVV Integrals",(char *)
    &(vOVVO[0][0]),naocc_*navir_*naocc_*navir_*sizeof(double));
 
  double **gOVVO = block_matrix(naocc_*navir_,naocc_*navir_);

  psio_->read_entry(DFCC_INT_FILE,"OVOV Integrals",(char *)
    &(gOVVO[0][0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  C_DSCAL(naocc_*navir_*naocc_*navir_,2.0,gOVVO[0],1);
  C_DAXPY(naocc_*navir_*naocc_*navir_,-1.0,vOVVO[0],1,gOVVO[0],1);

  psio_->write_entry(DFCC_INT_FILE,"G OVVO Integrals",(char *)
    &(gOVVO[0][0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  free_block(B_p_OO);
  free_block(vOVVO);
  free_block(gOVVO);

  double **vVVVV = block_matrix(navir_*navir_,navir_*navir_);

  for (int a=0,ab=0; a<navir_; a++) {
    for (int b=0; b<navir_; b++,ab++) {
      C_DGEMM('N','T',navir_,navir_,ndf_,1.0,B_p_VV[a*navir_],ndf_,
        B_p_VV[b*navir_],ndf_,0.0,vVVVV[ab],navir_);
  }}
    
  psio_->write_entry(DFCC_INT_FILE,"VVVV Integrals",(char *)
    &(vVVVV[0][0]),navir_*navir_*navir_*navir_*sizeof(double));

  free_block(B_p_VV);
  free_block(vVVVV);
}

double CCD::energy()
{
  double energy;

  psio_->read_entry(DFCC_INT_FILE,"G OVOV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  energy = C_DDOT(naocc_*navir_*naocc_*(long int) navir_,t2IAJB_,1,vIAJB_,1);

  C_DCOPY(naocc_*navir_*naocc_*(long int) navir_,t2IAJB_,1,tIAJB_,1);

  return(energy);
}

void CCD::apply_denom()
{
  for (int i=0,ia=0; i<naocc_; i++) {
  for (int a=0; a<navir_; a++,ia++) {
    for (int j=0,jb=0; j<naocc_; j++) {
    for (int b=0; b<navir_; b++,jb++) {
      long int iajb = ia*naocc_*(long int) navir_ + (long int) jb;
      t2IAJB_[iajb] /= evals_aoccp_[i] + evals_aoccp_[j] - evals_avirp_[a]
        - evals_avirp_[b];
    }}
  }}
} 

void CCD::symmetrize()
{
  for (int ia=0; ia<naocc_*navir_; ia++) {
    for (int jb=0; jb<ia; jb++) {
      long int iajb = ia*naocc_*(long int) navir_ + (long int) jb;
      long int jbia = jb*naocc_*(long int) navir_ + (long int) ia;
      double tval = t2IAJB_[iajb] + t2IAJB_[jbia];
      t2IAJB_[iajb] = 0.5*tval;
      t2IAJB_[jbia] = 0.5*tval;
  }}
}

void CCD::term_1()
{
  psio_->read_entry(DFCC_INT_FILE,"OVOV Integrals",(char *)
    &(t2IAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));
}

void CCD::term_2()
{
  psio_->read_entry(DFCC_INT_FILE,"G OVVO Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,2.0,
    tIAJB_,naocc_*navir_,vIAJB_,naocc_*navir_,1.0,t2IAJB_,naocc_*navir_);
}

void CCD::term_3()
{
  psio_->read_entry(DFCC_INT_FILE,"G OVOV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,1.0,
    tIAJB_,naocc_*navir_,vIAJB_,naocc_*navir_,0.0,xIAJB_,naocc_*navir_);

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,2.0,
    xIAJB_,naocc_*navir_,tIAJB_,naocc_*navir_,1.0,t2IAJB_,naocc_*navir_);
}

void CCD::term_4()
{
  double **xAB = block_matrix(navir_,navir_);

  C_DGEMM('T','N',navir_,navir_,naocc_*navir_*naocc_,1.0,vIAJB_,navir_,
    tIAJB_,navir_,0.0,xAB[0],navir_);

  C_DGEMM('N','N',naocc_*navir_*naocc_,navir_,navir_,-2.0,tIAJB_,navir_,
    xAB[0],navir_,1.0,t2IAJB_,navir_);

  free_block(xAB);
}

void CCD::term_5()
{
  double **xIJ = block_matrix(naocc_,naocc_);

  C_DGEMM('N','T',naocc_,naocc_,navir_*naocc_*navir_,1.0,vIAJB_,
    navir_*naocc_*navir_,tIAJB_,navir_*naocc_*navir_,0.0,xIJ[0],naocc_);

  C_DGEMM('T','N',naocc_,navir_*naocc_*navir_,naocc_,-2.0,xIJ[0],naocc_,
    tIAJB_,navir_*naocc_*navir_,1.0,t2IAJB_,navir_*naocc_*navir_);

  free_block(xIJ);
}

void CCD::term_6()
{
  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,-2.0,
    xIAJB_,naocc_*navir_,tIAJB_,naocc_*navir_,1.0,t2IAJB_,naocc_*navir_);
}

void CCD::term_7()
{
  psio_->read_entry(DFCC_INT_FILE,"OVOV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,1.0,
    tIAJB_,naocc_*navir_,vIAJB_,naocc_*navir_,0.0,xIAJB_,naocc_*navir_);

  C_DAXPY(naocc_*navir_*naocc_*(long int) navir_,-2.0,xIAJB_,1,t2IAJB_,1);

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,1.0,
    xIAJB_,naocc_*navir_,tIAJB_,naocc_*navir_,1.0,t2IAJB_,naocc_*navir_);
}

void CCD::term_8()
{
  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,1.0,
    tIAJB_,naocc_*navir_,vIAJB_,naocc_*navir_,0.0,xIAJB_,naocc_*navir_);

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,1.0,
    xIAJB_,naocc_*navir_,tIAJB_,naocc_*navir_,1.0,t2IAJB_,naocc_*navir_);
}

void CCD::term_9()
{
  psio_->read_entry(DFCC_INT_FILE,"OOVV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,-2.0,
    tIAJB_,naocc_*navir_,vIAJB_,naocc_*navir_,1.0,t2IAJB_,naocc_*navir_);
}

void CCD::term_10()
{
  double *vIJKL = init_array(naocc_*naocc_*naocc_*(long int) naocc_);

  psio_->read_entry(DFCC_INT_FILE,"OOOO Integrals",(char *)
    &(vIJKL[0]),naocc_*naocc_*naocc_*naocc_*sizeof(double));

  C_DGEMM('N','N',naocc_*naocc_,navir_*navir_,naocc_*naocc_,1.0,vIJKL,
    naocc_*naocc_,tIAJB_,navir_*navir_,1.0,t2IAJB_,navir_*navir_);

  free(vIJKL);
}

void CCD::term_11()
{
  double *vABCD = init_array(navir_*navir_*navir_*(long int) navir_);

  psio_->read_entry(DFCC_INT_FILE,"VVVV Integrals",(char *)
    &(vABCD[0]),navir_*navir_*navir_*navir_*sizeof(double));

  C_DGEMM('N','T',naocc_*naocc_,navir_*navir_,navir_*navir_,1.0,tIAJB_,
    navir_*navir_,vABCD,navir_*navir_,1.0,t2IAJB_,navir_*navir_);

  free(vABCD);
}

void CCD::term_12()
{
  psio_->read_entry(DFCC_INT_FILE,"OVOV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  iajb_ijab(vIAJB_);

  double *vIJKL = init_array(naocc_*naocc_*naocc_*(long int) naocc_);

  C_DGEMM('N','T',naocc_*naocc_,naocc_*naocc_,navir_*navir_,1.0,tIAJB_,
    navir_*navir_,vIAJB_,navir_*navir_,0.0,vIJKL,naocc_*naocc_);

  C_DGEMM('N','N',naocc_*naocc_,navir_*navir_,naocc_*naocc_,1.0,vIJKL,
    naocc_*naocc_,tIAJB_,navir_*navir_,1.0,t2IAJB_,navir_*navir_);

  free(vIJKL);

}

}}
