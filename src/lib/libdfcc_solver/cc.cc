#include "cc.h"

#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.hpp>
#include <libmints/mints.h>
#include <lib3index/3index.h>
#include <psifiles.h>

using namespace boost;
using namespace psi;

namespace psi { namespace dfcc {

CC::CC(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
  : Wavefunction(options, psio, chkpt)
{
  get_options();
  get_params();
  get_ribasis();
}

CC::~CC()
{
}

void CC::get_options()
{
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");
}

void CC::get_params()
{
  // Init with checkpoint until Rob gets bitchy
  // RP: I'm bitchy. But not a bitch unless you want chinese food.
  // MO basis info
  nirrep_ = chkpt_->rd_nirreps();

  if (nirrep_ != 1)
    throw PsiException("You want symmetry? Try ccenergy", __FILE__,
      __LINE__);

  int *clsdpi_ = new int[8];
  int *orbspi_ = new int[8];
  int *frzcpi_ = new int[8];
  int *frzvpi_ = new int[8];

  clsdpi_ = chkpt_->rd_clsdpi();
  orbspi_ = chkpt_->rd_orbspi();
  frzcpi_ = chkpt_->rd_frzcpi();
  frzvpi_ = chkpt_->rd_frzvpi();

  nso_ = chkpt_->rd_nso();
  nmo_ = orbspi_[0];
  nocc_ = clsdpi_[0];
  nvir_ = orbspi_[0]-clsdpi_[0];
  nfocc_ = frzcpi_[0];
  nfvir_ = frzvpi_[0];
  naocc_ = clsdpi_[0]-frzcpi_[0];
  navir_ = orbspi_[0]-clsdpi_[0]-frzvpi_[0];
  namo_ = navir_ + naocc_;

  delete[] clsdpi_;
  delete[] orbspi_;
  delete[] frzcpi_;
  delete[] frzvpi_;

  // Reference wavefunction info
  Eref_ = chkpt_->rd_escf();
  energies_["Reference Energy"] = Eref_;
  double* evals_t = chkpt_->rd_evals();
  double** C_t = chkpt_->rd_scf();

  evals_ = shared_ptr<Vector>(new Vector("Epsilon (full)",nmo_));
  evalsp_ = evals_->pointer();
  memcpy(static_cast<void*> (evalsp_), static_cast<void*> (evals_t), nmo_*sizeof(double));
  C_ = SharedMatrix(new Matrix("C (full)", nso_, nmo_));
  Cp_ = C_->pointer();
  memcpy(static_cast<void*> (Cp_[0]), static_cast<void*> (C_t[0]), nmo_*nso_*sizeof(double));

  // Convenience matrices (may make it easier on the helper objects)
  // ...because Rob is a pussy
  // asshole
  evals_aocc_ = shared_ptr<Vector>(new Vector("Epsilon (Active Occupied)",naocc_));
  evals_aoccp_ = evals_aocc_->pointer();
  evals_avir_ = shared_ptr<Vector>(new Vector("Epsilon (Active Virtual)",navir_));
  evals_avirp_ = evals_avir_->pointer();

  C_aocc_ = SharedMatrix(new Matrix("C (Active Occupied)", nso_, naocc_));
  C_aoccp_ = C_aocc_->pointer();
  C_avir_ = SharedMatrix(new Matrix("C (Active Virtual)", nso_, navir_));
  C_avirp_ = C_avir_->pointer();

  memcpy(static_cast<void*> (evals_aoccp_), static_cast<void*> (&evals_t[nfocc_]), naocc_*sizeof(double));
  memcpy(static_cast<void*> (evals_avirp_), static_cast<void*> (&evals_t[nocc_]), navir_*sizeof(double));

  for (int m = 0; m < nso_; m++) {
    memcpy(static_cast<void*> (C_aoccp_[m]), static_cast<void*> (&C_t[m][nfocc_]), naocc_*sizeof(double));
    memcpy(static_cast<void*> (C_avirp_[m]), static_cast<void*> (&C_t[m][nocc_]), navir_*sizeof(double));
  }

  free(evals_t);
  free_block(C_t);

  sss_ = options_.get_double("SCALE_SS");
  oss_ = options_.get_double("SCALE_OS");
  // We're almost always using Laplace, but, just in case
  denominator_algorithm_ = options_.get_str("DENOMINATOR_ALGORITHM");
  denominator_delta_ = options_.get_double("DENOMINATOR_DELTA");

  schwarz_cutoff_ = options_.get_double("INTS_TOLERANCE");
  fitting_condition_ = options_.get_double("FITTING_COND");
  fitting_algorithm_ = options_.get_str("FITTING_TYPE");

  doubles_ = memory_ / 8L;
}

void CC::get_ribasis()
{
  shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
  // If the user doesn't spec a basis name, pick it yourself
  // TODO: Verify that the basis assign does not messs this up
  if (options_.get_str("DF_BASIS_CC") == "") {
      basisset_->molecule()->set_basis_all_atoms(options_.get_str("BASIS") + "-RI", "DF_BASIS_CC");
      fprintf(outfile, "  No auxiliary basis selected, defaulting to %s-RI\n\n", options_.get_str("BASIS").c_str());
  }
  ribasis_ = BasisSet::construct(parser, molecule_, "DF_BASIS_CC");
  zero_ = BasisSet::zero_ao_basis_set();
  ndf_ = ribasis_->nbf();
}

void CC::print_header()
{
     fprintf(outfile,"    Orbital Information\n");
     fprintf(outfile,"  -----------------------\n");
     fprintf(outfile,"    NSO      = %8d\n",nso_);
     fprintf(outfile,"    NMO Tot  = %8d\n",nmo_);
     fprintf(outfile,"    NMO Act  = %8d\n",namo_);
     fprintf(outfile,"    NOCC Tot = %8d\n",nocc_);
     fprintf(outfile,"    NOCC Frz = %8d\n",nfocc_);
     fprintf(outfile,"    NOCC Act = %8d\n",naocc_);
     fprintf(outfile,"    NVIR Tot = %8d\n",nvir_);
     fprintf(outfile,"    NVIR Frz = %8d\n",nfvir_);
     fprintf(outfile,"    NVIR Act = %8d\n",navir_);
     fprintf(outfile,"    NDF      = %8d\n",ndf_);
     fprintf(outfile,"\n");
     fflush(outfile);
}

void CC::df_integrals()
{
  IntegralFactory rifactory_J(ribasis_, zero_, ribasis_, zero_);
  TwoBodyAOInt* Jint = rifactory_J.eri();

  double **J = block_matrix(ndf_,ndf_);
  double **J_mhalf = block_matrix(ndf_,ndf_);
  const double *Jbuffer = Jint->buffer();

  for (int MU=0; MU < ribasis_->nshell(); ++MU) {
    int nummu = ribasis_->shell(MU).nfunction();

    for (int NU=0; NU <= MU; ++NU) {
      int numnu = ribasis_->shell(NU).nfunction();

      Jint->compute_shell(MU, 0, NU, 0);

      int index = 0;
      for (int mu=0; mu < nummu; ++mu) {
        int omu = ribasis_->shell(MU).function_index() + mu;

        for (int nu=0; nu < numnu; ++nu, ++index) {
          int onu = ribasis_->shell(NU).function_index() + nu;

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
    int numPshell = ribasis_->shell(Pshell).nfunction();
    if (numPshell > maxPshell) maxPshell = numPshell;
  }

  double** AO_RI = block_matrix(maxPshell,nso_*nso_);
  double* halftrans = init_array(nso_*namo_);
  double** MO_RI = block_matrix(ndf_,namo_*namo_);

  for (int Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
    int numPshell = ribasis_->shell(Pshell).nfunction();
    for (int MU=0; MU < basisset_->nshell(); ++MU) {
      int nummu = basisset_->shell(MU).nfunction();
      for (int NU=0; NU <= MU; ++NU) {
        int numnu = basisset_->shell(NU).nfunction();

        eri->compute_shell(Pshell, 0, MU, NU);

        for (int P=0, index=0; P < numPshell; ++P) {

          for (int mu=0; mu < nummu; ++mu) {
            int omu = basisset_->shell(MU).function_index() + mu;

            for (int nu=0; nu < numnu; ++nu, ++index) {
              int onu = basisset_->shell(NU).function_index() + nu;

              AO_RI[P][omu*nso_+onu] = buffer[index];
              AO_RI[P][onu*nso_+omu] = buffer[index];
            }
          }
        }
      }
    }
    for (int P=0; P < numPshell; ++P) {
      int oP = ribasis_->shell(Pshell).function_index() + P;
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

void CC::mo_integrals()
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

  int virtri = navir_*(navir_+1)/2;
  int svirtri = navir_*(navir_-1)/2;
  double **VV = block_matrix(navir_,navir_);
  double *xVV = init_array(virtri);
  double *yVV = init_array(svirtri);

  zero_disk(DFCC_INT_FILE,"VVVV+ Integrals",(char *) &(xVV[0]),virtri,
    virtri);
  zero_disk(DFCC_INT_FILE,"VVVV- Integrals",(char *) &(yVV[0]),svirtri,
    svirtri);

  psio_address next_VVVVp = PSIO_ZERO;
  psio_address next_VVVVm = PSIO_ZERO;

  for (int a=0; a < navir_; a++) {
  for (int b=0; b <= a; b++) {

    C_DGEMM('N','T',navir_,navir_,ndf_,1.0,&(B_p_VV[a*navir_][0]),ndf_,
      &(B_p_VV[b*navir_][0]),ndf_,0.0,&(VV[0][0]),navir_);

    for (int c=0; c < navir_; c++) {
    for (int d=0; d <= c; d++) {
      int cd = INDEX(c,d);
      xVV[cd] = VV[c][d] + VV[d][c];
    }}
    psio_->write(DFCC_INT_FILE,"VVVV+ Integrals",(char *) &(xVV[0]),
      virtri*(ULI) sizeof(double),next_VVVVp,&next_VVVVp);

    if (a != b) {
      for (int c=0; c < navir_; c++) {
      for (int d=0; d < c; d++) {
        int cd = INDEX(c-1,d);
        yVV[cd] = VV[c][d] - VV[d][c];
      }}
      psio_->write(DFCC_INT_FILE,"VVVV- Integrals",(char *) &(yVV[0]),
        svirtri*(ULI) sizeof(double),next_VVVVm,&next_VVVVm);
    }

  }}

  free(xVV);
  free(yVV);
  free_block(VV);
  free_block(B_p_VV);
}

void CC::iajb_ibja(double *ijkl)
{
  double *X = init_array(naocc_*naocc_*navir_*navir_);

  for (int i=0; i<naocc_; i++) {
  for (int a=0; a<navir_; a++) {
    for (int j=0; j<naocc_; j++) {
    for (int b=0; b<navir_; b++) {
      int iajb = i*navir_*naocc_*navir_ + a*naocc_*navir_ + j*navir_ + b;
      int ibja = i*navir_*naocc_*navir_ + b*naocc_*navir_ + j*navir_ + a;
      X[ibja] = ijkl[iajb];
  }}}}

  C_DCOPY(naocc_*naocc_*navir_*navir_,X,1,ijkl,1);

  free(X);
/*
  double *X = init_array(navir_);

  for (int i=0; i<naocc_; i++) {
  for (int a=0; a<navir_; a++) {
    for (int j=0; j<naocc_; j++) {
    int ia = i*navir_ + a;
    int ja = j*navir_ + a;
    long int iajb = ia*naocc_*(long int) navir_ + (long int) j*navir_;
    long int ibja = i*navir_*naocc_*(long int) navir_ + (long int) ja;
    C_DCOPY(navir_,&(ijkl[iajb]),1,X,1);
    C_DCOPY(navir_,&(ijkl[ibja]),naocc_*navir_,&(ijkl[iajb]),1);
    C_DCOPY(navir_,X,1,&(ijkl[ibja]),naocc_*navir_);
  }}}

  free(X);
*/
}

void CC::iajb_ijab(double *ijkl)
{
  double *X = init_array(naocc_*naocc_*navir_*navir_);

  for (int i=0; i<naocc_; i++) {
  for (int a=0; a<navir_; a++) {
    for (int j=0; j<naocc_; j++) {
    for (int b=0; b<navir_; b++) {
      int iajb = i*navir_*naocc_*navir_ + a*naocc_*navir_ + j*navir_ + b;
      int ijab = i*naocc_*navir_*navir_ + j*navir_*navir_ + a*navir_ + b;
      X[ijab] = ijkl[iajb];
  }}}}

  C_DCOPY(naocc_*naocc_*navir_*navir_,X,1,ijkl,1);

  free(X);
/*
  double **X = block_matrix(navir_,naocc_);

  for (int i=0; i<naocc_; i++) {
  for (int b=0; b<navir_; b++) {
    long int iajb = i*navir_*naocc_*(long int) navir_ + (long int) b;
    C_DCOPY(naocc_*navir_,&(ijkl[iajb]),navir_,X[0],1);
    for (int j=0; j<naocc_; j++) {
      int ij = i*naocc_ + j;
      long int ijab = ij*navir_*(long int) navir_ + (long int) b;
      C_DCOPY(navir_,&(X[0][j]),naocc_,&(ijkl[ijab]),navir_);
    }
  }}

  free_block(X);
*/
}

void CC::ijab_iajb(double *ijkl)
{
  double *X = init_array(naocc_*naocc_*navir_*navir_);

  for (int i=0; i<naocc_; i++) {
  for (int a=0; a<navir_; a++) {
    for (int j=0; j<naocc_; j++) {
    for (int b=0; b<navir_; b++) {
      int iajb = i*navir_*naocc_*navir_ + a*naocc_*navir_ + j*navir_ + b;
      int ijab = i*naocc_*navir_*navir_ + j*navir_*navir_ + a*navir_ + b;
      X[iajb] = ijkl[ijab];
  }}}}

  C_DCOPY(naocc_*naocc_*navir_*navir_,X,1,ijkl,1);

  free(X);
/*
  double **X = block_matrix(naocc_,navir_);

  for (int i=0,ij=0; i<naocc_; i++) {
  for (int j=0,jb=0; j<naocc_; j++,ij++) {
    long int ijab = ij*navir_*(long int) navir_;
    C_DCOPY(naocc_*navir_,&(ijkl[ijab]),navir_,X[0],1);
    for (int b=0; b<navir_; b++,jb++) {
      long int iajb = i*navir_*naocc_*(long int) navir_ + (long int) jb;
      C_DCOPY(naocc_,&(X[0][b]),navir_,&(ijkl[iajb]),navir_);
    }
  }}

  free_block(X);
*/
}

void CC::zero_disk(int file, const char *array, char *zero, int nri, int ijmax)
{
  psio_address next_PSIF = PSIO_ZERO;

  for (int ij=0; ij<ijmax; ij++) {
    psio_->write(file,array,zero,sizeof(double)*(ULI) nri,next_PSIF,&next_PSIF);
  }
}

DFCCDIIS::DFCCDIIS(int diisfile, int length, int maxvec, shared_ptr<PSIO> psio)
  : psio_(psio)
{
    diis_file_ = diisfile;
    psio_->open(diis_file_,0);

    max_diis_vecs_ = maxvec;

    vec_length_ = length;

    curr_vec_ = 0;
    num_vecs_ = 0;
}

DFCCDIIS::~DFCCDIIS()
{
    psio_->close(diis_file_,0);
}

void DFCCDIIS::store_vectors(double *t_vec, double *err_vec)
{
    char *diis_vec_label = get_vec_label(curr_vec_);
    char *diis_err_label = get_err_label(curr_vec_);
    curr_vec_ = (curr_vec_+1)%max_diis_vecs_;
    num_vecs_++;
    if (num_vecs_ > max_diis_vecs_) num_vecs_ = max_diis_vecs_;

    psio_->write_entry(diis_file_,diis_vec_label,(char *) &(t_vec[0]),
      vec_length_*(ULI) sizeof(double));

    psio_->write_entry(diis_file_,diis_err_label,(char *) &(err_vec[0]),
      vec_length_*(ULI) sizeof(double));

    free(diis_vec_label);
    free(diis_err_label);
}

void DFCCDIIS::store_current_vector(char *t_vec)
{
    char *diis_vec_label = get_vec_label(curr_vec_);

    psio_->write_entry(diis_file_,diis_vec_label,t_vec,
        vec_length_*(ULI) sizeof(double));

    free(diis_vec_label);
}

void DFCCDIIS::store_error_vector(char *err_vec)
{
    char *diis_err_label = get_err_label(curr_vec_);

    psio_->write_entry(diis_file_,diis_err_label,err_vec,
        vec_length_*sizeof(double));

    free(diis_err_label);
}

void DFCCDIIS::increment_vectors()
{
    curr_vec_ = (curr_vec_+1)%max_diis_vecs_;
    num_vecs_++;
    if (num_vecs_ > max_diis_vecs_) num_vecs_ = max_diis_vecs_;
}

char *DFCCDIIS::get_last_vec_label()
{
    int vec_num = (curr_vec_+max_diis_vecs_-1)%max_diis_vecs_;
    return(get_vec_label(vec_num));
}

void DFCCDIIS::get_new_vector(double *vec_j, double *vec_i)
{
    int *ipiv;
    double *Cvec;
    double **Bmat;

    ipiv = init_int_array(num_vecs_+1);
    Bmat = block_matrix(num_vecs_+1,num_vecs_+1);
    Cvec = (double *) malloc((num_vecs_+1)*sizeof(double));

    for (int i=0; i<num_vecs_; i++) {
      char *err_label_i = get_err_label(i);
      psio_->read_entry(diis_file_,err_label_i,(char *) &(vec_i[0]),
        vec_length_*(ULI) sizeof(double));
      for (int j=0; j<=i; j++) {
        char *err_label_j = get_err_label(j);
        psio_->read_entry(diis_file_,err_label_j,(char *) &(vec_j[0]),
          vec_length_*(ULI) sizeof(double));
        Bmat[i][j] = Bmat[j][i] = C_DDOT(vec_length_,vec_i,1,vec_j,1);
        free(err_label_j);
      }
      free(err_label_i);
    }

    for (int i=0; i<num_vecs_; i++) {
      Bmat[num_vecs_][i] = -1.0;
      Bmat[i][num_vecs_] = -1.0;
      Cvec[i] = 0.0;
    }

    Bmat[num_vecs_][num_vecs_] = 0.0;
    Cvec[num_vecs_] = -1.0;

    C_DGESV(num_vecs_+1,1,&(Bmat[0][0]),num_vecs_+1,&(ipiv[0]),&(Cvec[0]),
      num_vecs_+1);

    memset(vec_j,'\0',sizeof(double)*vec_length_);

    for (int i=0; i<num_vecs_; i++) {
      char *vec_label_i = get_vec_label(i);
      psio_->read_entry(diis_file_,vec_label_i,(char *) &(vec_i[0]),
        vec_length_*(ULI) sizeof(double));
      C_DAXPY(vec_length_,Cvec[i],vec_i,1,vec_j,1);
      free(vec_label_i);
    }

    free(ipiv);
    free(Cvec);
    free_block(Bmat);
}

void DFCCDIIS::get_new_vector(double **vec_i, int cols)
{
    int *ipiv;
    double *Cvec;
    double **Bmat;

    double *temp = init_array(cols);
    int rows = vec_length_/cols;

    ipiv = init_int_array(num_vecs_+1);
    Bmat = block_matrix(num_vecs_+1,num_vecs_+1);
    Cvec = (double *) malloc((num_vecs_+1)*sizeof(double));

    for (int i=0; i<num_vecs_; i++) {
      char *err_label_i = get_err_label(i);
      psio_->read_entry(diis_file_,err_label_i,(char *) &(vec_i[0][0]),
        vec_length_*(ULI) sizeof(double));
      for (int j=0; j<=i; j++) {
        char *err_label_j = get_err_label(j);
        psio_address next_j = PSIO_ZERO;
        double Bval = 0.0;
        for (int n=0; n < rows; n++) {
            psio_->read(diis_file_,err_label_j,(char *) &(temp[0]),
                cols*(ULI) sizeof(double),next_j,&next_j);
            Bval = C_DDOT(cols,vec_i[n],1,temp,1);
        }
        Bmat[i][j] = Bmat[j][i] = Bval;
        free(err_label_j);
      }
      free(err_label_i);
    }

    for (int i=0; i<num_vecs_; i++) {
      Bmat[num_vecs_][i] = -1.0;
      Bmat[i][num_vecs_] = -1.0;
      Cvec[i] = 0.0;
    }

    Bmat[num_vecs_][num_vecs_] = 0.0;
    Cvec[num_vecs_] = -1.0;

    C_DGESV(num_vecs_+1,1,&(Bmat[0][0]),num_vecs_+1,&(ipiv[0]),&(Cvec[0]),
      num_vecs_+1);

    memset(vec_i[0],'\0',sizeof(double)*vec_length_);

    for (int i=0; i<num_vecs_; i++) {
        char *vec_label_i = get_vec_label(i);
        psio_address next_i = PSIO_ZERO;
        for (int n=0; n < rows; n++) {
            psio_->read(diis_file_,vec_label_i,(char *) &(temp[0]),
                cols*(ULI) sizeof(double),next_i,&next_i);
            C_DAXPY(cols,Cvec[i],temp,1,vec_i[n],1);
        }
        free(vec_label_i);
    }

    free(ipiv);
    free(Cvec);
    free(temp);
    free_block(Bmat);
}

char *DFCCDIIS::get_err_label(int num)
{
    char *label = (char *) malloc(16*sizeof(char));
    sprintf(label,"Error vector %2d",num);
    return(label);
}

char *DFCCDIIS::get_vec_label(int num)
{
    char *label = (char *) malloc(10*sizeof(char));
    sprintf(label,"Vector %2d",num);
    return(label);
}

}}
