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
  psio_->open(DFCC_INT_FILE,PSIO_OPEN_NEW);
  df_integrals();
}

CCD::~CCD()
{
  psio_->close(DFCC_INT_FILE,1);
}

double CCD::compute_energy()
{
  return(0.0);
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

}}
