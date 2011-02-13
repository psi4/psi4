#include "cc.h"

#include <libchkpt/chkpt.hpp>
#include <libmints/mints.h>

using namespace std;
using namespace psi;

namespace psi { namespace dfcc {

CC::CC(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
  : Wavefunction(options, psio, chkpt) 
{
  get_params();
  get_ribasis();
}

CC::~CC()
{
}

void CC::get_options()
{
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
  double* evals_t = chkpt_->rd_evals();
  double** C_t = chkpt_->rd_scf();

  evals_ = shared_ptr<Vector>(new Vector("Epsilon (full)",nmo_));
  evalsp_ = evals_->pointer();
  memcpy(static_cast<void*> (evalsp_), static_cast<void*> (evals_t), nmo_*sizeof(double));   
  C_ = shared_ptr<Matrix>(new Matrix("C (full)", nso_, nmo_));
  Cp_ = C_->pointer();
  memcpy(static_cast<void*> (Cp_[0]), static_cast<void*> (C_t[0]), nmo_*nso_*sizeof(double));   

  // Convenience matrices (may make it easier on the helper objects)
  // ...because Rob is a pussy
  // asshole
  evals_aocc_ = shared_ptr<Vector>(new Vector("Epsilon (Active Occupied)",naocc_));
  evals_aoccp_ = evals_aocc_->pointer();
  evals_avir_ = shared_ptr<Vector>(new Vector("Epsilon (Active Virtual)",navir_));
  evals_avirp_ = evals_avir_->pointer();

  C_aocc_ = shared_ptr<Matrix>(new Matrix("C (Active Occupied)", nso_, naocc_));
  C_aoccp_ = C_aocc_->pointer();
  C_avir_ = shared_ptr<Matrix>(new Matrix("C (Active Virtual)", nso_, navir_));
  C_avirp_ = C_avir_->pointer();

  memcpy(static_cast<void*> (evals_aoccp_), static_cast<void*> (&evals_t[nfocc_]), naocc_*sizeof(double));   
  memcpy(static_cast<void*> (evals_avirp_), static_cast<void*> (&evals_t[nocc_]), navir_*sizeof(double));  

  for (int m = 0; m < nso_; m++) { 
    memcpy(static_cast<void*> (C_aoccp_[m]), static_cast<void*> (&C_t[m][nfocc_]), naocc_*sizeof(double));   
    memcpy(static_cast<void*> (C_avirp_[m]), static_cast<void*> (&C_t[m][nocc_]), naocc_*sizeof(double));   
  }

  free(evals_t);
  free_block(C_t); 
}

void CC::get_ribasis()
{
  shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(
    options_.get_str("BASIS_PATH")));
  ribasis_ = BasisSet::construct(parser, molecule_, options_.get_str(
    "RI_BASIS_CC"));
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

}}
