/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/

#include <stdexcept>
#include <cstdlib>
#include <psifiles.h>
#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include "hfwfn.h"

using namespace std;
using namespace psi::dboc;

HFWavefunction::HFWavefunction()
{
  int unit_opened = 1;
  if (!psio_open_check(PSIF_CHKPT)) {
    chkpt_init(PSIO_OPEN_OLD);
    unit_opened = 0;
  }

  num_mo_ = chkpt_rd_nmo();
  num_so_ = chkpt_rd_nso();
  num_ao_ = chkpt_rd_nao();
  
  refnum_ = (reftype) chkpt_rd_ref();

  nirreps_ = chkpt_rd_nirreps();
  clsdpi_ = chkpt_rd_clsdpi();
  openpi_ = chkpt_rd_openpi();
  orbspi_ = chkpt_rd_orbspi();

  ndocc_ = nsocc_ = 0;
  for(int i=0; i<nirreps_; i++) {
    ndocc_ += clsdpi_[i];
    nsocc_ += openpi_[i];
  }
  
  aotoso_ = chkpt_rd_usotao();
  rref_ = chkpt_rd_rref();
  if (refnum_ == ref_uhf || refnum_ == ref_uks) {
    alpha_evec_ = chkpt_rd_alpha_scf();
    beta_evec_ = chkpt_rd_beta_scf();
  }
  else {
    alpha_evec_ = chkpt_rd_scf();
    beta_evec_ = NULL;
  }

  if (!unit_opened)
    chkpt_close();
}

HFWavefunction::~HFWavefunction()
{
  free(clsdpi_);
  free(openpi_);
  free(orbspi_);
  free_block(aotoso_);
  free_block(rref_);
  free_block(alpha_evec_);
  if (beta_evec_ != NULL) free_block(beta_evec_);
}

int
HFWavefunction::num_mo() { return num_mo_; }

int
HFWavefunction::num_so() { return num_so_; }

int
HFWavefunction::ndocc() { return ndocc_; }

int
HFWavefunction::nsocc() { return nsocc_; }

int
HFWavefunction::nalpha() { return ndocc_ + nsocc_; }

int
HFWavefunction::nbeta() { return ndocc_; }

int
HFWavefunction::num_ao() { return num_ao_; }

int
HFWavefunction::nirreps() { return nirreps_; }

int*
HFWavefunction::clsdpi() { return clsdpi_; }

int*
HFWavefunction::openpi() { return openpi_; }

int*
HFWavefunction::orbspi() { return orbspi_; }

double**
HFWavefunction::alpha_evec() { return alpha_evec_; }

double**
HFWavefunction::beta_evec() {
  if (refnum_ == ref_uhf || refnum_ == ref_uks)
    return beta_evec_;
  else
    return alpha_evec_;
}

double**
HFWavefunction::aotoso() { return aotoso_; }

double**
HFWavefunction::rref() { return rref_; }

void
HFWavefunction::set_rref(double** rref)
{
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      rref_[i][j] = rref[i][j];
}
