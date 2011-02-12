#include "cc.h"

#include <libchkpt/chkpt.hpp>
#include <libmints/basisset.h>
#include <libmints/basisset_parser.h>

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

  delete[] clsdpi_;
  delete[] orbspi_;
  delete[] frzcpi_;
  delete[] frzvpi_;

  // Reference wavefunction info
  Eref_ = chkpt_->rd_escf();
  evals_ = chkpt_->rd_evals();
  C_ = chkpt_->rd_scf();
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

}}
