#include "sapt0.h"

using namespace boost;

namespace psi { namespace sapt {

void SAPT0::disp20()
{
  shared_ptr<Vector> evals_aoccA(new Vector(aoccA_));
  shared_ptr<Vector> evals_virA(new Vector(nvirA_));
  shared_ptr<Vector> evals_aoccB(new Vector(aoccB_));
  shared_ptr<Vector> evals_virB(new Vector(nvirB_));

  for (int a=0; a<aoccA_; a++)
    evals_aoccA->set(0,a,evalsA_[a+foccA_]);
  for (int r=0; r<nvirA_; r++)
    evals_virA->set(0,r,evalsA_[r+noccA_]);
  for (int b=0; b<aoccB_; b++)
    evals_aoccB->set(0,b,evalsB_[b+foccB_]);
  for (int s=0; s<nvirB_; s++)
    evals_virB->set(0,s,evalsB_[s+noccB_]);

  denom_ = shared_ptr<SAPTLaplaceDenominator>(new SAPTLaplaceDenominator(
    evals_aoccA,evals_virA,evals_aoccB,evals_virB,
    options_.get_double("DENOMINATOR_DELTA"),debug_));

  denom_->debug();
}

}}
