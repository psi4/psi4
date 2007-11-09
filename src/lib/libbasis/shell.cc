/*! \file 
    \ingroup (BASIS)
    \brief Enter brief description of file here 
*/

#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>

#include "shell.h"


GaussianShell::GaussianShell(int nprims, int ncontr, int *am, bool puream, PSI_FLOAT *exps, PSI_FLOAT **ccoeffs, PSI_FLOAT origin[3]) :
  num_prims_(nprims), num_contr_(ncontr), puream_(puream)
{
  am_ = new int[num_contr_];
  for(int c=0; c<num_contr_; c++)
    am_[c] = am[c];
  min_am_ = am_[0];
  max_am_ = am_[0];
  for(int c=1; c<num_contr_; c++) {
    max_am_ = (max_am_ > am_[c] ? max_am_ : am_[c]);
    min_am_ = (min_am_ < am_[c] ? min_am_ : am_[c]);
  }

  num_ao_ = 0;
  for(int c=0; c<num_contr_; c++)
    num_ao_ += (am_[c]+1)*(am_[c]+2)/2;
  if (!puream_)
    num_bf_ = num_ao_;
  else {
    num_bf_ = 0;
    for(int c=0; c<num_contr_; c++)
      num_bf_ += 2*am_[c]+1;
  }
  
  exps_ = new PSI_FLOAT[num_prims_];
  ccoeffs_ = new PSI_FLOAT*[num_prims_];
  for(int p=0; p<num_prims_; p++) {
    ccoeffs_[p] = new PSI_FLOAT[num_contr_];
    exps_[p] = exps[p];
    for(int c=0; c<num_contr_; c++) {
      ccoeffs_[p][c] = ccoeffs[p][c];
    }
  }

  O_[0] = origin[0];
  O_[1] = origin[1];
  O_[2] = origin[2];
}

GaussianShell::GaussianShell(const GaussianShell& S) :
  max_am_(S.max_am_), min_am_(S.min_am_), num_bf_(S.num_bf_), num_ao_(S.num_ao_),
  num_prims_(S.num_prims_), num_contr_(S.num_contr_), puream_(S.puream_)
{
  am_ = new int[num_contr_];
  for(int c=0; c<num_contr_; c++)
    am_[c] = S.am_[c];

  exps_ = new PSI_FLOAT[num_prims_];
  ccoeffs_ = new PSI_FLOAT*[num_prims_];
  for(int p=0; p<num_prims_; p++) {
    ccoeffs_[p] = new PSI_FLOAT[num_contr_];
    exps_[p] = S.exps_[p];
    for(int c=0; c<num_contr_; c++) {
      ccoeffs_[p][c] = S.ccoeffs_[p][c];
    }
  }

  O_[0] = S.O_[0];
  O_[1] = S.O_[1];
  O_[2] = S.O_[2];
}

GaussianShell::~GaussianShell()
{
  delete[] am_;
  delete[] exps_;
  for(int p=0; p<num_prims_; p++)
    delete[] ccoeffs_[p];
  delete[] ccoeffs_;
}

void GaussianShell::check_contraction_index_(int ci) const
{
  if (ci < 0 || ci > num_contr_)
    throw std::runtime_error("ERROR: GaussianShell::check_contraction_index -- index out of bounds");
}

void GaussianShell::check_primitive_index_(int pi) const
{
  if (pi < 0 || pi > num_prims_)
    throw std::runtime_error("ERROR: GaussianShell::check_primitive_index -- index out of bounds");
}

int GaussianShell::am() const
{
  if (max_am_ !=  min_am_)
    throw std::runtime_error("ERROR: GaussianShell::am -- shell contains contractions of several angular momenta");
  return am_[0];
}

int GaussianShell::am(int ci) const
{
  check_contraction_index_(ci);
  return am_[ci];
}

int GaussianShell::max_am() const
{
  return max_am_;
}

int GaussianShell::min_am() const
{
  return min_am_;
}

PSI_FLOAT GaussianShell::exp(int pi) const
{
  check_primitive_index_(pi);
  return exps_[pi];
}

PSI_FLOAT GaussianShell::cc(int ci, int pi) const
{
  check_primitive_index_(pi);
  check_contraction_index_(ci);
  return ccoeffs_[pi][ci];
}

void GaussianShell::set_origin(PSI_FLOAT O[3])
{
  O_[0] = O[0];
  O_[1] = O[1];
  O_[2] = O[2];
}

void GaussianShell::print(int id, FILE* outfile) const {
  char indent1[] = "    ";
  char indent2[] = "      ";

  fprintf(outfile, "%sGaussian Shell %d\n",indent1,id);
  fprintf(outfile, "%sNumber of contractions        = %d\n",indent2,num_contr_);
  fprintf(outfile, "%sNumber of primitives          = %d\n",indent2,num_prims_);
  fprintf(outfile, "%sNumber of basis functions     = %d\n",indent2,num_bf_);
  fprintf(outfile, "%sNumber of Cartesian Gaussians = %d\n",indent2,num_ao_);
  fprintf(outfile, "%sSpherical Harmonics?          = %s\n",indent2,(puream_ ? "true" : "false"));
  if (max_am_ == min_am_)
    fprintf(outfile, "%sAngular momentum              = %d\n",indent2,max_am_);
  else {
    fprintf(outfile, "%sMax angular momentum          = %d\n",indent2,max_am_);
    fprintf(outfile, "%sMin angular momentum          = %d\n",indent2,min_am_);
  }

  fprintf(outfile, "%sExponent",indent2);
  for(int c=0; c<num_contr_; c++)
    fprintf(outfile, " Contr. %3d",c);
  fprintf(outfile,"\n");
  for(int p=0; p<num_prims_; p++) {
    fprintf(outfile, "%s%15.10lf",indent2,exps_[p]);
    for(int c=0; c<num_contr_; c++)
      fprintf(outfile, " %12.9lf",ccoeffs_[p][c]);
    fprintf(outfile,"\n");
  }
  fprintf(outfile,"\n");
}
