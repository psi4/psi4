/*! \file
    \ingroup BASIS
    \brief Enter brief description of file here 
*/

/*! \defgroup BASIS libbasis: The Basis Set Library */

#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <psifiles.h>
#include "basisset.h"

namespace psi {

BasisSet::BasisSet(int chkptfile)
{
  num_shells_ = chkpt_rd_nshell();
  num_prims_ = chkpt_rd_nprim();
  num_ao_ = chkpt_rd_nao();
  // Psi 3 only allows either all Cartesians or all Spherical harmonics only
  num_bf_ = chkpt_rd_nso();
  puream_ = chkpt_rd_puream();
  max_am_ = chkpt_rd_max_am();

  shells_ = new GaussianShell*[num_shells_];

/* shells_ = */ init_shells();
  /* shell_pairs = init_shell_pairs(); */
}

BasisSet::BasisSet(const BasisSet& S) :
  num_prims_(S.num_prims_), num_shells_(S.num_shells_), num_ao_(S.num_ao_),
  num_bf_(S.num_bf_), max_am_(S.max_am_), puream_(S.puream_)
{
  shells_ = new GaussianShell*[num_shells_];
  for(int s=0; s<num_shells_; s++)
    shells_[s] = new GaussianShell(*(S.shells_[s]));
    
  ncenters_ = S.ncenters_;
  coords_ = block_matrix(ncenters_,3);
  for(int c=0; c<ncenters_; c++)
    for(int xyz=0; xyz<3; xyz++)
      coords_[c][xyz] = S.coords_[c][xyz];
  
  shell_fbf_ = new int[num_shells_];
  shell_fao_ = new int[num_shells_];
  shell_center_ = new int[num_shells_];
  for(int s=0; s<num_shells_; s++) {
    shell_fbf_[s] = S.shell_fbf_[s];
    shell_fao_[s] = S.shell_fao_[s];
    shell_center_[s] = S.shell_center_[s];
  }
}

BasisSet::~BasisSet()
{
  //dealloc_pairs();
  for(int s=0; s<num_shells_; s++)
    shells_[s]->~GaussianShell();
  delete[] shells_;
  free(shell_center_);
  free(shell_fao_);
  free(shell_fbf_);
  free_block(coords_);
}

void BasisSet::init_shells()
{
   /*--- retrieve angular momentum of each shell (1=s, 2=p, 3=d, etc  ) ---*/
   int *shell_am = chkpt_rd_stype();

   /*--- retrieve number of primitives per shell ---*/
   int *shell_num_prims = chkpt_rd_snumg();

   /*--- retrieve exponents of primitive gaussians ---*/
   PSI_FLOAT *exponents = chkpt_rd_exps();

   /*--- retrieve coefficients of primitive gaussians ---*/
   PSI_FLOAT **ccoeffs = chkpt_rd_contr_full();

   /*--- retrieve pointer to first primitive in shell ---*/
   int *shell_fprim = chkpt_rd_sprim();

   /*--- retrieve pointer to first basisfn in shell ---*/
   shell_fbf_ = chkpt_rd_sloc_new();

   /*--- retrieve pointer to first AO in shell ---*/
   shell_fao_ = chkpt_rd_sloc();

   /*--- retrieve location of shells (which atom it's centered on) ---*/
   shell_center_ = chkpt_rd_snuc();

   /*--- retrieve number of centers ---*/
   ncenters_ = chkpt_rd_natom();

   /*--- retrieve geometry ---*/
   coords_ = chkpt_rd_geom();
   
   // Only segmented contractions can be handled in Psi 3 at present
   int ncontr = 1;
   int* am = new int[ncontr];
   for (int i=0; i<num_shells_; i++) {
     am[0] = shell_am[i]-1;
     int fprim = shell_fprim[i] - 1;
     int nprims = shell_num_prims[i];
     int center = shell_center_[i] - 1;
     PSI_FLOAT **cc = new PSI_FLOAT*[nprims];
     for(int p=0; p<nprims; p++) {
       cc[p] = new PSI_FLOAT[ncontr];
       cc[p][0] = ccoeffs[fprim+p][am[0]];
     }
     shells_[i] = new GaussianShell(nprims, ncontr, am, puream_, &(exponents[fprim]), cc, coords_[center]);
     for(int p=0; p<nprims; p++) {
       delete[] cc[p];
     }
     delete[] cc;
   }

   free_block(ccoeffs);
   free(exponents);
   free(shell_am);
   free(shell_num_prims);
   free(shell_fprim);
}

void BasisSet::check_shell_index(int si) const {
  if (si < 0 || si >= num_shells_)
    throw std::runtime_error("ERROR: BasisSet::check_shell_index -- shell index out of bounds");
}

int BasisSet::num_prims() const { return num_prims_; };
int BasisSet::num_shells() const { return num_shells_; };
int BasisSet::num_ao() const { return num_ao_; };
int BasisSet::num_bf() const { return num_bf_; };
int BasisSet::max_am() const { return max_am_; };

GaussianShell& BasisSet::shell(int si) const
{
  check_shell_index(si);
  return *(shells_[si]);
};

int BasisSet::first_bf(int si) const
{
  check_shell_index(si);
  return shell_fbf_[si] - 1;
};

int BasisSet::first_ao(int si) const
{
  check_shell_index(si);
  return shell_fao_[si] - 1;
};

int BasisSet::center(int si) const
{
  check_shell_index(si);
  return shell_center_[si] - 1;
}

void BasisSet::set_center(int ci, PSI_FLOAT O[3])
{
  for(int xyz=0; xyz<3; xyz++)
    coords_[ci][xyz] = O[xyz];
  for(int si=0; si<num_shells_; si++) {
    if (shell_center_[si] == ci + 1)
      shells_[si]->set_origin(O);
  }
}

PSI_FLOAT BasisSet::get_center(int ci, int i)
{
  return coords_[ci][i];
}

void BasisSet::print(char *id, FILE* outfile) const {
  char indent1[] = "  ";
  char indent2[] = "    ";

  fprintf(outfile, "%s-Basis Set %s\n",indent1,id);
  fprintf(outfile, "%sNumber of shells              = %d\n",indent2,num_shells_);
  fprintf(outfile, "%sNumber of basis functions     = %d\n",indent2,num_bf_);
  fprintf(outfile, "%sNumber of Cartesian Gaussians = %d\n",indent2,num_ao_);
  fprintf(outfile, "%sSpherical Harmonics?          = %s\n",indent2,(puream_ ? "true" : "false"));
  fprintf(outfile, "%sMax angular momentum          = %d\n\n",indent2,max_am_);

  fprintf(outfile, "%sShells:\n\n",indent2);
  for(int s=0; s<num_shells_; s++)
    shells_[s]->print(s,outfile);
}

}
