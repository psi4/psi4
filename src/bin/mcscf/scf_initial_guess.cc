#include <iostream>
#include <cstdlib>
#include <cstdio>

#include <libchkpt/chkpt.hpp>
#include <liboptions/liboptions.h>

#include "scf.h"

namespace psi{
  extern Chkpt* chkpt_;
}

namespace psi{ namespace mcscf{

void SCF::initial_guess()
{
  using namespace psi;

  bool read_MOs = false;
  double** saved_MOs = chkpt_->rd_scf();
  if(saved_MOs != NULL){
    free(saved_MOs);
    if(options_.get_bool("MO_READ"))
      read_MOs = true;
  }
  if(read_MOs){
    for(int h = 0; h < nirreps; ++h){
      if(sopi[h] > 0){
        double** block = chkpt_->rd_scf_irrep(h);
        for(int i = 0; i < sopi[h]; ++i){
          for(int j = 0; j < sopi[h]; ++j){
            C->set(h,i,j,block[i][j]);
          }
        }
        free(block);
      }
    }
    fprintf(outfile,"\n  Reading MOs from the checkpoint file.");
  }else{
    SBlockMatrix H_t("H_t",nirreps,sopi,sopi);
    SBlockVector eigenvectors("H_t",nirreps,sopi);

    transform(H,H_t,S_sqrt_inv);

    H_t.diagonalize(C_t,eigenvectors);

    C.multiply(false,false,S_sqrt_inv,C_t);

    epsilon = eigenvectors;

  }
  guess_occupation();
}

}} /* End Namespaces */
