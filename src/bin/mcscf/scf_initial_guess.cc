/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

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
    outfile->Printf("\n  Reading MOs from the checkpoint file.");
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
