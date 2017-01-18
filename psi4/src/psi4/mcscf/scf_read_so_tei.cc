/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#include <iostream>
#include <cmath>

#include "psi4/psifiles.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "scf.h"

#define MAX(i,j) ((i>j) ? i : j)
#define MIN(i,j) ((i>j) ? j : i)
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

using namespace std;

extern FILE* outfile;

namespace psi{ namespace mcscf{

extern MemoryManager* memory_manager;

void SCF::read_so_tei()
{
  generate_pairs();

  total_symmetric_block_size = INDEX(pairpi[0]-1,pairpi[0]-1)+1;

  size_t free_memory = memory_manager->get_FreeMemory();

  // Determine the number of matrix elements of the PK (and K) matrix to hold in core
  if(reference == rhf){
    nin_core        = min(free_memory / sizeof(double),total_symmetric_block_size);
  }else{
    nin_core        = min(free_memory / (2 * sizeof(double)),total_symmetric_block_size);
  }
  if(nin_core != total_symmetric_block_size)
    out_of_core = true;

  size_t total_symmetric_pairs = pairpi[0];

  // Determine how to split the two-electron operators
  nbatch            = 0;
  size_t pq_incore  = 0;
  size_t pqrs_index = 0;

  batch_pq_min[0]   = 0;
  batch_pq_max[0]   = 0;
  batch_index_min[0]= 0;
  batch_index_max[0]= 0;

  for(size_t pq = 0; pq < total_symmetric_pairs; ++pq){
    // Increment counters
    if(pq_incore + pq + 1 > nin_core){
      // The batch is full. Save info.
      batch_pq_max[nbatch]         = pq;
      batch_pq_min[nbatch + 1]     = pq;
      batch_index_max[nbatch]      = pqrs_index;
      batch_index_min[nbatch + 1]  = pqrs_index;
      pq_incore = 0;
      nbatch++;
    }
    pq_incore  += pq + 1;
    pqrs_index += pq + 1;
  }
  if(batch_pq_max[nbatch] != total_symmetric_pairs){
    batch_pq_max[nbatch]     = total_symmetric_pairs;
    batch_index_max[nbatch]  = total_symmetric_block_size;
    nbatch++;
  }

  for(int batch = 0; batch < nbatch; ++batch){
    batch_size[batch] = batch_index_max[ batch] - batch_index_min[batch];
    outfile->Printf("\n  batch %3d pq = [%8ld,%8ld] index = [%16ld,%16ld]",
        batch,
        batch_pq_min[batch],batch_pq_max[batch],
        batch_index_min[batch],batch_index_max[batch]);
  }


  // Allocate the PK matrix
  allocate1(double,PK,nin_core);
  for(size_t i=0; i < nin_core; i++)
    PK[i]    =0.0;
  outfile->Printf("\n\n  Allocated the PK matrix (%ld elements) ",(long int)nin_core);


  if(reference != rhf){
    // Allocate the K matrix
    allocate1(double,K,nin_core);
    for(size_t i=0; i < nin_core; i++)
      K[i]    =0.0;
    outfile->Printf("\n  Allocated the  K matrix (%ld elements) ",(long int)nin_core);

  }

  if(reference == rhf)
    read_so_tei_form_PK();
  else
    read_so_tei_form_PK_and_K();
}

void SCF::read_so_tei_form_PK()
{
  outfile->Printf("\n  Reading the two-electron integrals to form PK ... ");


  for(int batch = 0; batch < nbatch; ++batch){
    outfile->Printf("\n  batch %3d ... ",batch);

    // Compute the minimum and maximum indices
    size_t min_index   = batch_index_min[batch];
    size_t max_index   = batch_index_max[batch];
    size_t buffer_size = max_index - min_index;

    for(size_t pqrs = 0; pqrs < buffer_size; ++pqrs)
      PK[pqrs] = 0.0;

    double value;
    size_t p,q,r,s,four_index;
    int ilsti,nbuf,fi,index;

    IWL ERIIN(psio_.get(), PSIF_SO_TEI, 0.0, 1, 1);
    ERIIN.set_keep_flag(1);
    do {
      ilsti = ERIIN.last_buffer();
      nbuf  = ERIIN.buffer_count();

      fi = 0;
      for (index = 0; index < nbuf; ++index) {
        p = (ERIIN.labels()[fi] >= 0 ? ERIIN.labels()[fi] : -ERIIN.labels()[fi]);
        q = ERIIN.labels()[fi+1];
        r = ERIIN.labels()[fi+2];
        s = ERIIN.labels()[fi+3];
        value = ERIIN.values()[index];

        if(pair_sym[p][q] == 0)
        {
          four_index = INDEX(pair[p][q],pair[r][s]);
          if((four_index >= min_index) && (four_index < max_index)){
            PK[four_index - min_index]    += value;
          }
        }
        if(pair_sym[p][r] == 0)
        {
          four_index = INDEX(pair[p][r],pair[q][s]);
          if((four_index >= min_index) && (four_index < max_index)){
            if((p == r) || (q == s))
              PK[four_index - min_index]   -= 0.5 * value;
            else
              PK[four_index - min_index]   -= 0.25 * value;
          }
        }
        if(pair_sym[p][s] == 0){
          four_index = INDEX(pair[p][s],pair[q][r]);
          if((four_index >= min_index) && (four_index < max_index)){
            if((p != q) && (r != s)){
              if((p == s) || (q == r))
                PK[four_index - min_index] -= 0.5 * value;
              else
                PK[four_index - min_index] -= 0.25 * value;
            }
          }
        }
        fi += 4;
      }
      if (!ilsti)
        ERIIN.fetch();
    } while (!ilsti);

    // Half the diagonal elements held in core
    for(size_t pq = batch_pq_min[batch]; pq < batch_pq_max[batch]; ++pq){
      PK[INDEX(pq,pq) - min_index] *= 0.5;
    }

    // Write the PK matrix to disk
    write_Raffanetti("PK",PK,batch);

    outfile->Printf("done.");

  }
  outfile->Printf("\n");

}

void SCF::read_so_tei_form_PK_and_K()
{
  outfile->Printf("\n  Reading the two-electron integrals to form PK and K ... ");


  for(int batch = 0; batch < nbatch; ++batch){
    outfile->Printf("\n  batch %3d ... ",batch);

    // Compute the minimum and maximum indices
    size_t min_index   = batch_index_min[batch];
    size_t max_index   = batch_index_max[batch];
    size_t buffer_size = max_index - min_index;

    for(size_t pqrs = 0; pqrs < buffer_size; ++pqrs){
      PK[pqrs] = 0.0; K[pqrs] = 0.0;
    }

    double value;
    size_t p,q,r,s,four_index;
    int ilsti,nbuf,fi,index;

    IWL ERIIN(psio_.get(), PSIF_SO_TEI, 0.0, 1, 1);
    ERIIN.set_keep_flag(1);

    do {
      ilsti = ERIIN.last_buffer();
      nbuf  = ERIIN.buffer_count();

      fi = 0;
      for (index = 0; index < nbuf; ++index) {
        p = (ERIIN.labels()[fi] >= 0 ? ERIIN.labels()[fi] : -ERIIN.labels()[fi]);
        q = ERIIN.labels()[fi+1];
        r = ERIIN.labels()[fi+2];
        s = ERIIN.labels()[fi+3];
        value = ERIIN.values()[index];

        if(pair_sym[p][q] == 0)
        {
          four_index = INDEX(pair[p][q],pair[r][s]);
          if((four_index >= min_index) && (four_index < max_index)){
            PK[four_index - min_index]    += value;
          }
        }
        if(pair_sym[p][r] == 0)
        {
          four_index = INDEX(pair[p][r],pair[q][s]);
          if((four_index >= min_index) && (four_index < max_index)){
            if((p == r) || (q == s)){
              PK[four_index - min_index]   -= 0.5 * value;
              K[four_index - min_index]   -= 0.5 * value;
            }else{
              PK[four_index - min_index]   -= 0.25 * value;
              K[four_index - min_index]   -= 0.25 * value;
            }
          }
        }
        if(pair_sym[p][s] == 0){
          four_index = INDEX(pair[p][s],pair[q][r]);
          if((four_index >= min_index) && (four_index < max_index)){
            if((p != q) && (r != s)){
              if((p == s) || (q == r)){
                PK[four_index - min_index] -= 0.5 * value;
                K[four_index - min_index] -= 0.5 * value;
              }else{
                PK[four_index - min_index] -= 0.25 * value;
                K[four_index - min_index] -= 0.25 * value;
              }
            }
          }
        }
        fi += 4;
      }
      if (!ilsti)
        ERIIN.fetch();
    } while (!ilsti);

    // Half the diagonal elements held in core
    for(size_t pq = batch_pq_min[batch]; pq < batch_pq_max[batch]; ++pq){
      PK[INDEX(pq,pq) - min_index] *= 0.5;
      K[INDEX(pq,pq) - min_index] *= 0.5;
    }

    // Write the PK matrix to disk
    write_Raffanetti("PK",PK,batch);
    write_Raffanetti("K",K,batch);

    outfile->Printf("done.");

  }
  outfile->Printf("\n");

}

}} /* End Namespaces */
