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
#include <cstdio>
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "scf.h"

extern FILE* outfile;

namespace psi{ namespace mcscf{

extern MemoryManager* memory_manager;

void SCF::generate_pairs()
{
  npairs = 0;

  // Count the pairs
  for(int pq_sym = 0; pq_sym < nirreps; ++pq_sym){
    for(int p_sym = 0; p_sym < nirreps; ++p_sym){
      int q_sym = pq_sym ^ p_sym;
      if(p_sym >= q_sym){
        for(int p = 0; p < sopi[p_sym]; ++p){
          for(int q = 0; q < sopi[q_sym]; ++q){
            int p_abs = p + block_offset[p_sym];
            int q_abs = q + block_offset[q_sym];
            if(p_abs >= q_abs){
              pairpi[pq_sym]++;
              npairs++;
            }
          }
        }
      }
    }
  }

  allocate1(int,pairs,2*npairs);

  pair_offset[0] = 0;
  for(int h=1; h< nirreps; ++h)
    pair_offset[h] = pair_offset[h-1] + pairpi[h-1];

  // Store the pairs
  int k = 0;
  npairs = 0;
  for(int pq_sym = 0; pq_sym < nirreps; ++pq_sym){
    for(int p_sym = 0; p_sym < nirreps; ++p_sym){
      int q_sym = pq_sym ^ p_sym;
      if(p_sym >= q_sym){
        for(int p = 0; p < sopi[p_sym]; ++p){
          for(int q = 0; q < sopi[q_sym]; ++q){
            int p_abs = p + block_offset[p_sym];
            int q_abs = q + block_offset[q_sym];
            if(p_abs >= q_abs){
              pair[p_abs][q_abs]     = pair[q_abs][p_abs]     = npairs - pair_offset[pq_sym];
              pair_sym[p_abs][q_abs] = pair_sym[q_abs][p_abs] = pq_sym;
              pairs[k++] = p_abs;
              pairs[k++] = q_abs;
              npairs++;
            }
          }
        }
      }
    }
  }

  outfile->Printf("\n\n  Generated %d pairs\n  Distributed as ",npairs);
  for(int h=0; h< nirreps; ++h)
    outfile->Printf("[%d %s]",pairpi[h],moinfo_scf->get_irr_labs(h));
}

}} /* End Namespaces */
