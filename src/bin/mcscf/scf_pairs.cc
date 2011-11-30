#include <iostream>
#include <cstdio>

#include <libchkpt/chkpt.hpp>
#include <libmoinfo/libmoinfo.h>
#include <libutil/libutil.h>

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

  fprintf(outfile,"\n\n  Generated %d pairs\n  Distributed as ",npairs);
  for(int h=0; h< nirreps; ++h)
    fprintf(outfile,"[%d %s]",pairpi[h],moinfo_scf->get_irr_labs(h));
}

}} /* End Namespaces */
