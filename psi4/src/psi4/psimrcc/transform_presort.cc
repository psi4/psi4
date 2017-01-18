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

/**
 *  @file transform_presort.cc
 *  @ingroup (PSIMRCC)
*/

#include <cmath>
#include <cstdlib>

 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.h"
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/psifiles.h"

#define MAX(i,j) ((i>j) ? i : j)
#define MIN(i,j) ((i>j) ? j : i)
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define four(i,j,k,l) INDEX(INDEX(i,j),INDEX(k,l))

#include "algebra_interface.h"
#include "blas.h"
#include "matrix.h"
#include "index.h"
#include "transform.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;
    extern MemoryManager *memory_manager;

using namespace std;

/**
 * Reads and IWL buffer and sorts the two-electron integrals
 * (pq|rs) as p >= q, r >= s, and pq >= rs.
 * Each symmetry block (determined by pq or rs) is stored in
 * a separate file.  This routine makes psimrcc compatible with
 * transqt and transqt2.
 */
void CCTransform::presort_integrals()
{
  outfile->Printf("\n\n  Presorting two-electron integrals from IWL buffer");
  outfile->Printf("\n    Memory available                       = %14lu bytes",
                  (unsigned long)memory_manager->get_FreeMemory());

  size_t presort_memory = static_cast<size_t>(static_cast<double>(memory_manager->get_FreeMemory())*fraction_of_memory_for_presorting);
  outfile->Printf("\n    Memory available for presorting        = %14lu bytes (%.1f%%)",
                  (unsigned long)presort_memory,fraction_of_memory_for_presorting*100.0);


  // Get the indexing used to store the p >= q pairs
  std::vector<size_t> pairpi = tei_mo_indexing->get_pairpi();

  // Compute the memory for a full in-core presort
  size_t memory_required = 0;
  for(int h = 0; h < pairpi.size(); ++h){
    memory_required += (INDEX(pairpi[h]-1,pairpi[h]-1) + 1) * static_cast<size_t>(sizeof(double));
  }

  outfile->Printf("\n    Memory required for in-core presort    = %14lu bytes",
                  (unsigned long)memory_required);

  if(memory_required < static_cast<size_t>(3) * memory_manager->get_FreeMemory()){
    outfile->Printf("\n    Presorting is not required");
  }

  int first_irrep = 0;
  int last_irrep  = 0;
  bool done       = false;
  while(!done){
    // Determine the batch of irreps to process
    size_t available_presort_memory = presort_memory;

    for(int h = first_irrep; h < moinfo->get_nirreps(); ++h){
      size_t required_memory = (INDEX(pairpi[h]-1,pairpi[h]-1) + 1) * static_cast<size_t>(sizeof(double));
      if(required_memory < available_presort_memory){
        available_presort_memory -= required_memory;
        last_irrep = h + 1;
      }
    }

    // Read the integrals
    presort_blocks(first_irrep,last_irrep);

    // Check if we have done presorting all the irreps
    if(last_irrep >= moinfo->get_nirreps())  done = true;
    first_irrep = last_irrep;
  }

}

void CCTransform::presort_blocks(int first_irrep, int last_irrep)
{
  outfile->Printf("\n    Reading irreps %d -> %d",first_irrep,last_irrep - 1);


  CCIndex* pair_index = blas->get_index("[n>=n]");
  std::vector<size_t> pairpi = pair_index->get_pairpi();

  // Allocate the temporary space
  double** tei_mo;
  allocate1(double*,tei_mo,moinfo->get_nirreps());
  for(int h = first_irrep; h < last_irrep; ++h){
    allocate1(double,tei_mo[h],INDEX(pairpi[h]-1,pairpi[h]-1) + 1);
  }

  // Read all the (frozen + non-frozen) TEI in Pitzer order
  // and store them in a in-core block-matrix
  size_t lastbuf,inbuf,fi;
  size_t elements = 0;
  struct iwlbuf ERIIN;
  iwl_buf_init(&ERIIN,PSIF_MO_TEI,0.0,1,1);
    do{
      lastbuf = ERIIN.lastbuf;
      inbuf   = ERIIN.inbuf;
      fi      = 0;
      for(size_t index = 0; index < inbuf; index++){
        // Compute the [pq] index for this pqrs combination
        size_t p = std::abs(ERIIN.labels[fi]);
        size_t q = ERIIN.labels[fi+1];
        size_t r = ERIIN.labels[fi+2];
        size_t s = ERIIN.labels[fi+3];
        double value = ERIIN.values[index];
        int    irrep = pair_index->get_tuple_irrep(p,q);
        // Fill in only the blocks that fit
        if((first_irrep<=irrep) && (irrep<=last_irrep)){
          size_t pq    = pair_index->get_tuple_rel_index(p,q);
          size_t rs    = pair_index->get_tuple_rel_index(r,s);
          size_t pqrs  = INDEX(pq,rs);
          tei_mo[irrep][pqrs] = value;
        }
        fi += 4;
        elements++;
      }
      if(!lastbuf)
        iwl_buf_fetch(&ERIIN);
    } while(!lastbuf);
  iwl_buf_close(&ERIIN,1);

  outfile->Printf(" (%lu non-zero integrals)", (unsigned long)elements);


  // Write integrals to disk
  for(int h = first_irrep; h < last_irrep; ++h){
    char data_label[80];
    sprintf(data_label,"PRESORTED_TEI_IRREP_%d",h);
    _default_psio_lib_->write_entry(PSIF_PSIMRCC_INTEGRALS,data_label,(char*)&(tei_mo[h][0]),
                     static_cast<size_t>(INDEX(pairpi[h]-1,pairpi[h]-1) + 1) *
                     static_cast<size_t>(sizeof(double)));
  }

  // Deallocate the temporary space
  for(int h = first_irrep; h < last_irrep; ++h){
      release1(tei_mo[h]);
  }
  release1(tei_mo);
}

}} /* End Namespaces */
