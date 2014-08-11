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

#include <cmath>
#include <algorithm>

#include <libmoinfo/libmoinfo.h>
#include "transform.h"
#include "matrix.h"
#include <libutil/libutil.h>
#include "algebra_interface.h"
#include "blas.h"

#define CCTRANSFORM_USE_BLAS

#define MAX(i,j) ((i>j) ? i : j)
#define MIN(i,j) ((i>j) ? j : i)
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define four(i,j,k,l) INDEX(INDEX(i,j),INDEX(k,l))

#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libiwl/iwl.h>
#include "psifiles.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void CCTransform::read_integrals_mrpt2()
{
  read_oei_mo_integrals_mrpt2();
  read_tei_mo_integrals_mrpt2();
}

/**
 * Read the one electron MO integrals from an iwl buffer assuming Pitzer ordering and store them in oei_mo
 */
void CCTransform::read_oei_mo_integrals_mrpt2()
{
  read_oei_so_integrals();
  transform_oei_so_integrals();
}

/**
 * Read the two electron MO integrals from an iwl buffer assuming Pitzer ordering and store them in the packed array tei_mo
 */
void CCTransform::read_tei_mo_integrals_mrpt2()
{
  // Read all the (frozen + non-frozen) TEI in Pitzer order
  // and store them in a in-core block-matrix
//   CCIndex* mo_indexing = blas->get_index("[n>=n]");

  double value;
  size_t p,q,r,s,pq,rs,pqrs,irrep;
  size_t ilsti,nbuf,fi,index,elements;
  elements = 0;
  struct iwlbuf ERIIN;
  iwl_buf_init(&ERIIN,PSIF_MO_TEI,0.0,1,1);
    do{
      ilsti = ERIIN.lastbuf;
      nbuf  = ERIIN.inbuf;
      fi    = 0;
      for(index = 0; index < nbuf; index++){
        // Compute the [pq] index for this pqrs combination
        p = abs(ERIIN.labels[fi]);
        q = ERIIN.labels[fi+1];
        r = ERIIN.labels[fi+2];
        s = ERIIN.labels[fi+3];
        value = ERIIN.values[index];

//         outfile->Printf("\n  (%2d %2d|%2d %2d) = %20.12f",p,q,r,s,value);

        integral_map[four(p,q,r,s)]=value;
//         irrep = mo_indexing->get_tuple_irrep(p,q);
//         pq    = mo_indexing->get_tuple_rel_index(p,q);
//         rs    = mo_indexing->get_tuple_rel_index(r,s);
//         pqrs  = INDEX(pq,rs);
//         tei_mo[irrep][pqrs]=value;
        fi += 4;
        elements++;
      }
      if(!ilsti)
        iwl_buf_fetch(&ERIIN);
    } while(!ilsti);
  outfile->Printf("\n    CCTransform: read %lu non-zero integrals (MRPT2)", elements);
  
  iwl_buf_close(&ERIIN,1);
}

double CCTransform::tei_mrpt2(int p, int q, int r, int s)
{
//   outfile->Printf("\n  (%2d %2d|%2d %2d) = %20.15f",p,q,r,s,integral_map[four(p,q,r,s)]);
  return(integral_map[four(p,q,r,s)]);
}

}} /* End Namespaces */
