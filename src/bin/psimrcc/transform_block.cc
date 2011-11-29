#include <cmath>
#include <algorithm>

#include <boost/shared_ptr.hpp>
#include <libmoinfo/libmoinfo.h>
#include <libutil/libutil.h>

#define CCTRANSFORM_USE_BLAS

#define MAX(i,j) ((i>j) ? i : j)
#define MIN(i,j) ((i>j) ? j : i)
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define four(i,j,k,l) INDEX(INDEX(i,j),INDEX(k,l))

#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include "psifiles.h"

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
 * Read at least one block of the two electron MO integrals from an iwl buffer assuming Pitzer ordering and store them in the packed array tei_mo
 */
int CCTransform::read_tei_mo_integrals_block(int first_irrep)
{
  std::vector<size_t> pairpi = tei_mo_indexing->get_pairpi();
  int last_irrep = allocate_tei_mo_block(first_irrep);

  // Write integrals to disk
  for(int h = first_irrep; h < last_irrep; ++h){
    char data_label[80];
    sprintf(data_label,"PRESORTED_TEI_IRREP_%d",h);
    _default_psio_lib_->read_entry(PSIF_PSIMRCC_INTEGRALS,data_label,(char*)&(tei_mo[h][0]),
                     static_cast<size_t>(INDEX(pairpi[h]-1,pairpi[h]-1) + 1) *
                     static_cast<size_t>(sizeof(double)));
  }
  return(last_irrep);
}

/**
 * Allocate as many blocks of the tei_mo array and exit(EXIT_FAILURE) if there is not enough space
 */
int CCTransform::allocate_tei_mo_block(int first_irrep)
{
  if(first_irrep>moinfo->get_nirreps()){
    fprintf(outfile,"\n    Transform: allocate_tei_mo_block() was called with first_irrep > nirreps !");
    fflush(outfile);
    exit(EXIT_FAILURE);
  }

  size_t available_transform_memory = static_cast<size_t>(static_cast<double>(memory_manager->get_FreeMemory())*fraction_of_memory_for_presorting);


  int last_irrep = first_irrep;

  if(tei_mo == NULL){
    // Allocate the tei_mo matrix blocks
    allocate1(double*,tei_mo,moinfo->get_nirreps());
    for(int h = 0; h < moinfo->get_nirreps(); ++h)
      tei_mo[h] = NULL;
  }

  // Find how many irreps we can store in 95% of the free memory
  std::vector<size_t> pairpi = tei_mo_indexing->get_pairpi();
  for(int h = first_irrep; h < moinfo->get_nirreps(); ++h){
    size_t required_memory = (INDEX(pairpi[h]-1,pairpi[h]-1) + 1) * static_cast<size_t>(sizeof(double));
    if(required_memory != 0){
      if(required_memory < available_transform_memory){
        available_transform_memory -= required_memory;
        allocate1(double,tei_mo[h],INDEX(pairpi[h]-1,pairpi[h]-1) + 1)
        zero_arr(tei_mo[h],INDEX(pairpi[h]-1,pairpi[h]-1) + 1);
        last_irrep++;
      }
    }else{
      last_irrep++;
    }
  }
  fprintf(outfile,"\n    Integrals from irreps %d -> %d will be read in core",first_irrep,last_irrep-1);
  if(first_irrep == last_irrep){
    fprintf(outfile,"\n    CCTransform: allocate_tei_mo_block() has not enough memory!");
    fflush(outfile);
    exit(EXIT_FAILURE);
  }
  first_irrep_in_core = first_irrep;
  last_irrep_in_core  = last_irrep;
  return(last_irrep);
}

/**
 * Free the blocks included in the first_irrep->last_irrep range
 */
void CCTransform::free_tei_mo_integrals_block(int first_irrep, int last_irrep)
{
  for(int h = first_irrep; h < last_irrep; ++h){
    if(tei_mo[h] != NULL){
      release1(tei_mo[h]);
    }
  }
  if(last_irrep >= moinfo->get_nirreps()){
    release1(tei_mo);
    tei_mo = NULL;
  }
}

double CCTransform::tei_block(int p, int q, int r, int s)
{
  // Get the (pq|rs) integral
  int irrep(tei_mo_indexing->get_tuple_irrep(MAX(p,q),MIN(p,q)));
  if((first_irrep_in_core <= irrep) && (irrep < last_irrep_in_core))
    return(tei_mo[tei_mo_indexing->get_tuple_irrep(MAX(p,q),MIN(p,q))][INDEX(tei_mo_indexing->get_tuple_rel_index(MAX(p,q),MIN(p,q)),tei_mo_indexing->get_tuple_rel_index(MAX(r,s),MIN(r,s)))]);
  else
    return(0.0);
}

}} /* End Namespaces */


//  Old algorthm for reading ints
//  // Read all the (frozen + non-frozen) TEI in Pitzer order
//  // and store them in a in-core block-matrix
//  CCIndex* mo_indexing = blas->get_index("[n>=n]");
//
//  int last_irrep = allocate_tei_mo_block(first_irrep);
//
//  double value;
//  size_t p,q,r,s,pq,rs,pqrs,irrep;
//  size_t ilsti,nbuf,fi,index,elements;
//  elements = 0;
//  struct iwlbuf ERIIN;
//  iwl_buf_init(&ERIIN,PSIF_MO_TEI,0.0,1,1);
//    do{
//      ilsti = ERIIN.lastbuf;
//      nbuf  = ERIIN.inbuf;
//      fi    = 0;
//      for(index = 0; index < nbuf; index++){
//        // Compute the [pq] index for this pqrs combination
//        p = abs(ERIIN.labels[fi]);
//        q = ERIIN.labels[fi+1];
//        r = ERIIN.labels[fi+2];
//        s = ERIIN.labels[fi+3];
//        value = ERIIN.values[index];
//        irrep = mo_indexing->get_tuple_irrep(p,q);
//        // Fill in only the blocks that fit
//        if((first_irrep<=irrep) && (irrep<last_irrep)){
//          pq    = mo_indexing->get_tuple_rel_index(p,q);
//          rs    = mo_indexing->get_tuple_rel_index(r,s);
//          pqrs  = INDEX(pq,rs);
//          tei_mo[irrep][pqrs]=value;
//        }
//        fi += 4;
//        elements++;
//      }
//      if(!ilsti)
//        iwl_buf_fetch(&ERIIN);
//    } while(!ilsti);
//  fprintf(outfile,"\n    Read %d non-zero integrals", elements);
//  fflush(outfile);
//  iwl_buf_close(&ERIIN,1);
