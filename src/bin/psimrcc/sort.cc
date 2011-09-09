/**
 *  @file ccsort.cpp
 *  @ingroup (PSIMRCC)
*/
#include <cmath>
#include <algorithm>

#include <libmoinfo/libmoinfo.h>
#include <libutil/libutil.h>
#include <libutil/memory_manager.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include "psifiles.h"

#include "blas.h"
#include "sort.h"
#include "transform.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;

using namespace std;

CCSort::CCSort(SortAlgorithm algorithm):
  fraction_of_memory_for_sorting(0.5),nfzc(0),efzc(0.0),frozen_core(0)
{
  init();

  // Presort the integrals in the CCTransform class
  trans->presort_integrals();


  // Two algorithms for forming the integrals

  // 1. Full in-core algorithm: the transformed integrals and CCMatrix object fit in core
  //   build_integrals_in_core();

  // 2. Out-of-core algorithm: the transformed integrals and the CCMatrix(s) don't fit in
  //    core or they are requested to be out-of-core
  switch (algorithm) {
    case out_of_core_sort :
      build_integrals_out_of_core();
      break;
    case mrpt2_sort :
      build_integrals_mrpt2();
      break;
  }

  moinfo->set_fzcore_energy(efzc);
  fprintf(outfile,"\n\n    Frozen-core energy                     = %20.9f",efzc);
  fflush(outfile);
}

CCSort::~CCSort()
{
  cleanup();
}

/**
 * Initialize the CCSort class
 */
void CCSort::init()
{
  // Find the frozen core orbitals in Pitzer ordering
  nfzc        = moinfo->get_nfocc();
  intvec focc = moinfo->get_focc();
  intvec mopi = moinfo->get_mopi();
  allocate1(int,frozen_core,nfzc);
  int count1  = 0;
  int count2  = 0;
  for(int h = 0; h < moinfo->get_nirreps(); ++h){
    for(int i = 0; i < focc[h]; ++i)
      frozen_core[count1++] = count2 + i;
    count2 += mopi[h];
  }
}

/**
 * Clean the CCSort class
 */
void CCSort::cleanup()
{
  release1(frozen_core);
}

}} /* End Namespaces */
