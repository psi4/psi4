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

/**
 *  @file ccsort.cpp
 *  @ingroup (PSIMRCC)
*/
#include <cmath>
#include <algorithm>

#include <libmoinfo/libmoinfo.h>
#include <libtrans/integraltransform.h>
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
    extern MemoryManager* memory_manager;

using namespace std;

CCSort::CCSort(SortAlgorithm algorithm):
  fraction_of_memory_for_sorting(0.5),nfzc(0),efzc(0.0),frozen_core(0)
{

  // Use libtrans to generate MO basis integrals in Pitzer order
  boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
  std::vector<boost::shared_ptr<MOSpace> > spaces;
  spaces.push_back(MOSpace::all);
  IntegralTransform *ints = new IntegralTransform(wfn, spaces,
                                                  IntegralTransform::Restricted,
                                                  IntegralTransform::IWLOnly,
                                                  IntegralTransform::PitzerOrder,
                                                  IntegralTransform::None);
  ints->transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
  delete ints;
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
  outfile->Printf("\n\n    Frozen-core energy                     = %20.9f",efzc);
  
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
