#ifndef _psi_src_bin_psimrcc_ccsort_h
#define _psi_src_bin_psimrcc_ccsort_h

/**
 *  @file sort.h
 *  @ingroup (PSIMRCC)
*/

#include <iostream>
#include <map>
#include <vector>


namespace psi{ namespace psimrcc{

class CCMatrix;

#include "matrix_types.h"

#ifndef INDEX
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#endif
#define four(i,j,k,l) INDEX(INDEX(i,j),INDEX(k,l))

enum SortAlgorithm {out_of_core_sort,mrpt2_sort};

/**
 *  @class CCSort
 *  @brief Grabs the MO integrals from CCTransform and fills the CCMatrix objects in CCBLAS
*/
class CCSort{
public:
  CCSort(SortAlgorithm algorithm);
  ~CCSort();
private:
  void init();
  void cleanup();

  // In-core algorithm
  void   build_integrals_in_core();
  void   frozen_core_energy_in_core();
  void   sort_integrals_in_core();
  void   form_two_electron_integrals_in_core(MatrixMap::iterator& iter);
  void   form_fock_in_core(MatrixMap::iterator& iter);
  double add_fock_two_in_core(int p, int q, int k, bool exchange);

  // Out-of-core algorithm
  void   build_integrals_out_of_core();
  void   frozen_core_energy_out_of_core();
  void   sort_integrals_out_of_core(int first_irrep, int last_irrep, MatrixBlks& to_be_processed);
  void   form_fock_one_out_of_core(MatrixBlks& to_be_processed);
  void   form_fock_out_of_core(CCMatrix* Matrix, int h);
  void   form_two_electron_integrals_out_of_core(CCMatrix* Matrix, int h);
  double add_fock_two_out_of_core(int p, int q, int k, bool exchange);
  void   setup_out_of_core_list(MatMapIt& mat_it,int& mat_irrep,MatMapIt& mat_end,MatrixBlks&  to_be_processed);
  void   dump_integrals_to_disk(MatrixBlks& to_be_processed);

  // MRPT2 algorithm
  void   build_integrals_mrpt2();
  void   frozen_core_energy_mrpt2();
  void   allocate_and_sort_integrals_mrpt2();
  void   allocate_amplitudes_mrpt2();
  void   form_two_electron_integrals_mrpt2(MatrixMap::iterator& iter);
  void   form_fock_mrpt2(MatrixMap::iterator& iter);
  double add_fock_two_mrpt2(int p, int q, int k, bool exchange);

  // Data
  double    fraction_of_memory_for_sorting;
  int       nfzc;
  double    efzc;
  int*      frozen_core;
};

extern CCSort *sorter;

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_ccsort_h
