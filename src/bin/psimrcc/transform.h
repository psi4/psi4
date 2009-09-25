#ifndef _psi_src_bin_psimrcc_cctransform_h
#define _psi_src_bin_psimrcc_cctransform_h

#include <map>

namespace psi{ namespace psimrcc{

class CCIndex;

/**
	@author Francesco A. Evangelista and Andrew C. Simmonett <frank@ccc.uga.edu>
*/
class CCTransform{
public:
  CCTransform();
  ~CCTransform();
  void print();
  // Presorting
  void presort_integrals();
  void read_oei_from_transqt() {read_oei_mo_integrals();}
  void read_integrals_from_transqt() {read_mo_integrals();}
  void read_integrals_mrpt2();
  int  read_tei_mo_integrals_block(int first_irrep);
  void free_tei_mo_integrals_block(int first_irrep, int last_irrep);
  void free_memory();
  void transform_tei_integrals();
  double oei(int p, int q);
  double tei(int p, int q, int r, int s);
  double tei_block(int p, int q, int r, int s);
  double tei_mrpt2(int p, int q, int r, int s);
private:
  size_t*     ioff;
  double**    s_so;
  double**  oei_mo;
  double**  oei_so;
  double**  tei_so;
  double*** tei_half_transformed;
  double**  tei_mo;
  CCIndex*  oei_so_indexing;
  CCIndex*  tei_so_indexing;
  CCIndex*  tei_mo_indexing;

  void read_mo_integrals();
  void read_so_integrals();
  void read_oei_so_integrals();
  void read_oei_mo_integrals();
  void read_oei_mo_integrals_mrpt2();
  void read_tei_so_integrals();
  void read_tei_mo_integrals();
  void read_tei_mo_integrals_mrpt2();

  void transform_oei_so_integrals();
  void transform_tei_so_integrals();

  void allocate_oei_so();
  void allocate_oei_mo();
  void free_oei_mo();
  void free_oei_so();

  void allocate_tei_so();
  void allocate_tei_mo();
  void allocate_tei_half_transformed();

  void free_tei_mo();
  void free_tei_so();
  void free_tei_half_transformed();

  // Block
  int  first_irrep_in_core;
  int  last_irrep_in_core;
  int  allocate_tei_mo_block(int first_irrep);
  std::map <size_t,double> integral_map;

  void presort_blocks(int first_irrep, int last_irrep);

  double    fraction_of_memory_for_presorting;
};

extern CCTransform* trans;

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_cctransform_h
