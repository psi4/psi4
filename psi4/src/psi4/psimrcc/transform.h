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

#ifndef _psi_src_bin_psimrcc_cctransform_h
#define _psi_src_bin_psimrcc_cctransform_h

#include <map>

namespace psi{

class IntegralTransform;

namespace psimrcc{

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
  void read_integrals_mrpt2(IntegralTransform *ints);
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
  void read_tei_mo_integrals_mrpt2(IntegralTransform *ints);

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