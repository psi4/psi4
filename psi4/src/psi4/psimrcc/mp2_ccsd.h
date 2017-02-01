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

#ifndef _psi_src_bin_psimrcc_mp2_ccsd_h_
#define _psi_src_bin_psimrcc_mp2_ccsd_h_

#include "manybody.h"

namespace psi{ namespace psimrcc{

/**
	@author Francesco Evangelista <frank@ccc.uga.edu>
*/
class MP2_CCSD : public CCManyBody
{
public:
  MP2_CCSD(SharedWavefunction ref_wfn, Options &options);
  virtual ~MP2_CCSD();
  void compute_mp2_ccsd_energy();
private:
  void add_matrices();
  void read_mp2_ccsd_integrals();
  double compute_energy();
  void compute_mp2_components();
  void compute_mp2_ccsd_components();
  void synchronize_amps();

  /* AMPLITUDES EQUATIONS */
  void build_tau();
  void build_amplitudes();
  void build_t1_ia_amplitudes();
  void build_t1_IA_amplitudes();
  void build_t2_ijab_amplitudes();
  void build_t2_iJaB_amplitudes();
  void build_t2_IJAB_amplitudes();  
  void build_mp2_t2_iJaB_amplitudes();

  /* INTERMEDIATES */

  // F Intermediates
  void build_F_intermediates();
  void build_offdiagonal_F();
  void build_F_ae_intermediates();
  void build_F_AE_intermediates();
  void build_F_me_intermediates();
  void build_F_ME_intermediates();
  void build_F_mi_intermediates();
  void build_F_MI_intermediates();
  void build_F_prime_ae_intermediates();
  void build_F_prime_AE_intermediates();
  void build_F_prime_mi_intermediates();
  void build_F_prime_MI_intermediates();

  // W Intermediates
  void build_W_intermediates();
  void build_W_mNiJ_intermediates();
  void build_W_jbme_intermediates();
  void build_W_JBme_intermediates();
  void build_W_jBmE_intermediates();
  void build_W_jbME_intermediates();
  void build_W_JbMe_intermediates();
  void build_W_JBME_intermediates();

  // Z Intermediates
  void build_Z_intermediates();
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_mp2_ccsd_h_