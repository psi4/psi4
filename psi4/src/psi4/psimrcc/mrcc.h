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

#ifndef _psi_src_bin_psimrcc_mrcc_h_
#define _psi_src_bin_psimrcc_mrcc_h_

#include "psi4/liboptions/liboptions.h"
#include "manybody.h"
#include "heff.h"

namespace psi{ namespace psimrcc{

class Updater;

class CCMRCC : public CCManyBody
{
public:
  // Constructor and destructor
  CCMRCC(SharedWavefunction ref_wfn, Options &options);
  virtual ~CCMRCC();

  // CCSD
  void compute_energy(Updater* updater);

  // CCSD(T)
  void compute_perturbative_triples();

  // Perturbative correction for CBS
  void compute_first_order_amps();
  void perturbative_cbs();

private:
  Options &options_;

  bool ap_correction;

  void add_matrices();

  // These are used to call member functions

  Hamiltonian h_eff;

  void diis(int cycle);
  void diis_save_t_amps();

  void synchronize_amps();
  void compute_delta_amps();
  void add_ccsd_matrices();

  void build_tau_intermediates();

  void build_F_intermediates();
  void build_F_ae_intermediates();
  void build_F_AE_intermediates();
  void build_F_mi_intermediates();
  void build_F_MI_intermediates();
  void build_F_me_intermediates();
  void build_F_ME_intermediates();
  void build_F_prime_ae_intermediates();
  void build_F_prime_AE_intermediates();
  void build_F_prime_mi_intermediates();
  void build_F_prime_MI_intermediates();

  // Intermediates required by triples
  void build_F2_me_intermediates();
  void build_F2_ME_intermediates();


  void build_W_intermediates();
  void build_W_mnij_intermediates();
  void build_W_mNiJ_intermediates();
  void build_W_MNIJ_intermediates();
  void build_W_jbme_intermediates();
  void build_W_JBme_intermediates();
  void build_W_jBmE_intermediates();
  void build_W_jbME_intermediates();
  void build_W_JbMe_intermediates();
  void build_W_JBME_intermediates();

  void build_Z_intermediates();

  void build_t1_amplitudes();
  void build_t1_ia_amplitudes();
  void build_t1_IA_amplitudes();
  void build_t1_amplitudes_triples();
  void build_t1_ia_amplitudes_triples();
  void build_t1_IA_amplitudes_triples();

  void build_t2_amplitudes();
  void build_t2_ijab_amplitudes();
  void build_t2_iJaB_amplitudes();
  void build_t2_IJAB_amplitudes();
  void build_t2_amplitudes_triples();
  void build_t2_ijab_amplitudes_triples_diagram1();
  void build_t2_iJaB_amplitudes_triples_diagram1();
  void build_t2_IJAB_amplitudes_triples_diagram1();
  void build_t2_ijab_amplitudes_triples_diagram2();
  void build_t2_iJaB_amplitudes_triples_diagram2();
  void build_t2_IJAB_amplitudes_triples_diagram2();
  void build_t2_ijab_amplitudes_triples_diagram3();
  void build_t2_iJaB_amplitudes_triples_diagram3();
  void build_t2_IJAB_amplitudes_triples_diagram3();

  void form_similarity_transformed_hamiltonian();

  bool build_diagonalize_Heff(int cycle,double time);
  void build_Heff_diagonal();
  void build_Heff_offdiagonal();

  void update_amps();
  void update_t1_amps();
  void update_t2_amps();
  void update_t3_amps();
  void update_t3_ijkabc_amps();
  void update_t3_ijKabC_amps();
  void update_t3_iJKaBC_amps();
  void update_t3_IJKABC_amps();

  void update_amps_bwccsd();
  void update_t1_amps_bwccsd();
  void update_t2_amps_bwccsd();

  void update_amps_mkccsd();
  void update_t1_amps_mkccsd();
  void update_t2_amps_mkccsd();
  void update_t1_t2_amps_mkccsd();
  void update_t3_amps_mkccsd();
  void update_t3_ijkabc_amps_mkccsd();
  void update_t3_ijKabC_amps_mkccsd();
  void update_t3_iJKaBC_amps_mkccsd();
  void update_t3_IJKABC_amps_mkccsd();

  void update_amps_mkccsd_residual();
  void update_t1_amps_mkccsd_residual();
  void update_t2_amps_mkccsd_residual();
  void update_t1_t2_amps_mkccsd_residual();

  void print_mrccsd_energy(int cycle);

  // CCSD(T)
  std::vector<double>  compute_ooo_triples();
  std::vector<double>  compute_OOO_triples();
  std::vector<double>  compute_ooO_triples();
  std::vector<double>  compute_oOO_triples();


};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_mrcc_h_
