#ifndef _psi_src_bin_psimrcc_mrcc_h_
#define _psi_src_bin_psimrcc_mrcc_h_

#include "manybody.h"
#include "heff.h"

namespace psi{ namespace psimrcc{

class Updater;

class CCMRCC : public CCManyBody
{
public:
  // Constructor and destructor
  CCMRCC();
  virtual ~CCMRCC();

  // CCSD
  void compute_energy(Updater* updater);

  // CCSD(T)
  void compute_perturbative_triples();

  // Perturbative correction for CBS
  void compute_first_order_amps();
  void perturbative_cbs();

private:
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

  void build_W_T3_intermediates();
  void build_W_prime_abic_intermediates();
  void build_W_prime_aBIc_intermediates();
  void build_W_prime_AbiC_intermediates();
  void build_W_prime_ABIC_intermediates();

  void build_W_prime_ajki_intermediates();
  void build_W_prime_AjKi_intermediates();
  void build_W_prime_aJkI_intermediates();
  void build_W_prime_AJKI_intermediates();

  void build_W_kija_intermediates();
  void build_W_kiJA_intermediates();
  void build_W_KIja_intermediates();
  void build_W_KIJA_intermediates();

  void build_W_aibc_intermediates();
  void build_W_aIbC_intermediates();
  void build_W_AiBc_intermediates();
  void build_W_AIBC_intermediates();

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

  void build_t3_amplitudes();
  void build_t3_ijkabc_amplitudes_diagrams12();
  void build_t3_ijKabC_amplitudes_diagrams12();
  void build_t3_iJKaBC_amplitudes_diagrams12();
  void build_t3_IJKABC_amplitudes_diagrams12();

  void build_t3_ijkabc_amplitudes_diagram3();
  void build_t3_ijKabC_amplitudes_diagram3();
  void build_t3_iJKaBC_amplitudes_diagram3();
  void build_t3_IJKABC_amplitudes_diagram3();

  void build_t3_ijkabc_amplitudes_diagram4();
  void build_t3_ijKabC_amplitudes_diagram4();
  void build_t3_iJKaBC_amplitudes_diagram4();
  void build_t3_IJKABC_amplitudes_diagram4();

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
