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

#ifndef dfocc_h
#define dfocc_h

#include "psi4/libmints/wavefunction.h"
#include "psi4/psifiles.h"
#include "arrays.h"
#include "tensors.h"

using namespace std;

namespace psi{

class DIISManager;

namespace dfoccwave {

class DFOCC : public Wavefunction
{

    void common_init();

public:
    DFOCC(SharedWavefunction ref_wfn, Options &options);

    virtual ~DFOCC();

    virtual double compute_energy();

protected:
    void mem_release();
    void get_moinfo();
    void title();
    void title_grad();
    void lambda_title();
    void pt_title();
    void pat_title();
    void pdm_title();
    void ref_energy();
    void mp2_energy();
    void scs_mp2_energy();
    void mp2_direct();
    void t1_1st_sc();
    void t2_1st_sc();
    void t2_1st_gen();
    void t2_1st_scs_sc();
    void t2_1st_scs_gen();
    void Fint_zero();
    void fock();
    void separable_tpdm();
    void sep_tpdm_cc();
    void combine_ref_sep_tpdm();
    void tpdm_tilde();
    void back_trans();
    void tpdm_tilde_cc();
    void back_trans_cc();
    void dfgrad();
    void oei_grad();
    void tei_grad_ref();
    void tei_grad_corr();
    void gfock_oo();
    void gfock_vo();
    void gfock_ov();
    void gfock_vv();
    void gftilde_vv();
    void gfock_cc_oo();
    void gfock_cc_vo();
    void gfock_cc_ov();
    void gfock_cc_vv();
    void idp();
    void idp2();
    void idp_hf();
    void mograd();
    void occ_iterations();
    void kappa_orb_resp();
    void kappa_orb_resp_pcg();
    void orb_resp_pcg_rhf();
    void orb_resp_pcg_uhf();
    void kappa_diag_hess();
    void kappa_qchf();
    void update_mo();
    void update_hfmo();
    void semi_canonic();
    void canonic();
    void diagonal_mohess_vo();
    void diagonal_mohess_oo();
    void approx_diag_mohess_vo();
    void approx_diag_mohess_oo();
    void approx_diag_hf_mohess_vo();
    void approx_diag_hf_mohess_oo();
    void approx_diag_ekt_mohess_vo();
    void approx_diag_ekt_mohess_oo();
    void prepare4grad();
    void z_vector();
    void z_vector_pcg();
    void z_vector_cg();
    void z_vector_solver();
    void pcg_solver_rhf();
    void pcg_solver_uhf();
    void cg_solver();
    void zvec_solver_rhf();
    void zvec_solver_uhf();
    void effective_pdms();
    void effective_gfm();
    void effective_pdm_gfm();
    void effective_mograd();
    void fc_grad_terms();
    void z_vector_fc();
    void oeprop();
    void s2_response();
    void s2_lagrangian();
    void gwh();
    void qchf();
    void mo_coeff_blocks();
    void ccl_energy();
    void ccl_energy2();
    void save_mo_to_wfn();
    void remove_binary_file(int fileno);

    void diis(int dimvec, SharedTensor2d &vecs, SharedTensor2d &errvecs, SharedTensor1d &vec_new, SharedTensor1d &errvec_new);
    void sigma_rhf(SharedTensor1d& sigma, SharedTensor1d& p_vec);
    void sigma_uhf(SharedTensor1d& sigma_A, SharedTensor1d& sigma_B, SharedTensor1d& p_vecA, SharedTensor1d& p_vecB);
    void sigma_orb_resp_rhf(SharedTensor1d& sigma, SharedTensor1d& p_vec);
    void build_rhf_mohess(SharedTensor2d& Aorb_);
    void build_uhf_mohess(SharedTensor2d& Aorb_);
    void t2_rmp2_direct(SharedTensor2d& T);
    void u2_rmp2_direct(SharedTensor2d& U);
    void u2_rmp2_direct(SharedTensor2d& T, SharedTensor2d& U);
    void t2AA_ump2_direct(SharedTensor2d& T);
    void t2BB_ump2_direct(SharedTensor2d& T);
    void t2AB_ump2_direct(SharedTensor2d& T);

    // Conventional integrals for DF-BASIS-CC
    void tei_ijkl_chem();
    void tei_ijka_chem();
    void tei_ijab_chem();
    void tei_iajb_chem();
    void tei_oooo_chem();
    void tei_ooov_chem();
    void tei_oovv_chem();
    void tei_ovov_chem();

    void tei_ijkl_phys();
    void tei_ijka_phys();
    void tei_ijab_phys();
    void tei_iajb_phys();
    void tei_oooo_phys();
    void tei_ooov_phys();
    void tei_oovv_phys();
    void tei_ovov_phys();

    void tei_ijkl_anti_symm();
    void tei_ijka_anti_symm();
    void tei_ijab_anti_symm();
    void tei_iajb_anti_symm();
    void tei_oooo_anti_symm();
    void tei_ooov_anti_symm();
    void tei_oovv_anti_symm();
    void tei_ovov_anti_symm();

    // Conventional integrals for DF-BASIS-CC with direct algorithm
    // Integrals in chemist notations
    void tei_chem_direct(SharedTensor2d &K, SharedTensor2d &L, SharedTensor2d &M);
    void tei_ijkl_chem_directAA(SharedTensor2d &K);
    void tei_ijkl_chem_directBB(SharedTensor2d &K);
    void tei_ijkl_chem_directAB(SharedTensor2d &K);

    void tei_ijka_chem_directAA(SharedTensor2d &K);
    void tei_ijka_chem_directBB(SharedTensor2d &K);
    void tei_ijka_chem_directAB(SharedTensor2d &K);
    void tei_iajk_chem_directAB(SharedTensor2d &K);

    void tei_ijab_chem_directAA(SharedTensor2d &K);
    void tei_ijab_chem_directBB(SharedTensor2d &K);
    void tei_ijab_chem_directAB(SharedTensor2d &K);
    void tei_abij_chem_directAB(SharedTensor2d &K);

    void tei_iajb_chem_directAA(SharedTensor2d &K);
    void tei_iajb_chem_directBB(SharedTensor2d &K);
    void tei_iajb_chem_directAB(SharedTensor2d &K);

    void tei_oooo_chem_directAA(SharedTensor2d &K);
    void tei_oooo_chem_directBB(SharedTensor2d &K);
    void tei_oooo_chem_directAB(SharedTensor2d &K);

    void tei_ooov_chem_directAA(SharedTensor2d &K);
    void tei_ooov_chem_directBB(SharedTensor2d &K);
    void tei_ooov_chem_directAB(SharedTensor2d &K);
    void tei_ovoo_chem_directAB(SharedTensor2d &K);

    void tei_oovv_chem_directAA(SharedTensor2d &K);
    void tei_oovv_chem_directBB(SharedTensor2d &K);
    void tei_oovv_chem_directAB(SharedTensor2d &K);
    void tei_vvoo_chem_directAB(SharedTensor2d &K);

    void tei_ovov_chem_directAA(SharedTensor2d &K);
    void tei_ovov_chem_directBB(SharedTensor2d &K);
    void tei_ovov_chem_directAB(SharedTensor2d &K);

    // Integrals in physist notations
    void tei_phys_direct(SharedTensor2d &I, SharedTensor2d &K, SharedTensor2d &L, SharedTensor2d &M);
    void tei_ijkl_phys_directAA(SharedTensor2d &K);
    void tei_ijkl_phys_directBB(SharedTensor2d &K);
    void tei_ijkl_phys_directAB(SharedTensor2d &K);

    void tei_ijka_phys_directAA(SharedTensor2d &K);
    void tei_ijka_phys_directBB(SharedTensor2d &K);
    void tei_ijka_phys_directAB(SharedTensor2d &K);
    void tei_ijak_phys_directAB(SharedTensor2d &K);

    void tei_ijab_phys_directAA(SharedTensor2d &K);
    void tei_ijab_phys_directBB(SharedTensor2d &K);
    void tei_ijab_phys_directAB(SharedTensor2d &K);

    void tei_iajb_phys_directAA(SharedTensor2d &K);
    void tei_iajb_phys_directBB(SharedTensor2d &K);
    void tei_iajb_phys_directAB(SharedTensor2d &K);
    void tei_aibj_phys_directAB(SharedTensor2d &K);

    void tei_oooo_phys_directAA(SharedTensor2d &K);
    void tei_oooo_phys_directBB(SharedTensor2d &K);
    void tei_oooo_phys_directAB(SharedTensor2d &K);

    void tei_ooov_phys_directAA(SharedTensor2d &K);
    void tei_ooov_phys_directBB(SharedTensor2d &K);
    void tei_ooov_phys_directAB(SharedTensor2d &K);
    void tei_oovo_phys_directAB(SharedTensor2d &K);

    void tei_oovv_phys_directAA(SharedTensor2d &K);
    void tei_oovv_phys_directBB(SharedTensor2d &K);
    void tei_oovv_phys_directAB(SharedTensor2d &K);

    void tei_ovov_phys_directAA(SharedTensor2d &K);
    void tei_ovov_phys_directBB(SharedTensor2d &K);
    void tei_ovov_phys_directAB(SharedTensor2d &K);
    void tei_vovo_phys_directAB(SharedTensor2d &K);

    // Anti-symmetrized integrals
    void tei_pqrs_anti_symm_direct(SharedTensor2d &K, SharedTensor2d &L);
    void tei_pqrs2_anti_symm_direct(SharedTensor2d &K, SharedTensor2d &L);
    void tei_pqrs3_anti_symm_direct(SharedTensor2d &K, SharedTensor2d &L, SharedTensor2d &M);
    void tei_cs1_anti_symm_direct(SharedTensor2d &I, SharedTensor2d &J, SharedTensor2d &K);
    void tei_cs2_anti_symm_direct(SharedTensor2d &I, SharedTensor2d &J, SharedTensor2d &K);
    void tei_cs3_anti_symm_direct(SharedTensor2d &I, SharedTensor2d &J, SharedTensor2d &K);
    void tei_cs4_anti_symm_direct(SharedTensor2d &I, SharedTensor2d &J, SharedTensor2d &K);

    // Conventional integrals for DF-BASIS-SCF
    void tei_oooo_chem_ref();
    void tei_ooov_chem_ref();
    void tei_oovv_chem_ref();
    void tei_ovov_chem_ref();

    void tei_oooo_phys_ref();
    void tei_ooov_phys_ref();
    void tei_oovv_phys_ref();
    void tei_ovov_phys_ref();

    void tei_oooo_anti_symm_ref();
    void tei_ooov_anti_symm_ref();
    void tei_oovv_anti_symm_ref();
    void tei_ovov_anti_symm_ref();

    // Conventional integrals for DF-BASIS-SCF with direct algorithm
    // Integrals in chemist notations
    void tei_oooo_chem_ref_directAA(SharedTensor2d &K);
    void tei_oooo_chem_ref_directBB(SharedTensor2d &K);
    void tei_oooo_chem_ref_directAB(SharedTensor2d &K);

    void tei_ooov_chem_ref_directAA(SharedTensor2d &K);
    void tei_ooov_chem_ref_directBB(SharedTensor2d &K);
    void tei_ooov_chem_ref_directAB(SharedTensor2d &K);
    void tei_ovoo_chem_ref_directAB(SharedTensor2d &K);

    void tei_oovv_chem_ref_directAA(SharedTensor2d &K);
    void tei_oovv_chem_ref_directBB(SharedTensor2d &K);
    void tei_oovv_chem_ref_directAB(SharedTensor2d &K);
    void tei_vvoo_chem_ref_directAB(SharedTensor2d &K);

    void tei_ovov_chem_ref_directAA(SharedTensor2d &K);
    void tei_ovov_chem_ref_directBB(SharedTensor2d &K);
    void tei_ovov_chem_ref_directAB(SharedTensor2d &K);

    void tei_vovo_chem_ref_directAA(SharedTensor2d &K);
    void tei_vovo_chem_ref_directBB(SharedTensor2d &K);
    void tei_vovo_chem_ref_directAB(SharedTensor2d &K);

    // Integrals in physist notations
    void tei_oooo_phys_ref_directAA(SharedTensor2d &K);
    void tei_oooo_phys_ref_directBB(SharedTensor2d &K);
    void tei_oooo_phys_ref_directAB(SharedTensor2d &K);

    void tei_ooov_phys_ref_directAA(SharedTensor2d &K);
    void tei_ooov_phys_ref_directBB(SharedTensor2d &K);
    void tei_ooov_phys_ref_directAB(SharedTensor2d &K);
    void tei_oovo_phys_ref_directAB(SharedTensor2d &K);

    void tei_oovv_phys_ref_directAA(SharedTensor2d &K);
    void tei_oovv_phys_ref_directBB(SharedTensor2d &K);
    void tei_oovv_phys_ref_directAB(SharedTensor2d &K);

    void tei_ovov_phys_ref_directAA(SharedTensor2d &K);
    void tei_ovov_phys_ref_directBB(SharedTensor2d &K);
    void tei_ovov_phys_ref_directAB(SharedTensor2d &K);
    void tei_vovo_phys_ref_directAB(SharedTensor2d &K);

    // df
    void df();
    void df_corr();
    void df_ref();
    void trans_corr();
    void trans_ref();
    void trans_mp2();
    void formJ(std::shared_ptr<BasisSet> auxiliary_, std::shared_ptr<BasisSet> zero);
    void formJ_ref(std::shared_ptr<BasisSet> auxiliary_, std::shared_ptr<BasisSet> zero);
    void b_so(std::shared_ptr<BasisSet> primary_, std::shared_ptr<BasisSet> auxiliary_, std::shared_ptr<BasisSet> zero);
    void b_so_ref(std::shared_ptr<BasisSet> primary_, std::shared_ptr<BasisSet> auxiliary_, std::shared_ptr<BasisSet> zero);
    void b_oo();
    void b_ov();
    void b_vv();
    void b_ij();
    void b_ia();
    void b_ab();
    void c_oo();
    void c_ov();
    void c_vv();
    void c_ij();
    void c_ia();
    void c_ab();
    void b_oo_ref();
    void b_ov_ref();
    void b_vv_ref();
    void c_oo_ref();
    void c_ov_ref();
    void c_vv_ref();
    void trans_oei();
    void pair_index();
    void fock_so();
    void ref_grad();
    void b_so_non_zero(); // form non-zero so-basis df ints

    // Cholesky
    void cd_ints();
    void trans_cd();
    void trans_cd_mp2();
    void b_oo_cd();
    void b_ov_cd();
    void b_vv_cd();
    void b_ij_cd();
    void b_ia_cd();
    void b_ab_cd();
    void b_so_non_zero_cd();
    void cd_aob_cints();
    void cd_aob_xints();
    void cd_abcd_cints();
    void cd_abcd_xints();
    void ldl_abcd_ints();
    void ldl_pqrs_ints(int dim1, int dim2, SharedTensor2d &bQ);

    // OMP2
    void omp2_manager();
    void mp2_manager();
    void cd_omp2_manager();
    void cd_mp2_manager();
    void omp2_opdm();
    void omp2_tpdm();
    void mp2l_energy();

    // OMP3
    void omp3_manager();
    void mp3_manager();
    void omp3_manager_cd();
    void mp3_manager_cd();
    void omp3_opdm();
    void omp3_tpdm();
    void t2_2nd_sc();
    void t2_2nd_gen();
    void mp3_WmnijT2();
    void mp3_WmbejT2();
    void mp3_WabefT2();
    void mp3_WmnijT2AA();
    void mp3_WmnijT2BB();
    void mp3_WmnijT2AB();
    void mp3_WmbejT2AA();
    void mp3_WmbejT2BB();
    void mp3_WmbejT2AB();
    void mp3_WabefT2AA();
    void mp3_WabefT2BB();
    void mp3_WabefT2AB();
    void mp3_pdm_3index_intr();
    void mp3_t2_1st_sc();
    void mp3_t2_1st_gen();
    void mp3l_energy();

    // OMP2.5
    void omp2_5_manager();
    void mp2_5_manager();
    void omp2_5_manager_cd();
    void mp2_5_manager_cd();

    // OLCCD
    void olccd_manager();
    void lccd_manager();
    void olccd_manager_cd();
    void lccd_manager_cd();
    void olccd_tpdm();
    void lccd_iterations();
    void lccd_t2_1st_sc();
    void lccd_t2_amps();
    void lccd_WmnijT2();
    void lccd_WmbejT2();
    void lccd_WabefT2();
    void lccd_WmnijT2AA();
    void lccd_WmnijT2BB();
    void lccd_WmnijT2AB();
    void lccd_WmbejT2AA();
    void lccd_WmbejT2BB();
    void lccd_WmbejT2AB();
    void lccd_WabefT2AA();
    void lccd_WabefT2BB();
    void lccd_WabefT2AB();
    void lccd_pdm_3index_intr();
    void lccdl_energy();

    // CCSD
    void ccsd_manager();
    void ccsd_manager_cd();
    void ccsd_mp2();
    void ccsd_iterations();
    void ccsd_3index_intr();
    void ccsd_F_intr();
    void ccsd_WmnijT2();
    void ccsd_WijamT2();
    void ccsd_WijamT2_high_mem();
    void ccsd_WmbejT2();
    void ccsd_WabefT2();
    void ccsd_Wabef2T2();
    void ccsd_WabefT2_high_mem();
    void ccsd_WabefT2_ao_basis();
    void ccsd_WabefT2_ldl();
    void ccsd_WabefT2_cd();
    void ccsd_t1_amps();
    void ccsd_t2_amps();
    void ccsd_energy();
    void ccsd_u2_amps(SharedTensor2d &U, SharedTensor2d &T);
    void ccsd_u2_amps2(SharedTensor2d &U, SharedTensor2d &T);
    void ccsd_t2_prime_amps(SharedTensor2d &U, SharedTensor2d &T);
    void ccsd_t2_prime_amps2(SharedTensor2d &U, SharedTensor2d &T);
    void ccsd_tau_amps(SharedTensor2d &U, SharedTensor2d &T);
    void ccsd_tau_tilde_amps(SharedTensor2d &U, SharedTensor2d &T);
    void ccsd_mp2_low();
    void ccsd_iterations_low();
    void ccsd_3index_intr_low();
    void ccsd_F_intr_low();
    void ccsd_WmnijT2_low();
    void ccsd_WijamT2_low();
    void ccsd_WmbejT2_low();
    void ccsd_WabefT2_low();
    void ccsd_Wabef2T2_low();
    void ccsd_t1_amps_low();
    void ccsd_t2_amps_low();

    // CCSDL
    void ccsdl_l1_amps();
    void ccsdl_l2_amps();
    void ccsdl_iterations();
    void ccsdl_3index_intr();
    void ccsdl_Wmbej();         // OVVO
    void ccsdl_Wmbje();         // OVOV
    void ccsdl_Wmnie();         // OOOV
    void ccsdl_Wmbij();         // OVOO
    void ccsdl_Wmnij();         // OOOO
    void ccsdl_WmbejL2();
    void ccsdl_VmnijL2();
    void ccsdl_WijmnL2();
    void ccsdl_WabefL2();
    void ccsdl_WabefL2_high_mem();
    void ccsdl_Wmnie_direct(SharedTensor2d &W);
    void ccsdl_tau_amps(SharedTensor2d &U, SharedTensor2d &T);
    void ccsdl_LijmeL2_high_mem();

    // CCSD Density
    void ccsd_pdm_3index_intr();
    void ccsd_pdm_yQia();
    void ccsd_opdm();
    void ccsd_tpdm();

    // CCSD(T)
    void ccsd_canonic_triples();
    void ccsd_canonic_triples_hm();
    void ccsd_canonic_triples_disk();
    void ccsd_t_manager();
    void ccsd_t_manager_cd();

    // Lambda-CCSD(T)
    void ccsdl_canonic_triples_disk();
    void ccsdl_t_manager();
    void ccsdl_t_manager_cd();

    // CCD
    void ccd_manager();
    void ccd_manager_cd();
    void ccd_mp2();
    void ccd_iterations();
    void ccd_3index_intr();
    void ccd_F_intr();
    void ccd_WmnijT2();
    void ccd_WmbejT2();
    void ccd_WabefT2();
    void ccd_WabefT2_high_mem();
    void ccd_WabefT2_ldl();
    void ccd_t2_amps();
    void ccd_mp2_low();
    void ccd_iterations_low();
    void ccd_3index_intr_low();
    void ccd_F_intr_low();
    void ccd_WmnijT2_low();
    void ccd_WmbejT2_low();
    void ccd_WabefT2_low();
    void ccd_t2_amps_low();

    // CCDL
    void ccdl_l2_amps();
    void ccdl_iterations();
    void ccdl_3index_intr();
    void ccdl_Wmbej();         // OVVO
    void ccdl_Wmbje();         // OVOV
    void ccdl_Wmnij();         // OOOO
    void ccdl_WmbejL2();
    void ccdl_VmnijL2();
    void ccdl_WijmnL2();
    void ccdl_WabefL2();

    // CCD Density
    void ccd_pdm_3index_intr();
    void ccd_pdm_yQia();
    void ccd_opdm();
    void ccd_tpdm();

    // QCHF
    void qchf_manager();

    // orbital pairs
    int so_pair_idx(int i, int j);
    int mo_pair_idx(int i, int j);
    int oo_pair_idxAA(int i, int j);
    int ij_pair_idxAA(int i, int j);
    int vv_pair_idxAA(int i, int j);
    int ab_pair_idxAA(int i, int j);
    int ov_pair_idxAA(int i, int j);
    int ia_pair_idxAA(int i, int j);
    int vo_pair_idxAA(int i, int j);
    int ai_pair_idxAA(int i, int j);
    int get_rotation_block(string rotblock);

    // DIIS
    std::shared_ptr<DIISManager> ccsdDiisManager;
    std::shared_ptr<DIISManager> ccsdDiisManagerAA;
    std::shared_ptr<DIISManager> ccsdDiisManagerBB;
    std::shared_ptr<DIISManager> ccsdDiisManagerAB;
    std::shared_ptr<DIISManager> ccsdlDiisManager;

    // Gradients
    std::map<std::string, SharedMatrix> gradients;
    std::vector<std::string> gradient_terms;

     int natom;
     int nmo;		// Number of MOs
     int nao;		// Number of AOs
     int nao_nz;	// Number of non-zero AOs
     int ndf_nz;	// Number of non-zero DF ints in AO-basis
     int nso;		// Number of SOs
     int noccA;		// Number of alpha occupied orbitals
     int noccB;		// Number of beta occupied orbitals
     int nvirA;		// Number of alpha virtual orbitals
     int nvirB;		// Number of beta virtual orbitals
     int naoccA;	// Number of active alpha occupied orbitals
     int naoccB;	// Number of active beta occupied orbitals
     int namo;		// Number of active SOs
     int navirA;	// Number of active alpha virtual orbitals
     int navirB;	// Number of active beta virtual orbitals
     int nirreps;	// Number of irreducible representations
     int nshell;	// Number of shells
     int nfrzc; 	// Number of frozen cores
     int nfrzv; 	// Number of frozen virtuals
     int npop; 		// Number of populated orbitals: npop=nmo-nfrzv
     int dimtei;	// dimension of tei in pitzer order for all integrals
     int ntri; 		// square matrix dimension (nmo) -> pitzer order
     int ntri_so;	// square matrix dimension (nso) -> pitzer order
     int ntri_ijAA;
     int ntri_ijBB;
     int ntri_abAA;
     int ntri_abBB;
     int ntri_anti_ijAA;
     int ntri_anti_ijBB;
     int ntri_anti_abAA;
     int ntri_anti_abBB;
     int nQ;          // numer of aux-basis
     int nQ_ref;      // numer of aux-basis for DF_BASIS_SCF
     int nQ_cd;      // numer of aux-basis functions for LDL
     int nso2_;       // nso * nso
     int naocc2AA;     // # of active OO pairs
     int naocc2AB;     // # of active OO pairs
     int naocc2BB;     // # of active OO pairs
     int navir2AA;     // # of active VV pairs
     int navir2AB;     // # of active VV pairs
     int navir2BB;     // # of active VV pairs
     int navoAA;       // # of active VO pairs
     int navoAB;       // # of active VO pairs
     int navoBA;       // # of active VO pairs
     int navoBB;       // # of active VO pairs
     int nvoAA;        // # of all VO pairs
     int nvoAB;        // # of all VO pairs
     int nvoBA;        // # of all VO pairs
     int nvoBB;        // # of all VO pairs
     int nocc2AA;     // # of all OO pairs
     int nocc2AB;     // # of all OO pairs
     int nocc2BB;     // # of all OO pairs
     int nvir2AA;     // # of all VV pairs
     int nvir2AB;     // # of all VV pairs
     int nvir2BB;     // # of all VV pairs
     int nidpA;       // # of alpha independent pairs
     int nidpB;       // # of beta independent pairs
     int nidp;
     int nidp_tot;
     int idp_returnA;
     int idp_returnB;
     int nvar;
     int pcg_conver;

     int exp_cutoff;
     int exp_int_cutoff;
     int multp;
     int charge;
     int print_;
     int conver;
     int cc_maxiter;
     int mo_maxiter;
     int pcg_maxiter;
     int num_vecs;              // Number of vectors used in diis (diis order)
     int do_diis_;
     int itr_diis;
     int itr_occ;
     int itr_pcg;
     int time4grad;             // If 0 it is not the time for grad, if 1 it is the time for grad
     int cc_maxdiis_;           // MAX Number of vectors used in CC diis
     int cc_mindiis_;           // MIN Number of vectors used in CC diis
     int trans_ab;              // 0 means do not transform, 1 means do transform B(Q, ab)
     int mo_optimized;          // 0 means MOs are not optimized, 1 means Mos are optimized
     int orbs_already_opt;      // 0 false, 1 true
     int orbs_already_sc;       // 0 false, 1 true
     int nincore_amp;

     ULI memory;
     double memory_mb;
     double cost_ampAA;          // Mem required for the amplitudes
     double cost_ampBB;          // Mem required for the amplitudes
     double cost_ampAB;          // Mem required for the amplitudes
     double cost_amp;            // Mem required for the amplitudes
     double cost_df;             // Mem required for the df integrals
     double cost_3amp;
     double cost_4amp;
     double cost_5amp;
     double cost_ppl_hm;        // Mem req. for high mem evaluation of 4-virtuals exchange term
     double cost_triples_iabc;   // Mem req. for high mem evaluation of (ia|bc) used in (T)

     // Common
     double Enuc;
     double Eelec;
     double Escf;
     double Eref;
     double Etotal;
     double Emp2;
     double Emp2_t1;
     double Emp2BB;
     double Emp2AA;
     double Emp2AB;
     double Emp2L;
     double Emp2L_old;
     double Ecorr;
     double EcorrL;
     double EccL;
     double Ecc_rdm;
     double Escsmp2;
     double Escsmp2BB;
     double Escsmp2AA;
     double Escsmp2AB;
     double Esosmp2AB;
     double Esosmp2;
     double Escsnmp2;
     double Escsnmp2AA;
     double Escsnmp2BB;
     double Escsnmp2AB;
     double DE;
     double cutoff;
     double int_cutoff_;
     double tol_Eod;
     double tol_grad;
     double tol_t2;
     double tol_pcg;
     double tol_ldl;
     double step_max;
     double mograd_max;
     double biggest_mograd;
     double biggest_mogradA;
     double biggest_mogradB;
     double biggest_kappa;
     double biggest_kappaA;
     double biggest_kappaB;
     double lshift_parameter;
     double os_scale;
     double ss_scale;
     double sos_scale;
     double sos_scale2;
     double e3_scale;
     double rms_t2;
     double rms_t2AA;
     double rms_t2AB;
     double rms_t2BB;
     double rms_l2;
     double mu_ls;
     double sc_ls;
     double rms_wog;
     double rms_wogA;
     double rms_wogB;
     double rms_kappa;
     double rms_kappaA;
     double rms_kappaB;
     double msd_oo_scale;
     double reg_param;
     double s2_resp;
     double s2_proj;
     double s2_lag;
     double s2_ref;

     // OMP3
     double Emp3;
     double Emp3BB;
     double Emp3AA;
     double Emp3AB;
     double Emp3L;
     double Emp3L_old;
     double Escsmp3BB;
     double Escsmp3AA;
     double Escsmp3AB;
     double Escsmp3;
     double Esosmp3AB;
     double Esosmp3;

     // OLCCD
     double Elccd;
     double Elccd_old;
     double ElccdAA;
     double ElccdBB;
     double ElccdAB;
     double ElccdL;
     double ElccdL_old;

     // CCSD
     double Eccsd;
     double Eccsd_old;
     double EccsdAA;
     double EccsdBB;
     double EccsdAB;
     double rms_t1;
     double rms_t1A;
     double rms_t1B;
     double EccsdL_old;
     double EccsdL;
     double EccsdLAA;
     double EccsdLBB;
     double EccsdLAB;
     double Eccsd_t;
     double E_t;
     double Eccsd_at;
     double E_at;

     // CCD
     double Eccd;
     double Eccd_old;
     double EccdAA;
     double EccdBB;
     double EccdAB;
     double EccdL_old;
     double EccdL;
     double EccdLAA;
     double EccdLBB;
     double EccdLAB;

     string wfn;
     string reference;
     string reference_;
     string jobtype;
     string dertype;
     string basis;
     string level_shift;
     string lineq;
     string orth_type;
     string natorb;
     string semicanonic;
     string opt_method;
     string hess_type;
     string occ_orb_energy;
     string do_scs;		// Spin-Component-Scaling
     string do_sos;		// Spin-Opposite-Scaling
     string scs_type_;
     string sos_type_;
     string wfn_type_;
     string orb_resp_solver_;
     string pcg_beta_type_;
     string ekt_ip_;
     string orb_opt_;
     string rotation_blocks;
     string conv_tei_type;
     string regularization;
     string do_cd;
     string read_scf_3index;
     string freeze_core_;
     string oeprop_;
     string comput_s2_;
     string mp2_amp_type_;
     string guess_type_;
     string qchf_;
     string cc_lambda_;
     string Wabef_type_;
     string triples_iabc_type_;

     bool df_ints_incore;
     bool t2_incore;
     bool do_ppl_hm;
     bool do_triples_hm;

     double **C_pitzerA;
     double **C_pitzerB;
     double **J_mhalf;

     // Common
     SharedTensor2d CmoA;              // C(mu, p)
     SharedTensor2d CmoB;              // C(mu, p)
     SharedTensor2d Cmo_refA;          // Reference (initial) MOs
     SharedTensor2d Cmo_refB;          // Reference (initial) MOs
     SharedTensor2d CaoccA;            // C(mu, i) active
     SharedTensor2d CaoccB;            // C(mu, i) active
     SharedTensor2d CavirA;            // C(mu, a) active
     SharedTensor2d CavirB;            // C(mu, a) active
     SharedTensor2d CoccA;             // C(mu, i) all
     SharedTensor2d CoccB;             // C(mu, i) all
     SharedTensor2d CvirA;             // C(mu, a) all
     SharedTensor2d CvirB;             // C(mu, a) all
     SharedTensor2d HmoA;
     SharedTensor2d HmoB;
     SharedTensor2d FockA;
     SharedTensor2d FockB;
     SharedTensor2d Hso;
     SharedTensor2d Sso;
     SharedTensor2d Dso;
     SharedTensor2d DsoA;
     SharedTensor2d FsoA;
     SharedTensor2d FsoB;
     SharedTensor2d Wso;
     SharedTensor2d DQmatA;
     SharedTensor2d HooA;
     SharedTensor2d HooB;
     SharedTensor2d HovA;
     SharedTensor2d HovB;
     SharedTensor2d HvoA;
     SharedTensor2d HvoB;
     SharedTensor2d HvvA;
     SharedTensor2d HvvB;
     SharedTensor2d FooA;          // Fock OO block
     SharedTensor2d FooB;          // Fock oo block
     SharedTensor2d FovA;          // Fock OV block
     SharedTensor2d FovB;          // Fock ov block
     SharedTensor2d FvoA;          // Fock VO block
     SharedTensor2d FvoB;          // Fock vo block
     SharedTensor2d FvvA;          // Fock VV block
     SharedTensor2d FvvB;          // Fock vv block
     SharedTensor1d eigooA;
     SharedTensor1d eigooB;
     SharedTensor1d eigvvA;
     SharedTensor1d eigvvB;
     SharedTensor1d eps_orbA;
     SharedTensor1d eps_orbB;

     // DF Integrals
     SharedTensor2d Jmhalf;             // J Metric DF_BASIS_CC (RI)
     SharedTensor2d bQso;               // b(Q|mu nu) from DF_BASIS_CC (RI)
     SharedTensor2d bQnoA;              // b(Q|mu i)
     SharedTensor2d bQnoB;              // b(Q|mu i)
     SharedTensor2d bQnvA;              // b(Q|mu a)
     SharedTensor2d bQnvB;              // b(Q|mu a)
     SharedTensor2d bQooA;              // b(Q|i j) : all
     SharedTensor2d bQooB;              // b(Q|i j) : all
     SharedTensor2d bQovA;              // b(Q|i a) : all
     SharedTensor2d bQovB;              // b(Q|i a) : all
     SharedTensor2d bQvvA;              // b(Q|a b) : all
     SharedTensor2d bQvvB;              // b(Q|a b) : all
     SharedTensor2d bQijA;              // b(Q|i j) : active
     SharedTensor2d bQijB;              // b(Q|i j) : active
     SharedTensor2d bQiaA;              // b(Q|i a) : active
     SharedTensor2d bQiaB;              // b(Q|i a) : active
     SharedTensor2d bQabA;              // b(Q|a b) : active
     SharedTensor2d bQabB;              // b(Q|a b) : active

     SharedTensor2d cQso;               // c(Q|mu nu) from DF_BASIS_CC (RI)
     SharedTensor2d cQnoA;              // c(Q|mu i)
     SharedTensor2d cQnoB;              // c(Q|mu i)
     SharedTensor2d cQnvA;              // c(Q|mu a)
     SharedTensor2d cQnvB;              // c(Q|mu a)
     SharedTensor2d cQooA;              // c(Q|i j) : all
     SharedTensor2d cQooB;              // c(Q|i j) : all
     SharedTensor2d cQovA;              // c(Q|i a) : all
     SharedTensor2d cQovB;              // c(Q|i a) : all
     SharedTensor2d cQvvA;              // c(Q|a b) : all
     SharedTensor2d cQvvB;              // c(Q|a b) : all
     SharedTensor2d cQijA;              // c(Q|i j) : active
     SharedTensor2d cQijB;              // c(Q|i j) : active
     SharedTensor2d cQiaA;              // c(Q|i a) : active
     SharedTensor2d cQiaB;              // c(Q|i a) : active
     SharedTensor2d cQabA;              // c(Q|a b) : active
     SharedTensor2d cQabB;              // c(Q|a b) : active

     // DF OPDM
     SharedTensor2d G1c_ij;
     SharedTensor2d G1c_ab;
     SharedTensor2d G1c_oo;
     SharedTensor2d G1c_ov;
     SharedTensor2d G1c_vo;
     SharedTensor2d G1c_vv;
     SharedTensor2d G1c;               // Correlation part of OPDM
     SharedTensor2d G1;                // Full OPDM (MO)
     SharedTensor2d G1ao;              // Full OPDM (AO)
     SharedTensor2d G1c_ijA;
     SharedTensor2d G1c_ijB;
     SharedTensor2d G1c_abA;
     SharedTensor2d G1c_abB;
     SharedTensor2d G1c_ooA;
     SharedTensor2d G1c_ooB;
     SharedTensor2d G1c_ovA;
     SharedTensor2d G1c_ovB;
     SharedTensor2d G1c_voA;
     SharedTensor2d G1c_voB;
     SharedTensor2d G1c_vvA;
     SharedTensor2d G1c_vvB;
     SharedTensor2d G1cA;              // Correlation part of OPDM
     SharedTensor2d G1cB;              // Correlation part of OPDM
     SharedTensor2d G1A;                // Full OPDM
     SharedTensor2d G1B;                // Full OPDM
     SharedTensor2d GijA;
     SharedTensor2d GijB;
     SharedTensor2d GabA;
     SharedTensor2d GabB;
     SharedTensor2d GiaA;
     SharedTensor2d GiaB;
     SharedTensor2d GaiA;
     SharedTensor2d GaiB;
     SharedTensor2d GtijA;
     SharedTensor2d GtijB;
     SharedTensor2d GtabA;
     SharedTensor2d GtabB;

     // DF TPDM
     SharedTensor2d G2c_ij;
     SharedTensor2d G2c_ia;
     SharedTensor2d G2c_ab;
     SharedTensor2d G2c_oo;
     SharedTensor2d G2c_ov;
     SharedTensor2d G2c_vo;
     SharedTensor2d G2c_vv;
     SharedTensor2d G2c;               // Correlation part of TPDM
     SharedTensor2d G2c_ijA;
     SharedTensor2d G2c_ijB;
     SharedTensor2d G2c_iaA;
     SharedTensor2d G2c_iaB;
     SharedTensor2d G2c_abA;
     SharedTensor2d G2c_abB;
     SharedTensor2d G2c_ooA;
     SharedTensor2d G2c_ooB;
     SharedTensor2d G2c_ovA;
     SharedTensor2d G2c_ovB;
     SharedTensor2d G2c_voA;
     SharedTensor2d G2c_voB;
     SharedTensor2d G2c_vvA;
     SharedTensor2d G2c_vvB;
     SharedTensor2d G2cA;              // Correlation part of TPDM
     SharedTensor2d G2cB;              // Correlation part of TPDM
     SharedTensor1d Jc;                // Correlation Coulomb matrix
     SharedTensor1d g1Q;
     SharedTensor1d g1Qc;
     SharedTensor1d g1Qp;
     SharedTensor1d g1Qt;
     SharedTensor1d g1Qt2;

     // DF GFM
     SharedTensor2d GF;                // Full GFM (MO)
     SharedTensor2d GFao;              // Full GFM (AO)
     SharedTensor2d GFA;               // Full GFM
     SharedTensor2d GFB;               // Full GFM
     SharedTensor2d GFoo;
     SharedTensor2d GFov;
     SharedTensor2d GFvo;
     SharedTensor2d GFvv;
     SharedTensor2d GFooA;
     SharedTensor2d GFooB;
     SharedTensor2d GFovA;
     SharedTensor2d GFovB;
     SharedTensor2d GFvoA;
     SharedTensor2d GFvoB;
     SharedTensor2d GFvvA;
     SharedTensor2d GFvvB;

     SharedTensor2d GFtvv;            // Complement of GFM
     SharedTensor2d GFtvvA;           // Complement of GFM
     SharedTensor2d GFtvvB;           // Complement of GFM

     // MO gradient and Hessian
     SharedTensor2d Worb;              // MO gradient matrix
     SharedTensor2d WorbA;
     SharedTensor2d WorbB;
     SharedTensor2d Aorb;              // MO Hessian matrix
     SharedTensor2d AorbAA;
     SharedTensor2d AorbAB;
     SharedTensor2d AorbBB;
     SharedTensor2d Aaibj;
     SharedTensor2d Aaijk;
     SharedTensor2d Aijkl;
     SharedTensor2d Aoo;
     SharedTensor2d Avo;
     SharedTensor2d AooA;
     SharedTensor2d AooB;
     SharedTensor2d AvoA;
     SharedTensor2d AvoB;
     SharedTensor2d ZvoA;            // Zvector in matrix form
     SharedTensor2d ZvoB;            // Zvector in matrix form
     SharedTensor2d ZovA;            // Transpose of Zvector in matrix form
     SharedTensor2d ZovB;            // Transpose of Zvector in matrix form
     SharedTensor2d ZklA;            // AOCC-FC Zvector in matrix form
     SharedTensor2d ZklB;            // AOCC-FC Zvector in matrix form
     SharedTensor2d ZlkA;            // FC-AOCC Zvector in matrix form
     SharedTensor2d ZlkB;            // FC-AOCC Zvector in matrix form
     SharedTensor2d WvoA;            // Effective MO gradient VO block
     SharedTensor2d WvoB;            // Effective MO gradient VO block

     // Orbital rotations
     SharedTensor2d UorbA;           // MO rotation matrix: wrt reference MOs
     SharedTensor2d UorbB;
     SharedTensor2d KorbA;           // K matrix: wrt reference MOs
     SharedTensor2d KorbB;
     SharedTensor2d KsqrA;
     SharedTensor2d KsqrB;

     // MO rotation vectors
     SharedTensor1d wog;             // MO gradient vector
     SharedTensor1d wogA;            // MO gradient vector
     SharedTensor1d wogB;
     SharedTensor1d wog_intA;        // Interpolated MO gradient vector
     SharedTensor1d wog_intB;
     SharedTensor1d kappa;           // where kappa = kappaA + kappaB
     SharedTensor1d kappaA;          // vector of orb rot parameters: wrt previous MOS
     SharedTensor1d kappaB;
     SharedTensor1d kappa_barA;      // vector of orb rot parameters: wrt reference MOS
     SharedTensor1d kappa_barB;
     SharedTensor1d kappa_newA;
     SharedTensor1d kappa_newB;
     SharedTensor1d zvector;
     SharedTensor1d zvectorA;
     SharedTensor1d zvectorB;
     SharedTensor1d zvec_newA;
     SharedTensor1d zvec_newB;
     SharedTensor1d zvec_new;
     SharedTensor1d Wvo_vecA;            // Effective MO gradient vector VO block

     // PCG intermediates
     SharedTensor1d r_pcgA;
     SharedTensor1d r_pcgB;
     SharedTensor1d r_pcg;
     SharedTensor1d z_pcgA;
     SharedTensor1d z_pcgB;
     SharedTensor1d z_pcg;
     SharedTensor1d p_pcgA;
     SharedTensor1d p_pcgB;
     SharedTensor1d p_pcg;
     SharedTensor1d sigma_pcgA;
     SharedTensor1d sigma_pcgB;
     SharedTensor1d sigma_pcg;
     SharedTensor1d Minv_pcgA;
     SharedTensor1d Minv_pcgB;
     SharedTensor1d Minv_pcg;
     SharedTensor1d r_pcg_newA;
     SharedTensor1d r_pcg_newB;
     SharedTensor1d r_pcg_new;
     SharedTensor1d z_pcg_newA;
     SharedTensor1d z_pcg_newB;
     SharedTensor1d z_pcg_new;
     SharedTensor1d p_pcg_newA;
     SharedTensor1d p_pcg_newB;
     SharedTensor1d p_pcg_new;
     SharedTensor1d dr_pcgA;
     SharedTensor1d dr_pcgB;
     SharedTensor1d dr_pcg;
     SharedTensor1d residualA;
     SharedTensor1d residualB;
     SharedTensor1d residual;

     // Independent pairs
     SharedTensor1i idprowA;
     SharedTensor1i idprowB;
     SharedTensor1i idpcolA;
     SharedTensor1i idpcolB;

     // Diis
     SharedTensor2d vecsA;
     SharedTensor2d vecsB;
     SharedTensor2d errvecsA;
     SharedTensor2d errvecsB;

     // SO basis
     SharedTensor2d gQso;              // Gamma(Q|mu nu): 3-index TPDM
     SharedTensor2d gQso_ref;          // Gamma(Q|mu nu): 3-index TPDM
     SharedTensor2d gQoo;              // Gamma(Q|i i): 3-index TPDM
     SharedTensor2d gQoo_ref;          // Gamma(Q|i i): 3-index TPDM
     SharedTensor2d gQon_ref;          // Gamma(Q|i nu): 3-index TPDM
     SharedTensor2d Gaux;              // Gamma(P,Q): 2-index TPDM
     SharedTensor2d Gaux_ref;          // Gamma(P,Q): 2-index TPDM
     SharedTensor2d dQso;              // D(Q|mu nu): 3-index OPDM for REF WFN

     // Amplitudes
     SharedTensor2d t2_1;              // T_ij^ab(1)
     SharedTensor2d t2p_1;             // T'(ia,jb) = T_ij^ab(1)
     SharedTensor2d u2_1;              // 2*T_ij^ab(1) - T_ji^ab(1)
     SharedTensor2d u2p_1;             // U'(ia,jb) = 2*T_ij^ab(1) - T_ji^ab(1)
     SharedTensor2d t2p_1new;          // T'(ia,jb) = T_ij^ab(1)
     SharedTensor2d t2;                // T_ij^ab
     SharedTensor2d t2new;             // T_ij^ab
     SharedTensor2d l2;                // L_ij^ab
     SharedTensor2d l2new;             // L_ij^ab

     SharedTensor2d t2_1AA;            // T_ij^ab(1)
     SharedTensor2d t2_1AB;            // T_ij^ab(1)
     SharedTensor2d t2_1BB;            // T_ij^ab(1)
     SharedTensor2d t2_1newAA;         // T_ij^ab(1)
     SharedTensor2d t2_1newAB;         // T_ij^ab(1)
     SharedTensor2d t2_1newBB;         // T_ij^ab(1)

     SharedTensor2d t1A;               // T_i^a(1)
     SharedTensor2d t1B;               // T_i^a(1)
     SharedTensor2d t1newA;            // T_i^a(1)
     SharedTensor2d t1newB;            // T_i^a(1)
     SharedTensor1d T1c;               // T1_Q
     SharedTensor1d L1c;               // L1_Q
     SharedTensor2d l1A;               // T_i^a(1)
     SharedTensor2d l1B;               // T_i^a(1)
     SharedTensor2d l1newA;            // T_i^a(1)
     SharedTensor2d l1newB;            // T_i^a(1)
     SharedTensor1d gQ;                // G_Q
     SharedTensor1d gQp;               // G_Q'
     SharedTensor1d gQt;               // Gt_Q

     SharedTensor2d FijA;
     SharedTensor2d FijB;
     SharedTensor2d FabA;
     SharedTensor2d FabB;
     SharedTensor2d FiaA;
     SharedTensor2d FiaB;
     SharedTensor2d FtijA;
     SharedTensor2d FtijB;
     SharedTensor2d FtabA;
     SharedTensor2d FtabB;

     // Intermediates
     SharedTensor2d uQia;

     // Conventional integrals for the DF_BASIS_CC
     SharedTensor2d JijklAA;             // (ij|kl) (active)
     SharedTensor2d JijklAB;             // (ij|kl) (active)
     SharedTensor2d JijklBB;             // (ij|kl) (active)
     SharedTensor2d JijkaAA;             // (ij|ka) (active)
     SharedTensor2d JijkaAB;             // (ij|ka) (active)
     SharedTensor2d JijkaBB;             // (ij|ka) (active)
     SharedTensor2d JijabAA;             // (ij|ab) (active)
     SharedTensor2d JijabAB;             // (ij|ab) (active)
     SharedTensor2d JijabBB;             // (ij|ab) (active)
     SharedTensor2d JiajbAA;             // (ia|jb) (active)
     SharedTensor2d JiajbAB;             // (ia|jb) (active)
     SharedTensor2d JiajbBB;             // (ia|jb) (active)
     SharedTensor2d JiajkAB;             // (ka|ij) (active)
     SharedTensor2d JabijAB;             // (ab|ij) (active)

     SharedTensor2d JooooAA;             // (ij|kl) (all)
     SharedTensor2d JooooAB;             // (ij|kl) (all)
     SharedTensor2d JooooBB;             // (ij|kl) (all)
     SharedTensor2d JooovAA;             // (ij|ka) (all)
     SharedTensor2d JooovAB;             // (ij|ka) (all)
     SharedTensor2d JooovBB;             // (ij|ka) (all)
     SharedTensor2d JoovvAA;             // (ij|ab) (all)
     SharedTensor2d JoovvAB;             // (ij|ab) (all)
     SharedTensor2d JoovvBB;             // (ij|ab) (all)
     SharedTensor2d JovovAA;             // (ia|jb) (all)
     SharedTensor2d JovovBB;             // (ia|jb) (all)
     SharedTensor2d JovovAB;             // (ia|jb) (all)
     SharedTensor2d JovooAB;             // (ka|ij) (all)
     SharedTensor2d JvvooAB;             // (ab|ij) (all)

     SharedTensor2d IijklAA;             // <ij|kl> (active)
     SharedTensor2d IijklAB;             // <ij|kl> (active)
     SharedTensor2d IijklBB;             // <ij|kl> (active)
     SharedTensor2d IijkaAA;             // <ij|ka> (active)
     SharedTensor2d IijkaAB;             // <ij|ka> (active)
     SharedTensor2d IijkaBB;             // <ij|ka> (active)
     SharedTensor2d IijabAA;             // <ij|ab> (active)
     SharedTensor2d IijabAB;             // <ij|ab> (active)
     SharedTensor2d IijabBB;             // <ij|ab> (active)
     SharedTensor2d IiajbAA;             // <ia|jb> (active)
     SharedTensor2d IiajbAB;             // <ia|jb> (active)
     SharedTensor2d IiajbBB;             // <ia|jb> (active)
     SharedTensor2d IijakAB;             // <ij|ak> (active)
     SharedTensor2d IaibjAB;             // <ai|bj> (active)

     SharedTensor2d IooooAA;             // <ij|kl> (all)
     SharedTensor2d IooooAB;             // <ij|kl> (all)
     SharedTensor2d IooooBB;             // <ij|kl> (all)
     SharedTensor2d IooovAA;             // <ij|ka> (all)
     SharedTensor2d IooovAB;             // <ij|ka> (all)
     SharedTensor2d IooovBB;             // <ij|ka> (all)
     SharedTensor2d IoovvAA;             // <ij|ab> (all)
     SharedTensor2d IoovvAB;             // <ij|ab> (all)
     SharedTensor2d IoovvBB;             // <ij|ab> (all)
     SharedTensor2d IovovAA;             // <ia|jb> (all)
     SharedTensor2d IovovAB;             // <ia|jb> (all)
     SharedTensor2d IovovBB;             // <ia|jb> (all)
     SharedTensor2d IoovoAB;             // <ij|ak> (all)
     SharedTensor2d IvovoAB;             // <ai|bj> (all)

     SharedTensor2d AIijklAA;             // <ij||kl> (active)
     SharedTensor2d AIijklBB;             // <ij||kl> (active)
     SharedTensor2d AIijkaAA;             // <ij||ka> (active)
     SharedTensor2d AIijkaBB;             // <ij||ka> (active)
     SharedTensor2d AIijabAA;             // <ij||ab> (active)
     SharedTensor2d AIijabBB;             // <ij||ab> (active)
     SharedTensor2d AIiajbAA;             // <ia||jb> (active)
     SharedTensor2d AIiajbBB;             // <ia||jb> (active)

     SharedTensor2d AIooooAA;             // <ij||kl> (all)
     SharedTensor2d AIooooBB;             // <ij||kl> (all)
     SharedTensor2d AIooovAA;             // <ij||ka> (all)
     SharedTensor2d AIooovBB;             // <ij||ka> (all)
     SharedTensor2d AIoovvAA;             // <ij||ab> (all)
     SharedTensor2d AIoovvBB;             // <ij||ab> (all)
     SharedTensor2d AIovovAA;             // <ia||jb> (all)
     SharedTensor2d AIovovBB;             // <ia||jb> (all)

     // 1D-Tensors
     SharedTensor1d DQvecA;
     SharedTensor1d dQ;

     // Pair indices
     SharedTensor2i so_idx;             // Pair index for active SO
     SharedTensor2i ij_idxAA;           // Pair index for active OO
     SharedTensor2i ij_idxAB;           // Pair index for active OO
     SharedTensor2i ij_idxBA;           // Pair index for active OO
     SharedTensor2i ij_idxBB;           // Pair index for active OO
     SharedTensor2i ia_idxAA;           // Pair index for active OV
     SharedTensor2i ia_idxAB;           // Pair index for active OV
     SharedTensor2i ia_idxBA;           // Pair index for active OV
     SharedTensor2i ia_idxBB;           // Pair index for active OV
     SharedTensor2i ai_idxAA;           // Pair index for active VO
     SharedTensor2i ai_idxAB;           // Pair index for active VO
     SharedTensor2i ai_idxBA;           // Pair index for active VO
     SharedTensor2i ai_idxBB;           // Pair index for active VO
     SharedTensor2i ab_idxAA;           // Pair index for active VV
     SharedTensor2i ab_idxAB;           // Pair index for active VV
     SharedTensor2i ab_idxBA;           // Pair index for active VV
     SharedTensor2i ab_idxBB;           // Pair index for active VV

     SharedTensor2i oo_idxAA;           // Pair index for all OO
     SharedTensor2i oo_idxAB;           // Pair index for all OO
     SharedTensor2i oo_idxBB;           // Pair index for all OO
     SharedTensor2i ov_idxAA;           // Pair index for all OV
     SharedTensor2i ov_idxAB;           // Pair index for all OV
     SharedTensor2i ov_idxBB;           // Pair index for all OV
     SharedTensor2i vo_idxAA;           // Pair index for all VO
     SharedTensor2i vo_idxAB;           // Pair index for all VO
     SharedTensor2i vo_idxBB;           // Pair index for all VO
     SharedTensor2i vv_idxAA;           // Pair index for all VV
     SharedTensor2i vv_idxAB;           // Pair index for all VV
     SharedTensor2i vv_idxBB;           // Pair index for all VV

     SharedMatrix Tso_;
     SharedMatrix Vso_;
     SharedMatrix Hso_;
     SharedMatrix Sso_;
     SharedMatrix bQnn;         // b(Q|mu nu)
     SharedVector e_orbA;

};

} }

#endif // dfocc_h

