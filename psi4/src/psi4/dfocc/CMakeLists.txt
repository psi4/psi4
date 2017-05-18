set(sources_list
approx_diag_mohess_oo.cc  fock.cc                   mp2_direct.cc
approx_diag_mohess_vo.cc  fock_so.cc                occ_iterations.cc
get_moinfo.cc
cc_energy.cc              gfock_oo.cc               omp2_opdm.cc
ccl_energy.cc             gfock_ov.cc               omp2_tpdm.cc
cd_ints.cc                gfock_vo.cc               pair_index.cc
conv_mo_tei.cc            gfock_vv.cc               ref_grad.cc
conv_mo_tei_direct.cc     gftilde_vv.cc             semi_canonic.cc
idp.cc                    separable_tpdm.cc
conv_mo_tei_ref_direct.cc kappa_diag_hess.cc        t1_1st_sc.cc
df_corr.cc                kappa_orb_resp.cc         t2_1st_gen.cc
df_ref.cc                 kappa_orb_resp_pcg.cc     t2_1st_sc.cc
dfocc.cc                  main.cc
diagonal_mohess_oo.cc     manager.cc
diagonal_mohess_vo.cc     manager_cd.cc             tensors.cc
diis.cc                   mograd.cc                 update_mo.cc
tpdm_tilde.cc             back_trans.cc
oei_grad.cc               tei_grad_ref.cc           tei_grad_corr.cc
z_vector_pcg.cc           z_vector.cc               prepare4grad.cc
dfgrad.cc                 s2_response.cc
t2_mp2_direct.cc          effective_mograd.cc       qchf.cc
update_hfmo.cc            ccsd_iterations.cc        ccsd_3index_intr.cc
ccsd_F_intr.cc            ccsd_W_intr.cc            ccsd_t1_amps.cc
ccsd_t2_amps.cc           ccd_iterations.cc         ccd_3index_intr.cc
ccd_F_intr.cc             ccd_W_intr.cc             ccd_t2_amps.cc
ccsd_iterations_low.cc    ccsd_3index_intr_low.cc   ccsd_F_intr_low.cc
ccsd_W_intr_low.cc        ccsd_t1_amps_low.cc       ccsd_t2_amps_low.cc
ccd_iterations_low.cc     ccd_3index_intr_low.cc    ccd_F_intr_low.cc
ccd_W_intr_low.cc         ccd_t2_amps_low.cc        ccsdl_iterations.cc
ccsdl_3index_intr.cc      ccsdl_W_intr.cc           ccsdl_l1_amps.cc
ccsdl_l2_amps.cc          ccdl_iterations.cc        ccdl_3index_intr.cc
ccdl_W_intr.cc            ccdl_l2_amps.cc           ccsd_pdm_3index_intr.cc
ccsd_opdm.cc              ccsd_tpdm.cc              ccd_pdm_3index_intr.cc
ccd_opdm.cc               ccd_tpdm.cc               ccsd_triples.cc
mp3_W_intr.cc             t2_2nd_sc.cc              t2_2nd_gen.cc
omp3_opdm.cc              omp3_tpdm.cc              mp3_pdm_3index_intr.cc
lccd_iterations.cc        olccd_tpdm.cc             lccd_W_intr.cc
lccd_pdm_3index_intr.cc   lccd_t2_amps.cc
# arrays.cc olddf.cc t2_1st_scs_gen.cc t2_1st_scs_sc.cc z_vector_cg.cc
# combine_ref_sep_tpdm.cc conv_mo_tei_ref.cc
)
psi4_add_module(bin dfocc sources_list mints diis)
