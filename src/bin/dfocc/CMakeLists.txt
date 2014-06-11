set(SRC
approx_diag_mohess_oo.cc  fock.cc                   mp2_direct.cc
approx_diag_mohess_vo.cc  fock_so.cc                occ_iterations.cc
arrays.cc                 get_moinfo.cc             olddf.cc
cc_energy.cc              gfock_oo.cc               omp2_opdm.cc
ccl_energy.cc             gfock_ov.cc               omp2_tpdm.cc
cd_ints.cc                gfock_vo.cc               pair_index.cc
conv_mo_tei.cc            gfock_vv.cc               ref_grad.cc
conv_mo_tei_direct.cc     gftilde_vv.cc             semi_canonic.cc
conv_mo_tei_ref.cc        idp.cc                    separable_tpdm.cc
conv_mo_tei_ref_direct.cc kappa_diag_hess.cc        t1_1st_sc.cc
df_corr.cc                kappa_orb_resp.cc         t2_1st_gen.cc
df_ref.cc                 kappa_orb_resp_pcg.cc     t2_1st_sc.cc
dfocc.cc                  main.cc                   t2_1st_scs_gen.cc
diagonal_mohess_oo.cc     manager.cc                t2_1st_scs_sc.cc
diagonal_mohess_vo.cc     manager_cd.cc             tensors.cc
diis.cc                   mograd.cc                 update_mo.cc
combine_ref_sep_tpdm.cc   tpdm_tilde.cc             back_trans.cc
oei_grad.cc               tei_grad_ref.cc           tei_grad_corr.cc
z_vector_pcg.cc           z_vector.cc               prepare4grad.cc
)
add_library(dfocc ${SRC})
add_dependencies(dfocc mints)
