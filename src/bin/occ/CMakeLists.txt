set(SRC ekt_ea.cc               kappa_orb_resp_iter.cc  omp3_g_int.cc           trans_ints_rmp2.cc 
ekt_ip.cc               main.cc                 omp3_ip_poles.cc        trans_ints_uhf.cc 
arrays.cc               ep2_ip.cc               manager.cc              omp3_response_pdms.cc   trans_ints_ump2.cc 
fock_alpha.cc           mograd.cc               omp3_t2_1st_general.cc  update_mo.cc 
cc_energy.cc            fock_beta.cc            occ_iterations.cc       omp3_t2_1st_sc.cc       v_2nd_order.cc 
ccl_energy.cc           get_moinfo.cc           occwave.cc              semi_canonic.cc         v_int.cc 
cepa_iterations.cc      gfock.cc                t1_1st.cc               w_1st_order.cc 
coord_grad.cc           gfock_diag.cc           ocepa_g_int.cc          t2_2nd_general.cc       w_int.cc 
corr_tpdm.cc            gfock_ea.cc             ocepa_response_pdms.cc  t2_2nd_sc.cc            z_vector.cc 
idp.cc                  ocepa_t2_1st_sc.cc      t2_amps.cc              
diis.cc                 idp2.cc                 omp2_ip_poles.cc        tei_sort_iabc.cc        
dpd.cc                  kappa_msd.cc            omp2_response_pdms.cc   tpdm_ref_corr_opdm.cc   
kappa_orb_resp.cc       omp2_t2_1st.cc          trans_ints_rhf.cc)
add_library(occ ${SRC})
add_dependencies(occ mints)
