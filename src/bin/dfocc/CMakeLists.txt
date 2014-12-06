set(headers_list "")
# List of headers
list(APPEND headers_list defines.h dfocc.h tensors.h arrays.h )

# If you want to remove some headers specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM headers_list "")
else()
   list(REMOVE_ITEM headers_list "")
endif()
# Sort alphabetically
list(SORT headers_list)

set(sources_list "")
# List of sources
list(APPEND sources_list 
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
dfgrad.cc                 s2_response.cc            z_vector_cg.cc
t2_mp2_direct.cc          effective_mograd.cc       qchf.cc
update_hfmo.cc            ccsd_iterations.cc        ccsd_3index_intr.cc
ccsd_F_intr.cc            ccsd_W_intr.cc            ccsd_t1_amps.cc
ccsd_t2_amps.cc           ccd_iterations.cc         ccd_3index_intr.cc
ccd_F_intr.cc             ccd_W_intr.cc             ccd_t2_amps.cc
ccsd_iterations_low.cc    ccsd_3index_intr_low.cc   ccsd_F_intr_low.cc
ccsd_W_intr_low.cc        ccsd_t1_amps_low.cc       ccsd_t2_amps_low.cc
ccd_iterations_low.cc     ccd_3index_intr_low.cc    ccd_F_intr_low.cc
ccd_W_intr_low.cc         ccd_t2_amps_low.cc
)

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "")
else()
   list(REMOVE_ITEM sources_list "")
endif()

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(dfocc STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(dfocc mints)
set(libs_to_merge dfocc ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(dfocc custom_boost)
endif()
install(TARGETS dfocc ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libdfocc)
