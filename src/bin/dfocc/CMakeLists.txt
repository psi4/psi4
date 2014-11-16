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
list(APPEND sources_list mograd.cc dfgrad.cc mp2_direct.cc combine_ref_sep_tpdm.cc z_vector_pcg.cc gfock_vv.cc update_mo.cc main.cc kappa_diag_hess.cc tensors.cc pair_index.cc t2_1st_scs_sc.cc tei_grad_ref.cc z_vector_solver.cc tpdm_tilde.cc manager_cd.cc fock.cc conv_mo_tei_ref_direct.cc fock_so.cc omp2_tpdm.cc approx_diag_mohess_vo.cc olddf.cc t2_1st_scs_gen.cc idp.cc ref_grad.cc df_ref.cc kappa_orb_resp.cc conv_mo_tei_ref.cc diagonal_mohess_oo.cc conv_mo_tei.cc diis.cc effective_mograd.cc z_vector.cc s2_response.cc gfock_oo.cc arrays.cc occ_iterations.cc manager.cc t2_1st_gen.cc update_hfmo.cc approx_diag_mohess_oo.cc t1_1st_sc.cc get_moinfo.cc ccl_energy.cc semi_canonic.cc gfock_ov.cc prepare4grad.cc tei_grad_corr.cc dfocc.cc cd_ints.cc diagonal_mohess_vo.cc separable_tpdm.cc omp2_opdm.cc cc_energy.cc t2_mp2_direct.cc gftilde_vv.cc back_trans.cc oei_grad.cc z_vector_cg.cc conv_mo_tei_direct.cc qchf.cc df_corr.cc gfock_vo.cc t2_1st_sc.cc kappa_orb_resp_pcg.cc )

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
