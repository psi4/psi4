set(headers_list "")
# List of headers
list(APPEND headers_list occwave.h defines.h dpd.h arrays.h )

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
list(APPEND sources_list mograd.cc t1_1st.cc trans_ints_ump2.cc t2_2nd_general.cc update_mo.cc main.cc fock_alpha.cc tpdm_ref_corr_opdm.cc dpd.cc trans_ints_uhf.cc ocepa_g_int.cc ekt_ip.cc omp3_t2_1st_sc.cc t2_2nd_sc.cc occwave.cc ocepa_t2_1st_sc.cc idp.cc v_int.cc omp2_ip_poles.cc kappa_orb_resp.cc tei_sort_iabc.cc ocepa_response_pdms.cc coord_grad.cc omp2_response_pdms.cc diis.cc fock_beta.cc ep2_ip.cc trans_ints_rmp2.cc idp2.cc z_vector.cc arrays.cc occ_iterations.cc kappa_msd.cc omp3_ip_poles.cc manager.cc get_moinfo.cc ccl_energy.cc semi_canonic.cc gfock.cc ekt_ea.cc gfock_ea.cc kappa_orb_resp_iter.cc trans_ints_rhf.cc w_int.cc cc_energy.cc t2_amps.cc w_1st_order.cc corr_tpdm.cc gfock_diag.cc omp2_t2_1st.cc cepa_iterations.cc v_2nd_order.cc omp3_t2_1st_general.cc omp3_g_int.cc omp3_response_pdms.cc )

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
add_library(occ STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(occ mints)
set(libs_to_merge occ ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(occ custom_boost)
endif()
install(TARGETS occ ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libocc)
