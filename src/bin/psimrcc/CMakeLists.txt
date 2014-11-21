set(headers_list "")
# List of headers
list(APPEND headers_list git.h main.h updater.h integral_class.h debugging.h heff.h index_types.h special_matrices.h operation.h index_iterator.h index.h mrcc.h blas.h matrix_types.h types.h matrix_iterator.h sort.h algebra_interface_mangle.h idmrpt2.h matrixtmp.h manybody.h algebra_interface.h transform.h mrccsd_t.h mp2_ccsd.h psimrcc.h matrix.h )

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
list(APPEND sources_list mrccsd_t_compute_restricted.cc mrcc_f_int.cc idmrpt2_f_int.cc mp2_ccsd_t1_amps.cc idmrpt2_Heff.cc sort_mrpt2.cc idmrpt2_Heff_doubles.cc mrccsd_t_heff_a.cc mrccsd_t_setup.cc idmrpt2_add_matrices.cc blas_resorting.cc transform.cc main.cc mrccsd_t_heff_a_restricted.cc blas_interface.cc sort.cc algebra_interface.cc mrcc_add_matrices.cc blas.cc mrcc_t1_amps.cc transform_block.cc blas_solve.cc blas_diis.cc mrcc_pert_cbs.cc mrccsd_t_heff_ab_restricted.cc idmrpt2_t2_amps.cc debugging.cc matrixtmp.cc mrcc_Heff.cc transform_read_so.cc matrix_iterator.cc psimrcc.cc mrccsd_t_heff_b_restricted.cc mrcc_t2_amps.cc matrix_memory_and_io.cc mrccsd_t_heff_ab.cc mrcc_tau.cc mrccsd_t_compute_spin_adapted.cc manybody.cc blas_parser.cc mrcc_z_int.cc mrccsd_t_compute.cc index_iterator.cc mrcc_pert_triples.cc manybody_denominators.cc mrcc_energy.cc mrccsd_t_heff.cc mp2_ccsd.cc mrcc_t_amps.cc updater_bw.cc blas_algorithms.cc updater_sr.cc mp2_ccsd_amps.cc mrccsd_t_heff_restricted.cc operation_contraction.cc heff_diagonalize.cc operation.cc special_matrices.cc mrccsd_t_form_matrices.cc operation_sort.cc mrccsd_t_heff_b.cc mp2_ccsd_f_int.cc matrix_addressing.cc idmrpt2.cc mrcc_mkccsd.cc mp2_ccsd_add_matrices.cc sort_out_of_core.cc operation_compute.cc sort_in_core.cc mrccsd_t.cc matrix.cc mp2_ccsd_z_int.cc transform_presort.cc updater.cc mrcc_w_int.cc mrcc_compute.cc heff.cc idmrpt2_t1_amps.cc updater_mk.cc mp2_ccsd_t2_amps.cc idmrpt2_Heff_singles.cc blas_compatibile.cc mp2_ccsd_w_int.cc mrcc.cc mrcc_mp3.cc transform_mrpt2.cc index.cc )

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
add_library(psimrcc STATIC ${sources_list})
# Specify dependencies for the library (if any)
#add_dependencies(psimrcc )
set(libs_to_merge psimrcc ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(psimrcc custom_boost)
endif()
install(TARGETS psimrcc ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libpsimrcc)
