# List of headers (needed to avoid an "empty list error" from CMake)
set(headers_list "")

set(sources_list "")
# List of sources
list(APPEND sources_list erd__gener_eri_derv_batch.f erd__derv_int2d_to_e0cd.f erd__memory_eri_batch.f erd__e0f0_pcgto_block.f erd__cartesian_norms.f erd__hrr_transform.f erd__csgto.f erd__ssss_pcgto_block.f erd__sspp_pcgto_block.f erd__set_abcd.f erd__memory_eri_derv_batch.f erd__ctr_1st_half.f erd__1111_csgto.f erd__normalize_cartesian.f erd__ctr_4index_block.f erd__dsqmin_line_segments.f erd__rys_4_roots_weights.f erd__derv_int2d_to_0000.f erd__derv_def_blocks.f erd__ctr_2nd_half_new.f erd__rys_roots_weights.f erd__derv_int2d_to_e000.f erd__transpose_batch.f erd__print_batch.f erd__rys_1_roots_weights.f erd__derv_csgto.f erd__int2d_to_e000.f erd__set_derv_abcd.f erd__memory_1111_csgto.f erd__move_ry.f erd__hrr_matrix.f erd__ctr_4index_reorder.f erd__memory_hrr_correction.f erd__derv_int2d_to_a000.f erd__2d_coefficients.f erd__derv_int2d_to_abc0.f erd__sppp_pcgto_block.f erd__ctr_2nd_half_update.f erd__derv_pcgto_block.f erd__hrr_step.f erd__2d_atom_coefficients.f erd__gener_eri_batch.f erd__rys_3_roots_weights.f erd__2d_derv_pcd_integrals.f erd__xyz_to_ry_abcd.f erd__2d_derv_abcd_integrals.f erd__set_ij_kl_pairs.f erd__rys_2_roots_weights.f erd__derv_int2d_to_ab00.f erd__e0f0_def_blocks.f erd__2d_pq_integrals.f erd__atom_int2d_to_e0f0.f erd__derv_int2d_to_abcd.f erd__derv_int2d_to_e0c0.f erd__1111_def_blocks.f erd__set_derv_sequence.f erd__map_ijkl_to_ikjl.f erd__ctr_rs_expand.f erd__2d_abcd_integrals.f erd__xyz_to_ry_matrix.f erd__sssp_pcgto_block.f erd__atom_int2d_to_e000.f erd__2d_pcd_integrals.f erd__ctr_tu_expand.f erd__memory_derv_csgto.f erd__spherical_transform.f erd__memory_csgto.f erd__int2d_to_e0f0.f erd__rys_5_roots_weights.f erd__pppp_pcgto_block.f erd__2d_atom_pq_integrals.f erd__rys_x_roots_weights.f erd__prepare_ctr.F )

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
add_library(erd STATIC ${sources_list})
# Specify dependencies for the library (if any)
#add_dependencies(erd )
set(libs_to_merge erd ${libs_to_merge} PARENT_SCOPE)
install(TARGETS erd ARCHIVE DESTINATION lib)

