set(headers_list "")
# List of headers
list(APPEND headers_list git.h block_matrix.h sblock_vector.h scf.h sblock_matrix.h block_vector.h algebra_interface_mangle.h matrix_base.h vector_base.h mcscf.h algebra_interface.h )

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
list(APPEND sources_list mcscf.cc block_vector.cc sblock_vector.cc scf_energy.cc block_matrix.cc scf_print_eigenvectors_and_MO.cc scf_diis.cc scf_compute_energy.cc scf_Feff.cc algebra_interface.cc scf_check_orthonormality.cc scf_guess_occupation.cc scf_save_info.cc vector_base.cc scf.cc scf_G.cc scf_F.cc scf_canonicalize_MO.cc scf_initial_guess.cc sblock_matrix.cc scf_density_matrix.cc scf_read_so_tei.cc scf_iterate_scf_equations.cc scf_read_so_oei.cc matrix_base.cc scf_S_inverse_sqrt.cc scf_pairs.cc )

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
add_library(mcscf STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(mcscf mints)
set(libs_to_merge mcscf ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(mcscf custom_boost)
endif()
install(TARGETS mcscf ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libmcscf)
