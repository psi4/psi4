set(headers_list "")
# List of headers
list(APPEND headers_list print.h stre.h mem.h opt_params.h cov_radii.h frag.h IRC_data.h globals.h interfrag.h coordinates.h package.h molecule.h bend.h cart.h combo_coordinates.h opt_except.h io.h simple_base.h linear_algebra.h tors.h atom_data.h v3d.h physconst.h opt_data.h fb_frag.h )

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
list(APPEND sources_list molecule_sd_step.cc geom_gradients_io.cc lindh_guess.cc frag.cc print.cc getIntcoFileName.cc cart.cc linear_algebra.cc interfrag_orient.cc frag_apply_frozen_constraints.cc fb_frag.cc molecule_fragments.cc molecule_print.cc molecule_linesearch_step.cc mem.cc bend.cc frag_print.cc optking.cc molecule_read_coords.cc combo_coordinates.cc frag_H_guess.cc stre.cc molecule_irc_step.cc molecule_prfo_step.cc tors.cc opt_data_io.cc v3d.cc molecule_nr_step.cc molecule_rfo_step.cc molecule_tests.cc atom_data.cc molecule.cc set_params.cc molecule_backstep.cc opt_data.cc interfrag.cc frag_disp.cc )

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
add_library(optking STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(optking mints)
set(libs_to_merge optking ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(optking custom_boost)
endif()
install(TARGETS optking ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/liboptking)
