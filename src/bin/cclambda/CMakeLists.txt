set(headers_list "")
# List of headers
list(APPEND headers_list MOInfo.h globals.h Params.h Local.h cclambda.h )

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
list(APPEND sources_list cc3_t3z.cc local.cc WejabL2.cc status.cc check_sum.cc hbar_extra.cc WefabL2.cc cc2_Gai.cc Lmag.cc cc2_fmiL2.cc cache.cc spinad_amps.cc FmiL2.cc Lsave.cc WabeiL1.cc L3_AAB.cc overlap.cc cc2_L2.cc Lnorm.cc get_params.cc check_ortho.cc L3_AAA.cc G.cc converged.cc sort_amps.cc L1FL2.cc projections.cc overlap_LAMPS.cc ortho_Rs.cc cclambda.cc cc3_l3l1.cc diis.cc halftrans.cc init_amps.cc cc3_l3l2.cc L1.cc get_moinfo.cc cc2_hbar_extra.cc update.cc GL2.cc WijmbL2.cc FaeL2.cc c_clean.cc pseudoenergy.cc BL2_AO.cc WmbejL2.cc cc2_faeL2.cc cc2_L1.cc DL2.cc cc3_t3x.cc denom.cc Lamp_write.cc dijabL2.cc WijmnL2.cc L2.cc )

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
add_library(cclambda STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(cclambda mints)
set(libs_to_merge cclambda ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(cclambda custom_boost)
endif()
install(TARGETS cclambda ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libcclambda)
