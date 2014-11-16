set(headers_list "")
# List of headers
list(APPEND headers_list MOInfo.h globals.h Params.h Frozen.h )

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
list(APPEND sources_list sortI.cc build_A_ROHF.cc sort_ltd_uhf.cc ex_sort_td_uhf.cc zero_pdm.cc get_rho_params.cc build_Z_ROHF.cc x_Gciab.cc ex_td_cleanup.cc classify.cc relax_I_ROHF.cc build_A_RHF.cc relax_D.cc get_frozen.cc cache.cc Gabcd.cc build_A_UHF.cc dump_RHF.cc Iij.cc build_X.cc td_cleanup.cc ex_tdensity_intermediates.cc transL.cc build_ex_tdensity.cc Iia.cc onepdm.cc Gciab.cc x_xi1_uhf.cc x_te_intermediates.cc relax_I_RHF.cc ex_tdensity.cc sort_rtd_uhf.cc x_xi2.cc sort_ltd_rohf.cc x_te_intermediates_rhf.cc x_xi1_connected.cc get_params.cc x_xi_check.cc energy_RHF.cc ex_sort_td_rohf.cc rotational_strength.cc add_ref_UHF.cc G.cc x_Gijkl.cc G_norm.cc x_xi2_rhf.cc resort_tei.cc sortone.cc oscillator_strength.cc Gijkl.cc resort_gamma.cc x_onepdm.cc Iai.cc x_Gijka_uhf.cc x_Gibja.cc Gijka.cc get_td_params.cc ex_rotational_strength.cc add_core_ROHF.cc distribute.cc sortone_ROHF.cc relax_I.cc relax_I_UHF.cc ex_td_setup.cc Gijab_ROHF.cc idx_permute.cc deanti_ROHF.cc dump_ROHF.cc sort_rtd_rohf.cc build_Z_UHF.cc transp.cc dipole.cc get_moinfo.cc twopdm.cc file_build.cc sortI_RHF.cc Gijab_RHF.cc x_Gijab.cc Gibja.cc ccdensity.cc tdensity.cc x_xi2_uhf.cc add_ref.cc deanti_RHF.cc ael.cc sortI_UHF.cc x_oe_intermediates.cc fold_UHF.cc td_setup.cc energy_UHF.cc energy_ROHF.cc energy.cc x_Gibja_uhf.cc x_xi1_rhf.cc rtdensity.cc x_xi1.cc x_Gijka.cc ltdensity_intermediates.cc add_ref_RHF.cc sortI_ROHF.cc fold_ROHF.cc V.cc idx_error.cc x_V.cc Gijab_UHF.cc td_print.cc lag.cc kinetic.cc build_Z.cc Iab.cc x_Gabcd.cc add_core_UHF.cc norm.cc ex_oscillator_strength.cc x_Gijab_uhf.cc ltdensity.cc fold.cc x_oe_intermediates_rhf.cc x_onepdm_uhf.cc build_Z_RHF.cc deanti.cc densgrid_RHF.cc setup_LR.cc fold_RHF.cc add_ref_ROHF.cc dump_UHF.cc ex_td_print.cc deanti_UHF.cc build_A.cc transdip.cc sortone_RHF.cc sortone_UHF.cc x_Gciab_uhf.cc Gijab.cc x_xi_intermediates.cc )

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
add_library(ccdensity STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(ccdensity mints)
set(libs_to_merge ccdensity ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(ccdensity custom_boost)
endif()
install(TARGETS ccdensity ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libccdensity)
