set(sources_list transtwo_rhf.cc cache.cc get_params.cc transtwo_uhf.cc transqt.cc transone.cc get_moinfo.cc file_build_presort.cc idx_error.cc idx_permute_presort.cc )
psi4_add_binary(transqt2 sources_list mints)
target_include_directories(transqt2 PUBLIC ${Boost_INCLUDE_DIRS})
