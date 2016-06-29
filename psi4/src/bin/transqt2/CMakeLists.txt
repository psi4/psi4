set(sources_list transtwo_rhf.cc cache.cc get_params.cc transtwo_uhf.cc transqt.cc transone.cc get_moinfo.cc file_build_presort.cc idx_error.cc idx_permute_presort.cc )
psi4_add_binary(transqt2 sources_list)
target_link_libraries(transqt2 mints)
