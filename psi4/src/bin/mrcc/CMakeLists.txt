set(sources_list mrcc.cc )
psi4_add_binary(mrcc sources_list)
target_link_libraries(mrcc mints)