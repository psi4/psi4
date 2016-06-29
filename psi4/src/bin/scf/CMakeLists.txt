set(sources_list main.cc )
psi4_add_binary(scf sources_list)
target_link_libraries(scf mints)
