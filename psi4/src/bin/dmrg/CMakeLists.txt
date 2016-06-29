set(sources_list dmrgscf.cc)
psi4_add_binary(dmrg sources_list)
target_link_libraries(dmrg mints CHEMPS2::CHEMPS2)
