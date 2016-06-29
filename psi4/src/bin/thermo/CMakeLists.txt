set(sources_list thermo.cc )
psi4_add_binary(thermo sources_list)
target_link_libraries(thermo mints)
