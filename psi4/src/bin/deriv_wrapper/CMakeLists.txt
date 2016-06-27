set(sources_list main.cc )
psi4_add_binary(deriv_wrapper sources_list)
target_link_libraries(deriv_wrapper mints)