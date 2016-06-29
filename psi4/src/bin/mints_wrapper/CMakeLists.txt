set(sources_list mints.cc )
psi4_add_binary(mints_wrapper sources_list)
link_libraries(mints_wrapper mints)