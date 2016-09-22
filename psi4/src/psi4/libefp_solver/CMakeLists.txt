set(sources_list efp_solver.cc)
psi4_add_module(lib efp_solver sources_list mints filesystem)
target_link_libraries(efp_solver PUBLIC libefp::efp)
