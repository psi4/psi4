set(sources_list efp_solver.cc)
psi4_add_library(efp_solver sources_list mints)
target_link_libraries(efp_solver PUBLIC libefp::efp)

