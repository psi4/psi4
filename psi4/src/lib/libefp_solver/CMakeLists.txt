find_package(libefp REQUIRED)
set(sources_list efp_solver.cc
)
psi4_add_library(efp_solver sources_list)
target_include_directories(efp_solver PUBLIC efp)
target_link_libraries(efp_solver mints efp)
