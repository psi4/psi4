find_package(libefp REQUIRED)
set(sources_list efp_solver.cc
)
psi4_add_library(efp_solver sources_list mints ${LIBEFP_LIBRARIES})
target_include_directories(efp_solver PUBLIC ${LIBEFP_INCLUDE_DIRS})
