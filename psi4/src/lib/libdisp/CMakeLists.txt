set(sources_list dispersion.cc )
psi4_add_library(disp sources_list)
target_link_libraries(disp mints)
