set(sources_list csg.cc cubeprop.cc)
psi4_add_library(cubeprop sources_list)
target_link_libraries(cubeprop mints)

