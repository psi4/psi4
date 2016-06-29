set(sources_list diismanager.cc diisentry.cc )
psi4_add_library(diis sources_list)
target_link_libraries(diis mints)