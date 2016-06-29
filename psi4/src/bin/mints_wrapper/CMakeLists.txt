set(sources_list mints.cc )
psi4_add_binary(mints_wrapper sources_list)
target_include_directories(mints_wrapper PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(mints_wrapper mints)