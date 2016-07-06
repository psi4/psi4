set(sources_list mints.cc )
psi4_add_binary(mints_wrapper sources_list mints)
target_include_directories(mints_wrapper PUBLIC ${Boost_INCLUDE_DIRS})