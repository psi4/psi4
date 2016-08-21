set(sources_list dmrgscf.cc)
psi4_add_binary(dmrg sources_list mints)
target_include_directories(dmrg PRIVATE ${CHEMPS2_INCLUDE_DIRS})
target_link_libraries(dmrg PUBLIC ${CHEMPS2_LIBRARIES})
