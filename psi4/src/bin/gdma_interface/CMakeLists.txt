set(sources_list wrapper.cc)
psi4_add_binary(gdma_interface sources_list mints)
target_include_directories(gdma_interface PRIVATE ${GDMA_INCLUDE_DIRS})
