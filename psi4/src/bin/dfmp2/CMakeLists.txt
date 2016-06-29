set(sources_list mp2.cc corr_grad.cc wrapper.cc )
psi4_add_binary(dfmp2 sources_list)
target_link_libraries(dfmp2 mints)
