set(sources_list response.cc jk_grad.cc scf_grad.cc wrapper.cc )
psi4_add_binary(scfgrad sources_list)
target_link_libraries(scfgrad mints)

