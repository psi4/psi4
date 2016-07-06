set(sources_list response.cc jk_grad.cc scf_grad.cc wrapper.cc )
psi4_add_binary(scfgrad sources_list mints 3index)

