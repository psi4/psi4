set(sources_list response.cc jk_grad.cc scf_grad.cc wrapper.cc )
psi4_add_module(bin scfgrad sources_list mints 3index)
