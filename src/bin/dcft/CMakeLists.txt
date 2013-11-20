set(SRC AO_contribute.cc dcft.cc dcft_compute.cc dcft_density.cc dcft_energy.cc dcft_gradient.cc dcft_integrals.cc dcft_intermediates.cc dcft_lambda.cc dcft_memory.cc dcft_mp2.cc dcft_n_representability.cc dcft_oo.cc dcft_qc.cc dcft_scf.cc dcft_tau.cc half_transform.cc dcft_triples.cc main.cc)
add_library(dcft ${SRC})
add_dependencies(dcft mints)
