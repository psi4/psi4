set(sources_list construct_sigma.cc 
                 adc_main.cc 
                 adc.cc 
                 compute_energy.cc 
                 differentiation.cc 
                 diagonalize.cc 
                 amps_write.cc 
                 denominator.cc 
                 prepare_tensors.cc 
                 init_tensors.cc 
)
psi4_add_module(bin adc sources_list mints)
