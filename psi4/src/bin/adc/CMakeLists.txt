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
psi4_add_binary(adc sources_list)
target_link_libraries(adc mints)

