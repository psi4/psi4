list(APPEND sources
  adc.cc
  adc_main.cc
  amps_write.cc
  compute_energy.cc
  construct_sigma.cc
  denominator.cc
  diagonalize.cc
  differentiation.cc
  init_tensors.cc
  prepare_tensors.cc
  )
psi4_add_module(bin adc sources)
