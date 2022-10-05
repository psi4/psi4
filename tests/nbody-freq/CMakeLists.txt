include(TestingMacros)

add_regression_test(nbody-freq "psi;nbody;freq;hessian;gradient;opt")
set_tests_properties(
  nbody-freq
  PROPERTIES
    COST 1000
    TIMEOUT 3000
  )
