include(TestingMacros)

add_regression_test(nbody-convergence "psi;nbody;gradient")
set_tests_properties(
  nbody-convergence
  PROPERTIES
    COST 1500
    TIMEOUT 3000
  )
