include(TestingMacros)

add_regression_test(nbody-convergence "psi;nbody;gradient")
set_tests_properties(nbody-convergence PROPERTIES TIMEOUT 3000)
set_tests_properties(nbody-convergence PROPERTIES COST 1500)
