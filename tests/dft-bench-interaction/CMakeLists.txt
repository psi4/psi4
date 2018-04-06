include(TestingMacros)

add_regression_test(dft-bench-interaction "psi;dft;scf")
set_tests_properties(dft-bench-interaction PROPERTIES TIMEOUT 3000)
set_tests_properties(dft-bench-interaction PROPERTIES COST 2500)
