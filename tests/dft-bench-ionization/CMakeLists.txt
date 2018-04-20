include(TestingMacros)

add_regression_test(dft-bench-ionization "psi;dft;scf")
set_tests_properties(dft-bench-ionization PROPERTIES TIMEOUT 3000)
set_tests_properties(dft-bench-ionization PROPERTIES COST 1300)
