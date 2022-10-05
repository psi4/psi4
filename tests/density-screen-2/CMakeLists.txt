include(TestingMacros)

add_regression_test(density-screen-2 "psi;dft;scf")
set_tests_properties(density-screen-2 PROPERTIES COST 500)
