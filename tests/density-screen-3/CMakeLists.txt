include(TestingMacros)

add_regression_test(density-screen-3 "psi;dft;scf;medlong")
set_tests_properties(density-screen-3 PROPERTIES COST 500)
