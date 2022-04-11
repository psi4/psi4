include(TestingMacros)

add_regression_test(composite-2 "psi;dft;scf;composite-scf")
set_tests_properties(composite-2 PROPERTIES COST 800)
