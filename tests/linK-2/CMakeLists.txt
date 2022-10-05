include(TestingMacros)

add_regression_test(linK-2 "psi;dft;scf;direct-scf")
set_tests_properties(linK-2 PROPERTIES COST 800)
