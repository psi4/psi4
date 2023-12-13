include(TestingMacros)

add_regression_test(linK-2 "psi;dft;scf;direct-scf;medlong")
set_tests_properties(linK-2 PROPERTIES COST 800)
