include(TestingMacros)

add_regression_test(linK-3 "psi;scf;direct-scf")
set_tests_properties(linK-3 PROPERTIES COST 100)
