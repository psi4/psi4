include(TestingMacros)

add_regression_test(scf-bz2 "psi;scf")
set_tests_properties(scf-bz2 PROPERTIES COST 100)
