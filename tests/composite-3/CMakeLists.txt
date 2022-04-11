include(TestingMacros)

add_regression_test(composite-3 "psi;scf;composite-scf")
set_tests_properties(composite-3 PROPERTIES COST 100)
