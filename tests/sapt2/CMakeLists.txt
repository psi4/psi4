include(TestingMacros)

add_regression_test(sapt2 "psi;longtests;sapt;cart")
set_tests_properties(sapt2 PROPERTIES COST 100)
