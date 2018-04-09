include(TestingMacros)

add_regression_test(sapt4 "psi;longtests;sapt")
set_tests_properties(sapt4 PROPERTIES COST 100)
