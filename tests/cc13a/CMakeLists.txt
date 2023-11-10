include(TestingMacros)

add_regression_test(cc13a "psi;cc;medlong")
set_tests_properties(cc13a PROPERTIES COST 100)
