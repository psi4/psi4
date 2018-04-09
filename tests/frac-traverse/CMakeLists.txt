include(TestingMacros)

add_regression_test(frac-traverse "psi;misc")
set_tests_properties(frac-traverse PROPERTIES COST 100)
