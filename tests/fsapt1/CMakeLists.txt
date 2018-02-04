include(TestingMacros)

add_regression_test(fsapt1 "psi;sapt;cart")
set_tests_properties(fsapt1 PROPERTIES COST 400)
