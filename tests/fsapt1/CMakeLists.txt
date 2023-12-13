include(TestingMacros)

add_regression_test(fsapt1 "psi;sapt;cart;fsapt;medlong")
set_tests_properties(fsapt1 PROPERTIES COST 400)
