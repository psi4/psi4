include(TestingMacros)

add_regression_test(dct7 "psi;dct;noc1")
set_tests_properties(dct7 PROPERTIES COST 100)
