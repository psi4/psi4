include(TestingMacros)

add_regression_test(cghf2 "psi;quicktests;cghf")
set_tests_properties(cghf2 PROPERTIES WILL_FAIL true)
