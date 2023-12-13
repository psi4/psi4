include(TestingMacros)

add_regression_test(cbs-xtpl-func "psi;cbs;medlong")
set_tests_properties(cbs-xtpl-func PROPERTIES COST 200)
