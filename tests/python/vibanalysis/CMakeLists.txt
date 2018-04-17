include(TestingMacros)

add_regression_test(python-vibanalysis "psi;freq;python")
set_tests_properties(python-vibanalysis PROPERTIES COST 700)
