include(TestingMacros)

add_regression_test(fd-freq-gradient-large "psi;findif;cart;noc1;medlong")
set_tests_properties(fd-freq-gradient-large PROPERTIES COST 200)
