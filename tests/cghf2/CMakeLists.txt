include(TestingMacros)

if(ENABLE_Einsums)
    add_regression_test(cghf2 "psi;cghf")
    set_tests_properties(cghf2 PROPERTIES WILL_FAIL true)
endif()
