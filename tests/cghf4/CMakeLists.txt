include(TestingMacros)

if(ENABLE_Einsums)
    add_regression_test(cghf4 "psi;quicktests;cghf")
endif()
