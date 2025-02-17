include(TestingMacros)

if(ENABLE_Einsums)
    add_regression_test(mp2f12-1 "psi;f12")
endif()