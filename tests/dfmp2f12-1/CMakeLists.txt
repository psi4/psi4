include(TestingMacros)

if(ENABLE_Einsums)
    add_regression_test(dfmp2f12-1 "psi;f12;df;quicktests")
endif()
