include(TestingMacros)

if(ENABLE_Einsums)
    add_regression_test(dlpnocc-5 "psi;dlpno;cc;quick")
endif()
