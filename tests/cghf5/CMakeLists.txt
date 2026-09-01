include(TestingMacros)

if(ENABLE_Einsums)
    add_regression_test(cghf5 "psi;quicktests;cghf")
    file(COPY density.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
