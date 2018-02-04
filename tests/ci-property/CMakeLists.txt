include(TestingMacros)
file(COPY grid.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_regression_test(ci-property "psi;quicktests;ci;cas;properties;cart")
