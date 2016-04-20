include(TestingMacros)
file(COPY grid.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_regression_test(mp2_properties "psi;quicktests;mp2;dfmp2;properties")
