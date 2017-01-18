include(TestingMacros)
file(COPY ESP.cube.ref DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY Dt.cube.ref DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_regression_test(cubeprop-esp "psi;cubeprop")
