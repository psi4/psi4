include(TestingMacros)
file(COPY ciroot-1.dat.ref DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_regression_test(cidump-2acene "psi;cidump")
