include(TestingMacros)
file(COPY psiaux1 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY psiaux2 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_regression_test(dftd3-psithon2 "psi;longertests;dftd3;ext-exe-dep")
