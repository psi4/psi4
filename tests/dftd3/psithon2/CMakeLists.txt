include(TestingMacros)
file(COPY psiaux1 DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../dftd3-psithon2)
file(COPY psiaux2 DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../dftd3-psithon2)
add_regression_test(dftd3-psithon2 "psi;dftd3;addon")
