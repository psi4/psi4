include(TestingMacros)
file(COPY psiaux1 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY psiaux2 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_regression_test(psithon2 "psi;misc")
