include(TestingMacros)

add_regression_test(fisapt-siao1 "psi;sapt")

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../../psi4/share/psi4/fsapt/fsapt.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
