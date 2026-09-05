include(TestingMacros)

add_regression_test(fsapt-d4 "psi;quicktests;smoketests;sapt;cart;fsapt")
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../../psi4/driver/procrouting/sapt/fsapt.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
