include(TestingMacros)

add_regression_test(fsapt-ext-abc "psi;quicktests;sapt;cart;fsapt;extern")
message(STATUS "${CMAKE_CURRENT_SOURCE_DIR}/../../psi4/driver/procrouting/sapt/fsapt.py ${CMAKE_CURRENT_BINARY_DIR}")
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../../psi4/driver/procrouting/sapt/fsapt.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
