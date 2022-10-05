include(TestingMacros)

add_regression_test(fsapt-ext-abc2 "psi;quicktests;sapt;cart;fsapt;extern")
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../../psi4/share/psi4/fsapt/fsapt.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
