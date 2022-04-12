include(TestingMacros)

add_regression_test(fd-freq-energy-large "psi;longtests;findif;cart;noc1")
set_tests_properties(
  fd-freq-energy-large
  PROPERTIES
    COST 1000
    TIMEOUT 3000
  )
