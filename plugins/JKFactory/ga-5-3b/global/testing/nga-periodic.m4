divert(-1)
  # what kind of data type to test?
  define(m4_test_int, `yes')
  define(m4_test_dbl, `yes')
  define(m4_test_dcpl, `yes')

  # test dimension from which to which?
  define(m4_dim_from, 1)
  define(m4_dim_to, 7)

  # periodic functions
  define(m4_test_NGA_PERIODIC_GET, `yes') 
  define(m4_test_NGA_PERIODIC_PUT, `yes') 
  define(m4_test_NGA_PERIODIC_ACC, `yes') 

divert
include(`ngatest_src/ngatest.def')
