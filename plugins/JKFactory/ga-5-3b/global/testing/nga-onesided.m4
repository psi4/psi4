divert(-1)
  # what kind of data type to test?
  define(m4_test_int, `yes')
  define(m4_test_dbl, `yes')
  define(m4_test_dcpl, `yes')

  # test dimension from which to which?
  define(m4_dim_from, 1)
  define(m4_dim_to, 7)

  # functions to test
  define(m4_test_NGA_PUT, `yes')
  define(m4_test_NGA_GET, `yes')
  define(m4_test_NGA_ACC, `yes')
  define(m4_test_NGA_SCATTER, `yes')
  define(m4_test_NGA_SCATTER_ACC, `yes')
  define(m4_test_NGA_GATHER, `yes')
divert
include(`ngatest_src/ngatest.def')
