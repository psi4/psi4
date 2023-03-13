list(APPEND sources
  3d_array.cc
  blas_intfc.cc
  blas_intfc23.cc
  cc_excited.cc
  dirprd_block.cc
  dot_block.cc
  dx_read.cc
  dx_write.cc
  fill_sym_matrix.cc
  invert.cc
  lapack_intfc.cc
  mat_print.cc
  newmm_rking.cc
  pople.cc
  probabil.cc
  ras_set.cc
  reorder_qt.cc
  schmidt.cc
  solve_pep.cc
  timer.cc
  )

# cmake_symbol_exists() or check_cxx_source_runs() are preferable to check_function_exists()
#   according to CMake docs. But both require an explicit header, either in the call or in the code,
#   and that's tough with multiple BLAS/LAPACK backends. So, we'll see if this is robust enough.
include(CMakePushCheckState)
include(CheckFunctionExists)
cmake_push_check_state()
list(APPEND CMAKE_REQUIRED_LIBRARIES tgt::blas)
check_function_exists(daxpby _has_daxpby)
cmake_pop_check_state()

psi4_add_module(lib qt sources)
target_compile_definitions(qt
  PRIVATE
    $<$<BOOL:${_has_daxpby}>:BLAS_HAS_DAXPBY>
    FC_SYMBOL=${FC_SYMBOL}
  )
unset(_has_daxpby)
