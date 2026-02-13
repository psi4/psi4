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
list(APPEND CMAKE_REQUIRED_LIBRARIES tgt::lapack)
check_function_exists(daxpby _has_daxpby)
if(NOT _has_daxpby)
  message(WARNING "${Yellow}Your BLAS/LAPACK library does not seem to be providing the DAXPBY subroutine. DAXPBY calls will be re-routed through DSCAL and DAXPY. If you are seeing this message, and you are not building with Apple Accelerate, you may want to re-check that the correct BLAS/LAPACK libraries are being used by the build system. Redo CMake configuration for changes to register.${ColourReset}")
endif()

check_function_exists(dggsvd3 _has_dggsvd3)
check_function_exists(dggsvp3 _has_dggsvp3)
if(NOT (_has_dggsvd3 AND _has_dggsvp3))
  message(WARNING "${Yellow}Your BLAS/LAPACK library does not seem to be providing the DGGSVD3 and DGGSVP3 subroutines. No re-routing is available. If you are seeing this message, and you are not building with Apple Accelerate, you may want to re-check that the correct BLAS/LAPACK libraries are being used by the build system. Redo CMake configuration for changes to register.${ColourReset}")
endif()
cmake_pop_check_state()

psi4_add_module(lib qt sources)
target_compile_definitions(qt
  PRIVATE
    $<$<BOOL:${_has_daxpby}>:BLAS_HAS_DAXPBY>
    $<$<BOOL:${_has_dggsvd3}>:BLAS_HAS_DGGSVD3>
    $<$<BOOL:${_has_dggsvp3}>:BLAS_HAS_DGGSVP3>
    FC_SYMBOL=${FC_SYMBOL}
  )

unset(_has_daxpby)
unset(_has_dggsvd3)
unset(_has_dggsvp3)
