set(sources_list
  3d_array.cc
  blas_level1.cc
  blas_level2.cc
  blas_level3.cc
  cc_excited.cc
  david.cc
  dirprd_block.cc
  dot_block.cc
  dx_read.cc
  dx_write.cc
  fill_sym_matrix.cc
  invert.cc
  lapack.cc
  mat_print.cc
  newmm_rking.cc
  normalize.cc
  pople.cc
  probabil.cc
  ras_set.cc
  reorder_qt.cc
  schmidt.cc
  schmidt_add.cc
  solve_pep.cc
  timer.cc
)

include(CheckSymbolExists)
if(${isMKL} MATCHES "MKL")
  set(_include "mkl_lapacke.h")
else()
  set(_include "lapacke.h")
endif()
set(CMAKE_REQUIRED_LIBRARIES tgt::lapk tgt::blas)
# Check whether LAPACKE_dggsvd3 is available
check_symbol_exists(LAPACKE_dggsvd3 ${_include} _has_dggsvd3)
# Check whether LAPACKE_dggsvp3 is available
check_symbol_exists(LAPACKE_dggsvp3 ${_include} _has_dggsvp3)
unset(CMAKE_REQUIRED_LIBRARIES)
unset(_include)

psi4_add_module(lib qt sources_list psio ciomr mints)
target_compile_definitions(qt
  PRIVATE
    $<$<BOOL:${_has_dggsvd3}>:LAPACK_HAS_DGGSVD3>
    $<$<BOOL:${_has_dggsvp3}>:LAPACK_HAS_DGGSVP3>
  )
unset(_has_dggsvd3)
unset(_has_dggsvp3)
