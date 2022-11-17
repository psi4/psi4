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

if(WIN32)
  set(_has_dggsvd3 TRUE)
  set(_has_dggsvp3 TRUE)
else()
  if(Fortran_ENABLED)
    include(CheckFortranFunctionExists)
    set(CMAKE_REQUIRED_LIBRARIES tgt::lapk)
    # Check whether DGGSVD3 is available
    check_fortran_function_exists(DGGSVD3 _has_dggsvd3)
    # Check whether DGGSVP3 is available
    check_fortran_function_exists(DGGSVP3 _has_dggsvp3)
    unset(CMAKE_REQUIRED_LIBRARIES)
  endif()
endif()

psi4_add_module(lib qt sources)
target_compile_definitions(qt
  PRIVATE
    $<$<BOOL:${_has_dggsvd3}>:LAPACK_HAS_DGGSVD3>
    $<$<BOOL:${_has_dggsvp3}>:LAPACK_HAS_DGGSVP3>
    FC_SYMBOL=${FC_SYMBOL}
  )
unset(_has_dggsvd3)
unset(_has_dggsvp3)
