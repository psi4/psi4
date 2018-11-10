set(sources_list lapack_intfc.cc
                 dx_write.cc
                 dirprd_block.cc
                 pople.cc
                 schmidt_add.cc
                 mat_print.cc
                 ras_set.cc
                 fill_sym_matrix.cc
                 cc_excited.cc
                 probabil.cc
                 dot_block.cc
                 timer.cc
                 dx_read.cc
                 blas_intfc.cc
                 normalize.cc
                 newmm_rking.cc
                 3d_array.cc
                 blas_intfc23.cc
                 reorder_qt.cc
                 schmidt.cc
                 invert.cc
                 solve_pep.cc
                 david.cc
)
add_definitions("-DFC_SYMBOL=${FC_SYMBOL}")

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

psi4_add_module(lib qt sources_list psio ciomr mints)
target_compile_definitions(qt
  PRIVATE
    $<$<BOOL:${_has_dggsvd3}>:LAPACK_HAS_DGGSVD3>
    $<$<BOOL:${_has_dggsvp3}>:LAPACK_HAS_DGGSVP3>
  )
unset(_has_dggsvd3)
unset(_has_dggsvp3)
