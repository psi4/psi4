function(test_fortran_mangling CODE PREFIX ISUPPER POSTFIX DOC SUB RESULT)
  if(ISUPPER)
    string(TOUPPER "${SUB}" sub)
  else(ISUPPER) 
    string(TOLOWER "${SUB}" sub)
  endif(ISUPPER)
  set(FUNCTION "${PREFIX}${sub}${POSTFIX}")
  # create a fortran file with sub called sub
  # 
  set(TMP_DIR
    "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/CheckFortranLink")
  file(REMOVE_RECURSE "${TMP_DIR}")
  file(WRITE "${TMP_DIR}/test.f" "${CODE}"    )
  message(STATUS "checking Fortran ${DOC} linkage: ${FUNCTION}")
  file(WRITE "${TMP_DIR}/ctof.c"
    "
      extern ${FUNCTION}();
      int main() { ${FUNCTION}(); return 0;}
    "
    )
  file(WRITE "${TMP_DIR}/CMakeLists.txt"
    "
     project(testf C Fortran)
     if(APPLE)
         set(CMAKE_C_FLAGS -m32)
     endif()
     add_library(flib test.f)
     add_executable(ctof ctof.c)
     target_link_libraries(ctof flib)
    "
    )
  set(FORTRAN_NAME_MANGLE_TEST FALSE)
  try_compile(FORTRAN_NAME_MANGLE_TEST "${TMP_DIR}" "${TMP_DIR}"
    testf
    OUTPUT_VARIABLE output)
  #if(output)
  #  message(${output})
  #  endif()
  if(FORTRAN_NAME_MANGLE_TEST)
    set(${RESULT} TRUE PARENT_SCOPE)
  else()
    set(${RESULT} FALSE PARENT_SCOPE)
  endif()
endfunction(test_fortran_mangling)

function(get_fc_symbol FCSYMBOLOUT)
    set(TESTCODE 
    "
      subroutine sub
      end subroutine sub
    ")
    #test_fortran_mangling(    CODE    pre isUpper post print_test sub worked )
    test_fortran_mangling("${TESTCODE}" "" True    "_" "FUNCTION_" "sub" FC_LINK_WORKED)
    if(FC_LINK_WORKED)
        set(${FCSYMBOLOUT} 4 PARENT_SCOPE)
        message(STATUS "Upper case with underscore is used")
        return()
    endif()
    test_fortran_mangling("${TESTCODE}" "" False   "_" "function_" "sub" FC_LINK_WORKED)
    if(FC_LINK_WORKED)
        set(${FCSYMBOLOUT} 2 PARENT_SCOPE)
        message(STATUS "Lower case with underscore is used")
        return()
    endif()
    test_fortran_mangling("${TESTCODE}" "" True    ""  "FUNCTION"  "sub" FC_LINK_WORKED)
    if(FC_LINK_WORKED)
        set(${FCSYMBOLOUT} 3 PARENT_SCOPE)
        message(STATUS "Upper case (no underscore) is used")
        return()
    endif()
    test_fortran_mangling("${TESTCODE}" "" False   ""  "function"  "sub" FC_LINK_WORKED)
    if(FC_LINK_WORKED)
        set(${FCSYMBOLOUT} 1 PARENT_SCOPE)
        message(STATUS "Lower case (no underscore) is used")
        return()
    endif()
    message(FATAL_ERROR "Unable to detect Fortran name mangling")
endfunction(get_fc_symbol)
