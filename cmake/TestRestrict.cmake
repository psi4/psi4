# - Define macro to check restrict keyword
#
#  VARIABLE will be set to the keyword
#
macro(test_restrict variable)
    if(NOT DEFINED test_${variable})
        message(STATUS "Checking for restrict keyword")

        set(_testfl ${CMAKE_BINARY_DIR}/test_restrict.c)
        file(WRITE  ${_testfl} "
        int foo(int * TESTRESTRICTDEF i, int * TESTRESTRICTDEF j){return *i+*j;}
        int main(int argc, char *argv[]){ int i=0; int j=0; return foo(&i,&j); }
        ")

        # Start with __restrict__, since that is the C++ default keyword.
        foreach(keyword "__restrict__" "__restrict" "restrict")
            if(NOT test_${variable})
                try_compile(test_${variable} "${CMAKE_BINARY_DIR}"
                            "${CMAKE_BINARY_DIR}/test_restrict.c"
                            COMPILE_DEFINITIONS "-DTESTRESTRICTDEF=${keyword}" )
                set(last_restrict_keyword ${keyword})
            endif()
        endforeach()

        if(test_${variable})
            set(${variable} ${last_restrict_keyword} CACHE STRING "Restrict keyword")
            message(STATUS "   keyword found : ${last_restrict_keyword}")
        else()
            set(${variable} " " CACHE STRING "Restrict keyword")
            message(STATUS "   keyword NOT found")
        endif()
        file(REMOVE ${_testfl})
    endif()
endmacro()
