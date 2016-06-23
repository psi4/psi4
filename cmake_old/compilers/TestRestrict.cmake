# - Define macro to check restrict keyword
#
#  VARIABLE will be set to the keyword
#
macro(test_restrict variable)
    if(NOT DEFINED test_${variable})
       message(STATUS "Checking for restrict keyword")

       # Start with __restrict__, since that is the C++ default keyword. 
       foreach(keyword "__restrict__" "__restrict" "restrict")
           if(NOT test_${variable})
               try_compile(test_${variable} "${CMAKE_BINARY_DIR}"    
                           "${CMAKE_SOURCE_DIR}/cmake/compilers/test_restrict.c"
                           COMPILE_DEFINITIONS "-DTESTRESTRICTDEF=${keyword}" )
               set(last_restrict_keyword ${keyword})
           endif(NOT test_${variable})
       endforeach(keyword)

       if(test_${variable})
           set(${variable} ${last_restrict_keyword} CACHE STRING "Restrict keyword")
           message(STATUS "   keyword found : ${last_restrict_keyword}")
       else(test_${variable})
           set(${variable} " " CACHE STRING "Restrict keyword")
           message(STATUS "   keyword NOT found")
       endif(test_${variable})
    endif(NOT DEFINED test_${variable})        
endmacro(test_restrict variable)
