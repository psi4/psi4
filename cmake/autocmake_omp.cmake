# Downloaded from
#   https://github.com/coderefinery/autocmake/blob/master/modules/omp.cmake
# * moved option up
# * toggled option default to ON
# * reorganized logic for Fortran + C/CXX, see https://github.com/coderefinery/autocmake/issues/177

#.rst:
#
# Enables OpenMP support.
#
# Variables used::
#
#   ENABLE_OPENMP
#   OPENMP_FOUND
#
# Variables modified (provided the corresponding language is enabled)::
#
#   CMAKE_Fortran_FLAGS
#   CMAKE_C_FLAGS
#   CMAKE_CXX_FLAGS
#
# autocmake.yml configuration::
#
#   docopt: "--omp Enable OpenMP parallelization [default: False]."
#   define: "'-DENABLE_OPENMP={0}'.format(arguments['--omp'])"

if(ENABLE_OPENMP)

    if(NOT OPENMP_FOUND)
        find_package(OpenMP)
    endif()

    if(OPENMP_FOUND)
        add_definitions(-DHAVE_OPENMP)
        if(DEFINED CMAKE_C_COMPILER_ID)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        endif()
        if(DEFINED CMAKE_CXX_COMPILER_ID)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        endif()
        if(DEFINED CMAKE_Fortran_COMPILER_ID)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
        endif()
    endif()

    if(DEFINED CMAKE_Fortran_COMPILER_ID AND NOT DEFINED OpenMP_Fortran_FLAGS)
        # we do this in a pedestrian way because the Fortran support is relatively recent
        if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
        endif()
        if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
            if(WIN32)
                set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Qopenmp")
            elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND
                   "${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS "15.0.0.20140528")
                set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp")
            else()
                set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp")
            endif()
        endif()
        if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mp")
        endif()
        if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qsmp")
        endif()
        if(CMAKE_Fortran_COMPILER_ID MATCHES Cray)
            # do nothing in this case
        endif()
        set(OPENMP_FOUND TRUE)
    endif()
endif()
