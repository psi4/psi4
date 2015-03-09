# <<<  PCMSolver Roberto Di Remigio  >>>
# PCMSolver uses the git submodules + CMake ExternalProject_Add combo
if(ENABLE_PCMSOLVER)
    include(ConfigExternal)
    set(PCMSOLVER_TESTS ON)
    get_filename_component(ZLIB_ROOT ${ZLIB_LIBRARIES} PATH)
    # PCMSolver does not know profile
    if(CMAKE_BUILD_TYPE MATCHES "profile")
        set(PCM_BUILD_TYPE "release")
    else()
        set(PCM_BUILD_TYPE ${CMAKE_BUILD_TYPE})
    endif()
    list(APPEND PCMSolverCMakeArgs
        -DCMAKE_BUILD_TYPE=${PCM_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${EXTERNAL_PROJECT_INSTALL_PREFIX}
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DEXTRA_Fortran_FLAGS=${PCM_EXTRA_Fortran_FLAGS}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DEXTRA_C_FLAGS=${PCM_EXTRA_C_FLAGS}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DEXTRA_CXX_FLAGS=${PCM_EXTRA_CXX_FLAGS}
	    -DENABLE_CXX11_SUPPORT=${ENABLE_CXX11_SUPPORT}
        -DPARENT_INCLUDE_DIR=${PROJECT_BINARY_DIR}/external/include
        -DPARENT_DEFINITIONS=${PARENT_DEFINITIONS}
        -DEIGEN3_ROOT=${EIGEN3_ROOT}
        -DDISABLE_EIGEN_OWN=${DISABLE_EIGEN_OWN}
        -DENABLE_EIGEN_MKL=${ENABLE_EIGEN_MKL}
        -DPARENT_BINARY_DIR=${PROJECT_BINARY_DIR}
	    -DBOOST_INCLUDEDIR=${BOOST_INCLUDE_DIRS}
	    -DBOOST_LIBRARYDIR=${BOOST_LIBRARIES}
        -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
        -DENABLE_TESTS=${PCMSOLVER_TESTS}
        -DZLIB_ROOT=${ZLIB_ROOT}
        -DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE}
     )

    if(PCMSOLVER_TESTS)
       add_external(pcmsolver "${PCMSolverCMakeArgs}" ${CMAKE_CTEST_COMMAND})
    else()
       add_external(pcmsolver "${PCMSolverCMakeArgs}")
    endif()

    include_directories(${EXTERNAL_PROJECT_INSTALL_PREFIX}/include
        ${EXTERNAL_PROJECT_INSTALL_PREFIX}/include/eigen3/include/eigen3)
    set(PCMSOLVER_PARSE_DIR ${EXTERNAL_PROJECT_INSTALL_PREFIX}/bin)
    configure_file(${CMAKE_SOURCE_DIR}/lib/python/pcm_placeholder.py.in
	           ${CMAKE_SOURCE_DIR}/lib/python/pcm_placeholder.py)
endif()

set(headers_list "")
# List of headers
list(APPEND headers_list psipcm.h )

# If you want to remove some headers specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM headers_list "")
else()
   list(REMOVE_ITEM headers_list "")
endif()
# Sort alphabetically
list(SORT headers_list)

set(sources_list "")
# List of sources
list(APPEND sources_list psipcm.cc )

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "")
else()
   list(REMOVE_ITEM sources_list "")
endif()

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(psipcm STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(psipcm pcmsolver)
# PCMSolver needs Fortran libraries
if(CMAKE_CXX_COMPILER_ID MATCHES GNU OR CMAKE_CXX_COMPILER_ID MATCHES Clang)
  list(APPEND PCMSOLVER_LIBS pcm getkw gfortran quadmath ${ZLIB_LIBRARIES})
else() # Intel case
  list(APPEND PCMSOLVER_LIBS pcm getkw ifcore imf ${ZLIB_LIBRARIES})
endif()
set_property(GLOBAL APPEND PROPERTY LIBLIST psipcm ${PCMSOLVER_LIBS})
if(BUILD_CUSTOM_BOOST)
   add_dependencies(pcmsolver custom_boost)
   add_dependencies(psipcm custom_boost)
endif()

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libpsipcm)
