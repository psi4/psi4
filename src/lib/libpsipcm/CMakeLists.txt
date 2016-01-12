# <<<  PCMSolver Roberto Di Remigio  >>>
# PCMSolver uses the git submodules + CMake ExternalProject_Add combo
if(ENABLE_PCMSOLVER)
    include(ConfigExternal)
    get_filename_component(ZLIB_ROOT ${ZLIB_LIBRARIES} PATH)
    # PCMSolver does not know profile
    if(CMAKE_BUILD_TYPE MATCHES "profile")
        set(PCM_BUILD_TYPE "release")
    else()
        set(PCM_BUILD_TYPE ${CMAKE_BUILD_TYPE})
    endif()
    list(APPEND PCMSolverCMakeArgs
      -DCMAKE_BUILD_TYPE=${PCM_BUILD_TYPE}
      -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/interfaces
      -DSUBMODULES_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/interfaces
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DEXTRA_Fortran_FLAGS=${PCM_EXTRA_Fortran_FLAGS}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DEXTRA_C_FLAGS=${PCM_EXTRA_C_FLAGS}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DEXTRA_CXX_FLAGS=${PCM_EXTRA_CXX_FLAGS}
      -DENABLE_CXX11_SUPPORT=${ENABLE_CXX11_SUPPORT}
      -DBOOST_INCLUDEDIR=${BOOST_INCLUDE_DIRS}
      -DBOOST_LIBRARYDIR=${BOOST_LIBRARIES}
      -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
      -DENABLE_TESTS=OFF
      -DENABLE_LOGGER=OFF
      -DENABLE_TIMER=OFF
      -DBUILD_STANDALONE=OFF
      -DENABLE_FORTRAN_API=OFF
      -DENABLE_GENERIC=ON
      -DZLIB_ROOT=${ZLIB_ROOT}
      -DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE}
     )

    add_external(pcmsolver "${PCMSolverCMakeArgs}")

    include_directories(${PROJECT_SOURCE_DIR}/interfaces/pcmsolver/api)
    set(PCMSOLVER_PARSE_DIR ${PROJECT_BINARY_DIR}/interfaces/bin)
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


# Build static library
add_library(psipcm STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(psipcm pcmsolver)
list(APPEND PCMSOLVER_LIBS pcm "${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
set_property(GLOBAL APPEND PROPERTY LIBLIST psipcm ${PCMSOLVER_LIBS})
if(BUILD_CUSTOM_BOOST)
   add_dependencies(pcmsolver custom_boost)
   add_dependencies(psipcm custom_boost)
endif()

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libpsipcm)
