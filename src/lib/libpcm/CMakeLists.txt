# <<<  PCMSolver Roberto Di Remigio  >>>
# PCMSolver uses the git submodules + CMake ExternalProject_Add combo
if(ENABLE_PCMSOLVER)
    set(PCMSOLVER_TESTS ON)
    get_filename_component(ZLIB_ROOT ${ZLIB_LIBRARIES} PATH)
    set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
	-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DPARENT_INCLUDE_DIR=${PROJECT_BINARY_DIR}/external/include
        -DPARENT_DEFINITIONS=${PARENT_DEFINITIONS}
        -DEIGEN3_ROOT=${EIGEN3_ROOT}
	-DDISABLE_EIGEN_OWN=${DISABLE_EIGEN_OWN}
	-DENABLE_EIGEN_MKL=${ENABLE_EIGEN_MKL}
        -DPARENT_BINARY_DIR=${PROJECT_BINARY_DIR}
	-DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS}
	-DBOOST_LIBRARYDIR=${Boost_LIBRARIES}
	-DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE}
	-DENABLE_TESTS=${PCMSOLVER_TESTS}
	-DZLIB_ROOT=${ZLIB_ROOT}
	)
    
    if(PCMSOLVER_TESTS)
	    add_external(pcmsolver ${CMAKE_CTEST_COMMAND})
    else()
	    add_external(pcmsolver)
    endif()


    include_directories(${PROJECT_BINARY_DIR}/external/include
	                ${PROJECT_BINARY_DIR}/external/include/eigen3/include/eigen3)
endif()

set(headers_list "")
# List of headers
list(APPEND headers_list psipcm.h)

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
list(APPEND sources_list psipcm.cc)

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
add_library(psipcm STATIC ${sources_list} ${headers_list})
# Specify dependencies for the library (if any)
add_dependencies(psipcm pcmsolver)
set(libs_to_merge psipcm ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(pcmsolver custom_boost)
   add_dependencies(psipcm custom_boost)
endif()
install(TARGETS psipcm ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libpcm)
