# Based on Ryan M. Richard work in PSI4
# First of all get BOOSTVER underscore separated
string(REGEX REPLACE "\\." "_" BOOSTVER ${BOOSTVER})
# Transform the ;-separated list to a ,-separated list (digested by the Boost build toolchain!)
string(REPLACE ";" "," b2_needed_components "${needed_components}")

set(BOOST_BUILD_DIR ${PROJECT_BINARY_DIR}/boost/boost_${BOOSTVER})

# Prepare stuff to be written to user-config.jam and project-config.jam
# We are going to build only what is necessary. We need to create a custom user-config.jam file
# To get boost to compile MPI we need to append "using mpi ;" to the end of the
# user-config.jam file. MPI_SENT will be the command we append
set(MPI_SENT "")
if(ENABLE_MPI)
   list(APPEND needed_components mpi)
   set(MPI_SENT "using mpi \n : ${CMAKE_CXX_COMPILER} \n;")
endif()
if(ENABLE_UNIT_TESTS)
   # Replace unit_test_framework (used by CMake's find_package) with test (understood by Boost build toolchain)
   string(REPLACE "unit_test_framework" "test" b2_needed_components "${b2_needed_components}")
endif()
# For a prettier printout...
string(REPLACE ";" ", " printout "${needed_components}")
message(STATUS "  Libraries to be built: ${printout}")

# To get boost to be built with our compiler of choice we need these logic
# statements, then we write the options to user-config.jam
set(toolset   "")
if(HAS_CXX11_SUPPORT)
   set(cxxflags "${cxxflags} ${CXX11_COMPILER_FLAGS}") 
endif()	  
set(linkflags "")
if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
	set(toolset "intel-linux")
elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  set(toolset "clang")
  if(CMAKE_SYSTEM_NAME MATCHES Darwin) # Read here: http://stackoverflow.com/a/19774902/2528668
  	set(cxxflags "${cxxflags} -stdlib=libc++")
	set(linkflags "${linkflags} -stdlib=libc++")
  endif()	  
else()
  if(CMAKE_SYSTEM_NAME MATCHES Darwin)
     set(toolset "darwin")
  else()
     set(toolset "gcc")
  endif()
endif()
message(STATUS "  Toolset to be used: ${toolset}")

# Write the user-config.jam file
file(WRITE  ${BOOST_BUILD_DIR}/user-config.jam 
	    "using ${toolset}\n : \n : ${CMAKE_CXX_COMPILER} \n : <cxxflags>${cxxflags} \n : <linkflags>${linkflags} \n;") 
file(APPEND ${BOOST_BUILD_DIR}/user-config.jam "${MPI_SENT}") 

# Unpack Boost and copy boostcpp.jam
add_custom_command(
    OUTPUT boost.unpacked 
    COMMAND ${CMAKE_COMMAND} -E tar xjf ${PROJECT_SOURCE_DIR}/boost/boost_${BOOSTVER}.tar.bz2 
    # Copy custom boostcpp.jam (to get versioned format as I want it!)
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/boost/boostcpp.jam ${BOOST_BUILD_DIR}
    # Append to project-config.jam 
    # It can't be done with file(APPEND ...) because it will be destroyed when unpacking
    COMMAND ${CMAKE_COMMAND} -E echo "${MPI_SENT}" >> ${BOOST_BUILD_DIR}/project-config.jam
    COMMAND touch boost.unpacked
    DEPENDS ${PROJECT_SOURCE_DIR}/boost/boost_${BOOSTVER}.tar.bz2
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/boost
    COMMENT "Unpacking Boost")

string(TOLOWER ${CMAKE_BUILD_TYPE} type)
# Run bootstrap.sh to configure the build. We will install in ${PROJECT_BINARY_DIR}/boost
add_custom_command(
    OUTPUT boost.configured
    COMMAND ./bootstrap.sh --with-toolset=${toolset}
            --with-libraries="${b2_needed_components}"
	    --with-python=${PYTHON_EXECUTABLE}
	    --prefix=${CUSTOM_BOOST_LOCATION} 1> boost.configured.log 2> boost.configured.err
    COMMAND ${CMAKE_COMMAND} -E touch ${PROJECT_BINARY_DIR}/boost/boost.configured	    
    WORKING_DIRECTORY ${BOOST_BUILD_DIR}
    DEPENDS boost.unpacked
    COMMENT "Configuring Boost")

# Build Boost
add_custom_command(
    OUTPUT boost.built
    COMMAND ./b2 variant=${type} link=static threading=multi 1> boost.built.log 2> boost.built.err
    COMMAND ${CMAKE_COMMAND} -E touch ${PROJECT_BINARY_DIR}/boost/boost.built
    WORKING_DIRECTORY ${BOOST_BUILD_DIR}
    DEPENDS boost.configured 
    COMMENT "Building Boost")

# Install Boost
add_custom_command(
    OUTPUT boost.installed
    COMMAND ./b2 install variant=${type} link=static threading=multi 1> boost.installed.log 2> boost.installed.err	   
    COMMAND ${CMAKE_COMMAND} -E touch ${PROJECT_BINARY_DIR}/boost/boost.installed
    WORKING_DIRECTORY ${BOOST_BUILD_DIR}
    DEPENDS boost.built 
    COMMENT "Installing Boost")

# Clean-up
add_custom_command(
    OUTPUT boost.cleanedup
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${BOOST_BUILD_DIR}
    COMMAND ${CMAKE_COMMAND} -E touch ${PROJECT_BINARY_DIR}/boost/boost.cleanedup
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/boost
    DEPENDS boost.installed 
    COMMENT "Clean-up Boost")

add_custom_target(custom_boost DEPENDS boost.cleanedup)
