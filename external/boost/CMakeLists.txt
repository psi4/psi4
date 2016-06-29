include(ExternalProject)
set(BOOST_COMPONENTS filesystem python regex serialization system timer chrono
                     thread
)

find_package(Boost 1.57.0 QUIET COMPONENTS ${BOOST_COMPONENTS})

if(NOT Boost_FOUND)
   message(STATUS "A suitable version of Boost could not be located.  One will"
                  " be built for you.")
   string(REPLACE ";" "," B2_BOOST_COMPONENTS "${BOOST_COMPONENTS}")
   
   string(TOLOWER ${CMAKE_BUILD_TYPE} BOOST_BUILD_TYPE)

   if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
      set(BOOST_COMPILER_TYPE "intel-linux")
   elseif(CMAKE_CXX_COMPILER_ID MATCHES GNU)
      set(BOOST_COMPILER_TYPE "gcc")
   elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
      set(BOOST_COMPILER_TYPE "clang")
   elseif(CMAKE_SYSTEM_NAME MATCHES Darwin)
      set(BOOST_COMPILER_TYPE "darwin")
   endif()

   set(BOOST_PREFIX 
       ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/external/boost)
   ExternalProject_Add(boost_external
      PREFIX ${CMAKE_BINARY_DIR}/external/boost/boost-src/boost_1_57_0
      URL ${CMAKE_CURRENT_SOURCE_DIR}/boost-src/boost_1_57_0.tar.bz2
      BUILD_IN_SOURCE 1
      CONFIGURE_COMMAND pwd && ./bootstrap.sh --with-toolset=${BOOST_COMPILER_ID}
                                       --with-libraries=${B2_BOOST_COMPONENTS}
                                       --with-python=${PYTHON_EXECUTABLE}
                                       --prefix=${BOOST_PREFIX}
      BUILD_COMMAND ./b2 threading=multi link=static
      INSTALL_DIR ${BOOST_PREFIX}
      INSTALL_COMMAND ./b2 install threading=multi link=static
   )
else()
    add_library(boost_external INTERFACE)
endif()