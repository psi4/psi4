include(ExternalProject)

### >> Boost version to use internally <<
set(Boost_INTERNAL_VERSION 1.57.0)

# Replace . with _ in version.
string(REPLACE "." "_" Boost_INTERNAL_VERSION_UNDERSCORED "${Boost_INTERNAL_VERSION}")

### >> Boost components currently required <<
set(BOOST_COMPONENTS filesystem python regex serialization system timer chrono
                     thread
)

find_package(Boost ${Boost_INTERNAL_VERSION} QUIET COMPONENTS ${BOOST_COMPONENTS})

if(NOT Boost_FOUND)
   message(STATUS "A suitable version of Boost could not be located.  One will"
                  " be built for you.")
   string(REPLACE ";" "," B2_BOOST_COMPONENTS "${BOOST_COMPONENTS}")
   string(TOLOWER "${CMAKE_BUILD_TYPE}" BOOST_BUILD_TYPE)

   if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
      set(BOOST_COMPILER_TYPE "intel-linux")
   elseif(CMAKE_CXX_COMPILER_ID MATCHES GNU)
      set(BOOST_COMPILER_TYPE "gcc")
   elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
      set(BOOST_COMPILER_TYPE "clang")
   elseif(CMAKE_CXX_COMPILER_ID MATCHES Darwin)
      set(BOOST_COMPILER_TYPE "darwin")
   endif()
   set(BOOST_PREFIX ${CMAKE_BINARY_DIR}/external/boost/boost-src/boost_${Boost_INTERNAL_VERSION_UNDERSCORED})
   set(BOOST_INSTALL_PREFIX 
       ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/external/boost)
   #Can't change the fact that boost has warnings....   
   set(BOOST_FLAGS "${CMAKE_CXX_FLAGS} -w")
   if(BUILD_FPIC)
      set(BOOST_FLAGS "${BOOST_FLAGS} -fPIC")
   endif()
   set(BOOST_LINK "static")
   if(BUILD_SHARED_LIBS)
      set(BOOST_LINK "shared")
   endif()

   set(Boost_USER_CONFIG_JAMFILE "${CMAKE_CURRENT_BINARY_DIR}/user-config.jam")
   file(WRITE "${Boost_USER_CONFIG_JAMFILE}" "using ${BOOST_COMPILER_TYPE} : : ${CMAKE_CXX_COMPILER} ;")
   ExternalProject_Add(boost_external
      PREFIX ${BOOST_PREFIX}
      URL ${CMAKE_CURRENT_SOURCE_DIR}/boost-src/boost_${Boost_INTERNAL_VERSION_UNDERSCORED}.tar.bz2
      URL_HASH MD5=1BE49BEFBDD9A5CE9DEF2983BA3E7B76
      BUILD_IN_SOURCE 1
      CONFIGURE_COMMAND
         ./bootstrap.sh --with-toolset=${BOOST_COMPILER_TYPE}
                                  --with-libraries=${B2_BOOST_COMPONENTS}
                                  --with-python=${PYTHON_EXECUTABLE}
                                  #--with-python-root=${PYTHON_INCLUDE_DIR}
                                  --prefix=${BOOST_INSTALL_PREFIX}

      BUILD_COMMAND ./b2 threading=multi 
                         link=${BOOST_LINK}   
                         cxxflags=${BOOST_FLAGS}
                         --user-config=${Boost_USER_CONFIG_JAMFILE}
      INSTALL_DIR ${BOOST_INSTALL_PREFIX}
      INSTALL_COMMAND ./b2 install threading=multi 
                                   link=${BOOST_LINK} 
                                   cxxflags=${BOOST_FLAGS}
                                   --user-config=${Boost_USER_CONFIG_JAMFILE}
   )
else()
    add_library(boost_external INTERFACE)
endif()
