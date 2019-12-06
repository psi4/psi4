#find_package(Libint2 2.6 CONFIG QUIET COMPONENTS ${MAX_AM_ERI})
find_package(Libint2 CONFIG COMPONENTS gss)

if(${Libint2_FOUND})
    # CHOOSE! `make export`
    #get_property(_loc TARGET Libint2::int2 PROPERTY LOCATION)
    # CHOOSE! pure cmake
    get_property(_loc TARGET Libint2::cxx PROPERTY LOCATION)
    message(STATUS "${Cyan}Found Libint2 ${Libint_MAX_AM_ERI}${ColourReset}: ${_loc} (found version ${Libint2_VERSION})")
    add_library(libint2_external INTERFACE)  # dummy
else()
    if(${CMAKE_INSIST_FIND_PACKAGE_Libint2})
        message(FATAL_ERROR "Suitable Libint2 could not be externally located as user insists")
    endif()

    include(ExternalProject)

    if(${BUILD_SHARED_LIBS})
        set(_a_only  OFF)
        set(_so_only ON)
    else()
        set(_a_only  ON)
        set(_so_only OFF)
    endif()

    if(NOT ${BUILD_Libint2_GENERATOR})
        message(STATUS "Suitable Libint2 could not be located, ${Magenta}Building Libint2${ColourReset} instead.")

        ExternalProject_Add(libint2_external
#            #URL https://github.com/evaleev/libint/archive/b13e71d3cf9960460c4019e5ecf2546a5f361c71.zip
#            #URL https://github.com/evaleev/libint/releases/download/v2.6.0/libint-2.6.0.tgz
#            URL https://github.com/loriab/libint/archive/l2cmake.zip
            URL "/home/psilocaluser/tars/Libint2-export.tgz"
            CMAKE_ARGS -GNinja
                       -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DBOOST_ROOT=${BOOST_ROOT}
                       -DEigen3_DIR=${Eigen3_DIR}
                       -DMPFR_ROOT=${MPFR_ROOT}
                       -DLIBINT2_SHGAUSS_ORDERING=gaussian
                       -DBUILD_SHARED=${_so_only}
                       -DBUILD_STATIC=${_a_only}
            CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS})

    else()
        message(STATUS "Suitable Libint2 could not be located, ${Magenta}Building Libint2${ColourReset} from generator source instead.")

        ExternalProject_Add(libint2_external
            # "git checkout" fails on Windows, because of "*" in filenames (e.g. basis set files)
            #URL https://github.com/evaleev/libint/archive/b13e71d3cf9960460c4019e5ecf2546a5f361c71.zip
            #URL https://github.com/evaleev/libint/releases/download/v2.6.0/libint-2.6.0.tgz
            URL https://github.com/loriab/libint/archive/l2cmake3.zip
            CMAKE_ARGS -GNinja
                       -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DENABLE_ONEBODY=2  # ?
                       -DDISABLE_ONEBODY_PROPERTY_DERIVS=ON
                       -DMULTIPOLE_MAX_ORDER=4
                       -DENABLE_ERI=2
                       -DENABLE_ERI3=2
                       -DENABLE_ERI2=2
                       -DBOOST_ROOT=${BOOST_ROOT}
                       -DEigen3_DIR=${Eigen3_DIR}
                       -DMPFR_ROOT=${MPFR_ROOT}
                       -DERI3_PURE_SH=OFF
                       -DERI2_PURE_SH=OFF
                       -DLIBINT2_SHGAUSS_ORDERING=gaussian
                       -DBUILD_SHARED=${_so_only}
                       -DBUILD_STATIC=${_a_only}
            CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DWITH_ERI_MAX_AM:STRING=5;5;5
                             -DWITH_ERI3_MAX_AM:STRING=6
                             -DWITH_ERI2_MAX_AM:STRING=6
        )

    endif()

    # CHOOSE! `make export`
    #set(Libint2_DIR ${STAGED_INSTALL_PREFIX}/lib/cmake/libint2 CACHE PATH "path to internally built libint2-config.cmake" FORCE)
    # CHOOSE! pure cmake
    set(Libint2_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/Libint2 CACHE PATH "path to internally built Libint2Config.cmake" FORCE)

endif()

