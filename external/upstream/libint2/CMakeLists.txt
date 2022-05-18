# NOTES on Libint2 build options for Psi4
# ---------------------------------------
# * In order for Psi4 to at all work with a libint2 build, need LIBINT2_CARTGAUSS_ORDERING=standard *at generation time*. This is already the default.
# * In order for Psi4 to at all work with a libint2 build, need LIBINT2_SHELL_SET=standard *at generation time*. This is already the default.
# * In order for Psi4 to at all work with a libint2 build, need LIBINT2_SHGAUSS_ORDERING=gaussian *at library build time*.
#   This is not the default, though it may become runtime selectable in future https://github.com/evaleev/libint/issues/190#issuecomment-691485225 .
# * In order for any DF to work, need ENABLE_ERI3 and ENABLE_ERI2 >=0
# * In order for gradient and Hessian tests to not segfault, need ENABLE_ERI and ENABLE_ERI3 and ENABLE_ERI2 =2
# * (new c. Jan 2022) In order for Psi4 to at all work with a libint2 build for one-electron ints, need ENABLE_ONEBODY =2
# * The `find_package(Libint2 ... COMPONENTS ...)` commands will check the above requirements met, provided the libint2
#   was built with cmake+cmake, not libtool+cmake.
# * Provided you know what you're doing wrt the above options and know what subset of Psi4 calculations you're targeting
#   (e.g., no Hessians or non-df only) and recongnize that DF is everywhere (e.g., SAD guess and DF_SCF_GUESS) and may
#   need further suppression, it's fine to reduce or eliminate the `find_package` components.
# * Psi4 can run a libtool+cmake, but one needs to disable the `find_package` components check and correct for some
#   Libint2 vs. libint2 spelling differences. See "CHOOSE!" toggles here and in <psi4>/psi4/CMakeLists.txt
# * In order for DF tests to yield right values when cartesian basis sets (6d) used, need ERI3_PURE_SH and ERI2_PURE_SH =OFF
# * In order for DF integrals to be present as nominally expected, need WITH_MAX_AM to mirror WITH_ERI3_MAX_AM
# * In order for DF calculations to not be limited by cheap integrals, need WITH_MAX_AM_ERI2 >= WITH_MAX_AM_ERI3
# * In order for DF fallback fitting to work with def2-universal on all elements, need WITH_ERI3_MAX_AM >=4 for energies
# * Since non-DF calculations often use DF guess or preiterations and DF calculations use 4-center integrals for sieving,
#   a natural relationship between ERI and ERI3 for energy integrals is N and N+1.
# * The combination WITH_MAX_AM_ERI="5;4;3" and WITH_MAX_AM_ERI3="6;5;4" permits any 5Z energy, QZ gradient, and TZ Hessian and is sufficient for the full test suite.
# * Psi4 doesn't know the limits of the Libint2 library to which it's linked. `libint2::Engine::lmax_exceeded` errors will be thrown from l2; nothing from Psi4.
# * If ever use a libint2 with AM>7, edit the `init_spherical_harmonics` line in libmints/integral.cc
# * In selecting AM values, before `grep`ing lots of basis set files, consult the guide at end of this file.

# CHOOSE! `make export`
#find_package(Libint2 CONFIG)
# CHOOSE! pure cmake
find_package(Libint2 CONFIG COMPONENTS gss CXX_ho impure_sh "eri_c4_d0_l${MAX_AM_ERI}" eri_c4_d1_l2 eri_c4_d2_l2 eri_c3_d0_l4 eri_c3_d1_l3 eri_c3_d2_l3 eri_c2_d0_l4 eri_c2_d1_l3 eri_c2_d2_l3 onebody_d0_l4 onebody_d1_l3 onebody_d2_l3)

if(${Libint2_FOUND})
    # CHOOSE! `make export`
    #get_property(_loc TARGET Libint2::int2 PROPERTY LOCATION)
    # CHOOSE! pure cmake
    get_property(_loc TARGET Libint2::int2 PROPERTY LOCATION)
    message(STATUS "${Cyan}Found Libint2 ${Libint2_MAX_AM_ERI}${ColourReset}: ${_loc} (found version ${Libint2_VERSION})")
    add_library(libint2_external INTERFACE)  # dummy
else()
    if(${CMAKE_INSIST_FIND_PACKAGE_Libint2})
        message(FATAL_ERROR "Suitable Libint2 could not be externally located as user insists")
    endif()

    include(ExternalProject)

    if(NOT ${BUILD_Libint2_GENERATOR})
        if (${MAX_AM_ERI} GREATER 5)
            set(_url_am_src "7-7-4-7-7-5")  # 22920  # passes E&G 7z conv and 6z df, H Qz conv and df. covers all internal basis sets for energy and gradient
        elseif (${MAX_AM_ERI} EQUAL 5)
            set(_url_am_src "5-4-3-6-5-4")  # 12383  # passes full. production Mac, temp production Linux_7
        elseif (${MAX_AM_ERI} LESS_EQUAL 4)
            set(_url_am_src "4-3-3-5-4-4")  # 10498  # passes quick tests with Bi->Sb mod
          # set(_url_am_src "3-2-2-4-3-3")  #  5485  # passes find_package(Libint2 COMPONENTS ...)
        endif()

        if(DEFINED ENV{LGTM_SRC})
          # set(_url_am_src "2-0-0-0-0-0")  #   389  # LGTM_SRC bulids
            set(_url_am_src "1-0-0-0-0-0")  #    99  # LGTM_SRC bulids
        endif()

        #TODO set(_url_l2_tarball "https://github.com/loriab/libint/releases/download/v0.1/Libint2-export-${_url_am_src}_1.tgz")
        set(_url_l2_tarball "https://github.com/loriab/libint/releases/download/v0.1/Libint2-export-5-4-3-6-5-4_mm4f12ob2.tgz")

        message(STATUS "Suitable Libint2 could not be located, ${Magenta}Building Libint2 ${_url_am_src}${ColourReset} instead.")

        if (MSVC)
            # Windows shared (dll) can't work
            set(_build_shared_libs "OFF")
        else()
            # stopgap to avert a report that static lib not working on Linux
            set(_build_shared_libs "ON")
        endif()

        ExternalProject_Add(libint2_external
            URL ${_url_l2_tarball}
            CMAKE_ARGS -GNinja
                       -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       #-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_BUILD_TYPE=Release
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DBOOST_ROOT=${BOOST_ROOT}
                       -DEigen3_DIR=${Eigen3_DIR}
                       -DLIBINT2_SHGAUSS_ORDERING=gaussian
                       #-DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -DBUILD_SHARED_LIBS=${_build_shared_libs}
                       -DREQUIRE_CXX_API=ON
                       -DREQUIRE_CXX_API_COMPILED=OFF
                       -DENABLE_FORTRAN=OFF
                       -DENABLE_XHOST=${ENABLE_XHOST}
                       -DBUILD_FPIC=${BUILD_FPIC}
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS})

    else()
        message(STATUS "Suitable Libint2 could not be located, ${Magenta}Building Libint2${ColourReset} from generator source instead.")
        message(WARNING "${Yellow}Libint2 generator source build requires Ninja (not Makefile) builder and Boost, Eigen, and MPFR dependencies. The build will be very long for non-trivial AM, and parallelism won't help for parts. Edit external/upstream/libint2/CMakeLists.txt for AM settings.${ColourReset}")

        if (MSVC)
            message(FATAL_ERROR "${Red}Libint2 generator source build faulty on Windows.${ColourReset}")
        endif()

        # NOTE: active AM settings below will pass few tests, as they're set for 1st row triple-zeta ene and double-zeta ene/grad/hess for both conv/df. Please edit for desired use.

        ExternalProject_Add(libint2_external
            URL https://github.com/loriab/libint/archive/new-cmake-harness-lab-rb1.zip
            CMAKE_ARGS -GNinja
                       -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DENABLE_ONEBODY=2
                       -DENABLE_ERI=2
                       -DENABLE_ERI3=2
                       -DENABLE_ERI2=2
                       -DBOOST_ROOT=${BOOST_ROOT}
                       -DEigen3_DIR=${Eigen3_DIR}
                       -DEigen3_ROOT=${Eigen3_ROOT}
                       -DMultiprecision_ROOT=${Multiprecision_ROOT}
                       -DERI3_PURE_SH=OFF
                       -DERI2_PURE_SH=OFF
                       -DLIBINT2_SHGAUSS_ORDERING=gaussian
                       -DLIBINT2_CARTGAUSS_ORDERING=standard
                       -DLIBINT2_SHELL_SET=standard
                       -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -DREQUIRE_CXX_API=ON
                       -DREQUIRE_CXX_API_COMPILED=OFF
                       -DENABLE_FORTRAN=OFF
                       -DBUILD_TESTING=${BUILD_TESTING}
                       -DENABLE_XHOST=${ENABLE_XHOST}
                       -DBUILD_FPIC=${BUILD_FPIC}
            CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DWITH_ERI_MAX_AM:STRING=3;2;2
                             -DWITH_ERI3_MAX_AM:STRING=4;3;3
                             -DWITH_ERI2_MAX_AM:STRING=4;3;3
                             -DWITH_MAX_AM:STRING=4;3;3
        )
    endif()

    # CHOOSE! `make export`
    #set(Libint2_DIR ${STAGED_INSTALL_PREFIX}/lib/cmake/libint2 CACHE PATH "path to internally built libint2-config.cmake" FORCE)
    # CHOOSE! pure cmake
    set(Libint2_DIR ${STAGED_INSTALL_PREFIX}/lib/cmake/libint2 CACHE PATH "path to internally built libint2-config.cmake" FORCE)

endif()


# NOTES Psi4 internal basis set library AM levels
#
# -------------          -------   ----
# Orbital Basis          1st Row   File
# -------------          -------   ----
# cc-pVDZ                   D2      F3
# cc-pVTZ                   F3      G4
# cc-pVQZ                   G4      H5
# cc-pV5Z                   H5      I6
# cc-pV6Z                   I6      I6
# lesser zapa follow Dunning
# 7zapa-nr                   7       7
# def2-svp*                 F3      F3
# def2-tzvp*                G4      G4
# def2-qzvp*                G4      G4
# -------------          -------   ----
#
# -------------          -------   ----
# Fitting Basis          1st Row   File
# -------------          -------   ----
# cc-pVDZ-RI                F3      G4
#        -JKFIT             F3      F3
# cc-pVTZ-RI                G4      H5
#        -JKFIT             G4      G4
# cc-pVQZ-RI                H5      I6
#        -JKFIT             H5      H5
# cc-pV5Z-RI                I6      I6
#        -JKFIT             I6      I6
# cc-pV6Z-RI                 7       7
#        -JKFIT             --      --
# def2-universal-jkfit      G4      I6
# sad-fitting               F3      F3
# -------------          -------   ----
