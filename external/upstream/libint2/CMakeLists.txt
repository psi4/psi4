# NOTES on Libint2 build options for Psi4
# ---------------------------------------
# * In order for Psi4 to at all work with a libint2 build, need LIBINT2_CARTGAUSS_ORDERING=standard *at generation time*. This is already the default.
# * In order for Psi4 to at all work with a libint2 build, need LIBINT2_SHELL_SET=standard *at generation time*. This is already the default.
# * ~In order for Psi4 to at all work with a libint2 build, need LIBINT2_SHGAUSS_ORDERING=gaussian *at library build time*.
#   This is not the default, though it may become runtime selectable in future GH:evaleev/libint/issues/190#issuecomment-691485225 .~ (see [Feb 2023] update)
# * In order for any DF to work, need ENABLE_ERI3 and ENABLE_ERI2 >=0
# * In order for gradient and Hessian tests to not segfault, need ENABLE_ERI and ENABLE_ERI3 and ENABLE_ERI2 =2
# * (new c. Jan 2022) In order for Psi4 to at all work with a libint2 build for one-electron ints, need ENABLE_ONEBODY =2
# * The `find_package(Libint2 ... COMPONENTS ...)` commands will check the above requirements met, provided the libint2
#   was built with cmake+cmake, not libtool+cmake.
# * Provided you know what you're doing wrt the above options and know what subset of Psi4 calculations you're targeting
#   (e.g., no Hessians or non-df only) and recognize that DF is everywhere (e.g., SAD guess and DF_SCF_GUESS) and may
#   need further suppression, it's fine to reduce or eliminate the `find_package` components.
# * ~Psi4 can run a libtool+cmake, but one needs to disable the `find_package` components check and correct for some
#   Libint2 vs. libint2 spelling differences. See "CHOOSE!" toggles here and in <psi4>/psi4/CMakeLists.txt~ (see [Feb 2023] update)
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
# [Feb 2023]
# * Upon additions to the Libint2 API to allow runtime SHGAUSS switching, Psi4 can use LIBINT2_SHGAUSS_ORDERING=standard|gaussian
#   *at Libint source generation time*. At Psi4 build time, it's going to be looking for standard (sss), so you'll have to
#   edit the CMake to use a gss.
# * Upon additions to the Libint2 API to allow runtime SHGAUSS switching, Psi4 can no longer be adapted to use a libtool+cmake source.
# * Libraries of different AM are not hot-swappable in a Psi4 build, but they only require a relinking.
# [Jan 2024]
# * As of v1.9, Psi4 can use an upstream v2.8 Libint2. This is a libtool+cmake build, so `find_package(Libint2 ... COMPONENTS ...)` has no components to check.
# [Apr 2024]
# * As of now, Psi4 needs v2.8.1 or higher for relativistic integral support.

find_package(Libint2 2.8.1 CONFIG  # edit in codedeps
#  COMPONENTS
#    sss
#    CXX_ho
#    impure_sh
#    "eri_c4_d0_l${MAX_AM_ERI}" eri_c3_d0_l4 eri_c2_d0_l4 onebody_d0_l4
#     eri_c4_d1_l2              eri_c3_d1_l3 eri_c2_d1_l3 onebody_d1_l3
#                                                         onebody_d2_l3
#  OPTIONAL_COMPONENTS
#     eri_c4_d2_l2              eri_c3_d2_l3 eri_c2_d2_l3
  )

# If L2 v2.8.1 (1) exported CMake components and (2) libint2-config.cmake actually
#   checked, the below is the new find_package call.
#if (${MAX_AM_ERI} GREATER 8)
#    set(amchar "llll")
#elseif (${MAX_AM_ERI} EQUAL 7)
#    set(amchar "kkkk")
#elseif (${MAX_AM_ERI} EQUAL 6)
#    set(amchar "iiii")
#elseif (${MAX_AM_ERI} EQUAL 5)
#    set(amchar "hhhh")
#elseif (${MAX_AM_ERI} LESS_EQUAL 4)
#    set(amchar "gggg")
#endif()
#find_package(
#  Libint2
#  2.8.1
#  CONFIG
#  COMPONENTS
#    ss
#    CXX_ho
#    impure_sh
#    "eri_${amchar}_d0" eri_ggg_d0 eri_gg_d0 onebody_g_d0
#     eri_dddd_d1       eri_fff_d1 eri_ff_d1 onebody_f_d1
#                                            onebody_f_d2
#  OPTIONAL_COMPONENTS
#     eri_dddd_d2       eri_fff_d2 eri_ff_d2
#  )

if (TARGET Libint2::cxx)
    get_target_property(_loc Libint2::int2 LOCATION)
    # defer until upstream provides full-dress targets get_target_property(Libint2_VERSION Libint2::int2 Libint2_VERSION)
    # defer until upstream provides full-dress targets get_target_property(Libint2_MAX_AM_ERI Libint2::int2 Libint2_MAX_AM_ERI)
    get_target_property(Libint2_CONFIGURATION Libint2::int2 Libint2_CONFIGURATION)
    message(STATUS "${Cyan}Found Libint2 ${Libint2_MAX_AM_ERI}${ColourReset}: ${_loc} (found version ${Libint2_VERSION})")
    # TODO revive this message. For now, can't tell if cmake build missing Hessian ints or libtool build w/o component info.
    # if (NOT eri_c4_d2_l2 IN_LIST Libint2_CONFIGURATION)
    #     message(WARNING "Libint2 detected but missing ERI 2nd derivative integrals (components eri_c4_d2_l2 eri_c3_d2_l3 eri_c2_d2_l3). Psi4 will automatically fall back to finite difference.")
    # endif()
    add_library(libint2_external INTERFACE)  # dummy

    if(${CMAKE_INSIST_FIND_PACKAGE_Libint2})
        message(VERBOSE "Suitable Libint2 located externally as user insists")
    endif()

else()
    if(${CMAKE_INSIST_FIND_PACKAGE_Libint2})
        message(FATAL_ERROR "Suitable Libint2 could not be externally located as user insists")
    endif()

    # stopgap to avert a report that static lib not working on Linux
    set(_build_shared_libs "ON")

    # reportedly Debug takes forever to build on Mac
    set(_cmake_build_type "Release")

    include(ExternalProject)

    if(NOT ${BUILD_Libint2_GENERATOR})
        if (${MAX_AM_ERI} GREATER 5)
            set(_url_am_src "7-7-4-12-7-5")  # passes E&G 7z conv and 6z df, H Qz conv and df. covers all internal basis sets for energy and gradient
        elseif (${MAX_AM_ERI} LESS_EQUAL 5)
            set(_url_am_src "5-4-3-6-5-4")  # passes full. production Mac, temp production Linux
        endif()

        set(_url_l2_tarball "https://github.com/loriab/libint/releases/download/v0.1/libint-2.8.1-${_url_am_src}_mm10f12ob2_0.tgz")

        message(STATUS "Suitable Libint2 could not be located, ${Magenta}Building Libint2 ${_url_am_src}${ColourReset} instead.")

        if (MSVC)
            # Windows shared (dll) can't work
            set(_build_shared_libs "OFF")
        endif()

        ExternalProject_Add(libint2_external
            URL ${_url_l2_tarball}
            CMAKE_ARGS -GNinja
                       -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${_cmake_build_type}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DBoost_DIR=${Boost_DIR}
                       -DBOOST_ROOT=${BOOST_ROOT}
                       -DEigen3_DIR=${Eigen3_DIR}
                       -DEigen3_ROOT=${Eigen3_ROOT}
                       -DBUILD_SHARED_LIBS=${_build_shared_libs}
                       -DREQUIRE_CXX_API=ON
                       -DENABLE_FORTRAN=OFF
                       -DLIBINT2_PYTHON=OFF
                       -DENABLE_XHOST=${ENABLE_XHOST}
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1}
                             -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1})

                       #-DLIBINT2_REQUIRE_CXX_API=ON
                       #-DLIBINT2_REQUIRE_CXX_API_COMPILED=OFF
                       #-DLIBINT2_ENABLE_FORTRAN=OFF
                       #-DLIBINT2_ENABLE_PYTHON=OFF
    else()
        message(STATUS "Suitable Libint2 could not be located, ${Magenta}Building Libint2${ColourReset} from generator source instead.")
        message(WARNING "${Yellow}Libint2 generator source build requires Ninja (not Makefile) builder and Boost, Eigen, and GMPXX dependencies. The build will be very long for non-trivial AM, and parallelism won't help for the generation phase. Edit external/upstream/libint2/CMakeLists.txt for AM, derivative, and integral classes settings. This is also building v2.7, not v2.8 Libint2.${ColourReset}")

        if (MSVC)
            message(FATAL_ERROR "${Red}Libint2 generator source build faulty on Windows.${ColourReset}")
        endif()

        # NOTE: active AM settings below will pass few tests, as they're set for 1st row triple-zeta ene and double-zeta ene/grad/hess for both conv/df. Please edit for desired use.

        ExternalProject_Add(libint2_external
            URL https://github.com/loriab/libint/archive/new-cmake-2023-take2-b.zip
            CMAKE_ARGS -GNinja
                       -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${_cmake_build_type}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DENABLE_ONEBODY=2
                       -DENABLE_ERI=2
                       -DENABLE_ERI3=2
                       -DENABLE_ERI2=2
                       -DENABLE_G12=1
                       -DDISABLE_ONEBODY_PROPERTY_DERIVS=ON
                       -DMULTIPOLE_MAX_ORDER=4
                       -DBoost_DIR=${Boost_DIR}
                       -DBOOST_ROOT=${BOOST_ROOT}
                       -DEigen3_DIR=${Eigen3_DIR}
                       -DEigen3_ROOT=${Eigen3_ROOT}
                       -DMultiprecision_ROOT=${Multiprecision_ROOT}
                       -DERI3_PURE_SH=OFF
                       -DERI2_PURE_SH=OFF
                       -DLIBINT2_SHGAUSS_ORDERING=standard
                       -DLIBINT2_CARTGAUSS_ORDERING=standard
                       -DLIBINT2_SHELL_SET=standard
                       -DBUILD_SHARED_LIBS=${_build_shared_libs}
                       -DLIBINT2_REQUIRE_CXX_API=ON
                       -DLIBINT2_REQUIRE_CXX_API_COMPILED=OFF
                       -DLIBINT2_ENABLE_FORTRAN=OFF
                       -DLIBINT2_ENABLE_PYTHON=OFF
                       -DBUILD_TESTING=${BUILD_TESTING}
                       -DENABLE_XHOST=${ENABLE_XHOST}
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1}
                             -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
                             -DWITH_ERI_MAX_AM:STRING=3;2;2
                             -DWITH_ERI3_MAX_AM:STRING=4;3;3
                             -DWITH_ERI2_MAX_AM:STRING=4;3;3
                             -DWITH_MAX_AM:STRING=4;3;3
                             -DWITH_G12_MAX_AM:STRING=3
        )
    endif()

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
