#       This is the top-level CMakeLists.txt for Psi4.  
# 
# This is version 2.0.0 of the build system.  Version 1 was a swell first
# pass, but quickly became too complicated.  This is a purging of that system
# and hopefully should be a lot cleaner and easier to extend.  Please try to
# keep it clean.  If you don't know what you are doing asking for help.
#
# Quick overview of how the new build works.
# For projects of this size, with lots of little (hopefully soon) independent 
# pieces, the general CMake strategy is known as a superbuild.  In such a build 
# all pieces of a package, including the dependencies, optional modules, and 
# the main project itself are treated as external projects.  What
# this does is allow us to use the same mechanism (namely find_package) to
# find dependencies and configure each part.  The main psi4 source lives in psi4
# whereas optional modules or dependencies live in external.  When trying to
# figure out where something goes, basically, if you personally didn't write the
# code, or it is intended to be used with many other packages, it probably falls
# under the category of going in the external folder.  Ideally the interfaces
# to these external codes will migrate to external as well, but I haven't done 
# that yet.
#
# The central idea is now that each external project (including the main one)
# should be capable of standing on its own.  Thus internally they should call
# find_package for every library they want.  It is the responsibility of the
# projects above it to set paths/flags appropriately so that find_package finds
# what it is looking for.  This also means that if you need access to a library
# that CMake doesn't know how to find you need to write a FindLibrary.cmake
# file for it or have the one in the external folder create a 
# LibraryConf.cmake file (neither file is terribly hard to make).
#
# The other common theme to this rewrite is cleanliness.  I have attempted to
# remove extraneous commands and have introduced macros for common operations.
# Much if/else logic has been replaced by simply having NULL variables and NULL
# targets.  For example, we mark Psi4 as always depending on all contents of the
# external folder; if those external targets aren't built they show up as NULL
# targets (targets which don't change the build in the slightest, but allow us
# to avoid if/else logic).  If they are built then their libraries are
# appended to the PSI4_MODULES variable.
#
#
# Big changes:
#   - Libraries that are not part of Psi4 proper are now in external
#     - Although allowed by the more permissive licenses, it is somewhat dubious
#       to mix external libraries in with ones own source code.  Separating them
#       makes it clear to others that these are not our codes, but used with
#       permission
#     - This right now applies to: boost, erd, libint, efp, dkh, and gdma
#     - The primary advantage to this is that we can keep our source code free
#       of the dependencies of codes we interface to (i.e. we no longer depend
#       on Fortran)
#   - All include files are relative to root
#     - In practice my find and replace missed headers that are in the same
#       directory and were thus accessible like: header.h vs. lib/header.h
#     - The reason for this is we want to transition to a modular nature and
#       it is easiest for the users of these modules if they only have to
#       specify the root dir as an include and not each subdir.  Note that this
#       is how nearly every professional C++ library I am aware of does this 
#       (Boost, Eigen, TBB, ...)
#     - The second (and arguably better reason) is doing it any other way makes
#       the superbuild a lot harder when something is staged.
#   - I have purged mints.h
#     - Convenience headers are really for users using your library
#     - Using them as a developer is evil for several reasons:
#        - Longer compile times, the reason you have to recompile half of Psi4
#          when you change something in libmints is partly because of mints.h
#          (the other part is every library was marked as depending on mints
#           even though they only need to be linked against it)
#        - Harder to know exactly what you brought in.  The vast majority of
#          people have no idea what is in mints.h or typedef.h as seen by their
#          redoing of work done by those headers
#        - Can destroy distinction between interface/implementation classes if
#          done incorrectly
#        - Splitting the header is a nightmare (as I can attest to) and is 
#          required if code gets factored, which really should happen to mints.
#   - I have pulled most of the forward declarations out of libmints/tyepdefs.h.
#     I saw many cases where people were forward declaring them again.  The ones
#     that remain are those required for the actual typedefs.  There really
#     isn't anything bad about this code wise, but it makes us look like we 
#     don't understand our own code.  Forward declarations are used incorrectly
#     all over the code and some clean-up is really needed.
#   - psiconfig.h is gone.  The macros in it are easily replace by CMake's
#     add_definitions function and doing so clears up the macro namespace
#
# Small changes:
#    - Building static vs. shared libraries can be toggled at compile time
#    - For the time being: cut out setup.py it just wrapped the CMake call My
#      understanding is this was meant to ease the transition to CMake from
#      autoconf, well we're transitioned
#    - The vast majority of includes that incorrectly used angled-brackets have 
#      been changed to quotes.  Although not exactly correct, for all
#      intents and purposes, use quotes to indicate that a header file is part 
#      of the current project and angled-brackets indicate it is an external 
#      dependency.  This is not just semantics, it changes how the compiler
#      behaves
#    - exception.h has been moved into libpsi4util to be with its source file
#    - process.h   has been moved into libparallel to be with its source file
#    - There is a lot of copy/paste going on in the preamble of source files
#      much of which makes you look like you don't understand C++ (e.g.
#      using namespace psi; namespace psi{}. You need one or the other...)
#      I have cleaned some of it up, but plenty remains.
#    - Most ConfigX files are gone.  From what I can tell, these files are what
#      should be in the CMakeLists.txt wrapping project X.  That's where
#      they now reside.
#    - Some options have moved from this top-level CMakeLists.txt:
#      - ENABLE_BOUNDS_CHECK is a Fortran thing and now appears before external
#        Fortran projects
#    - Some options are now gone
#      - ENABLE_MPI Psi4 proper contains no MPI now, if external projects want
#                   to include them that's fine
#      - ENABLE_OMP I never understood this one, if you don't want openMP for
#                   some reason set OpenMP_FOUND to FALSE
#      - ENABLE_AUTO_BLAS this was never used
#      - ENABLE_AUTO_LAPACK again never used
#      - ENABLE_ACCELERATE FindLAPACK knows about accelerate, if you want to use
#                          append the path to LAPACK_DIR
#      - ENABLE_CSR Wasn't used
#      - ENABLE_SCALAPACK again no MPI so don't need this
#      - ENABLE_SCALASCA  no scalapack no scalapack profiling
#      - ENABLE_CXX11_SUPPORT is not an option anymore
#      - ENABLE_DUMMY_PLUGIN not sure what this was supposed to do
#    - Other options that are in limbo
#      - ENABLE_UNIT_TESTS this was for turning on Boost's unit tests, if we
#        want it back put it in external/boost/CMakeLists.txt
#      - ENABLE_GPU_DFCC I don't think this worked, but in theory we want it to
#      - ENABLE_AMBIT Right now nothing is using it, if we are now requiring it
#        (I think that's the result of the tensor meeting) should this option go
#        the way of CXX11_SUPPORT?
#    - ConfigExplicitLibs is gone as it seemed tied to ENABLE_ACCELERATE
#    - ConfigExternal is gone, I think that was a poor man's attempt at
#      emulating ExternalProject_Add
#    - ConfigMPI is gone because again there is no more MPI
#    - ConfigOMP is gone because it was just a wrapper basically (Intel claims
#      there is no need to separately specify MKL and OMP threads as MKL will
#      determine if it is in a parallel region and not thread itself.  Google
#      mkl threads and openmp and see top hit).
#    - ConfigPython is gone in favor of using CMake's find_package
#    - The compiler/PerTargetCompilerFlags file is gone as it contained flags
#      for libint/libderiv and should be handled there
#    - I am somewhat uneasy with the the cmake/compilers cmake files, but don't
#      understand what the various flags do to weigh in.  The logic seems overly
#      complicated and at least on my box does next to nothing.  They are still
#      there though.  I'd like to see most of that go through the
#      optional_add_flags macro in Psi4Macros.cmake
#    - Targets no longer depend on header files.  CMake uses the compiler
#      generated dependency information to track changes to them.  Adding them
#      clutters the build
#    - src/bin/psi4 is now src/bin/psi4_main so that the target it builds is
#      the same as the directory name
#    General Comments:
#    - I updated the best practices wiki on good header file practices.  I 
#      encourage all users to read it so that hopefully our header files never
#      get so bad again.
#    - The header files in libmints need at least another pass to solidify the
#      actual dependencies
#    - Conventions need established.  *.h is usually used for header files in
#      Psi4, but *.hpp appears as well (*.h is by convention for C header files
#      so I recommend switching to *.hpp)
#    - Admittedly I am one of the worst culprits, but case semantics for file
#      names (and eventually the code itself) need to be established.  There
#      are some situations where snake_case is used, others where CamelCase is
#      used, and others wherethisisused (where this is used).

# Useful CMake variables you should consider setting:
#    - CMAKE_C_COMPILER The C compiler
#    - CMAKE_CXX_COMPILER The C++ compiler
#    - Fortran_COMPILER The Fortran compiler (required for some external projs)
#    - CMAKE_C_FLAGS Any extra C flags you want to pass
#    - CMAKE_CXX_FLAGS Any extra C++ flags you want to pass
#    - Fortran_FLAGS Any extra Fortran flags you want to pass
#    - CMAKE_BINARY_DIR The build directory (default: build)
#    - CMAKE_INSTALL_PREFIX The install directory (default: /usr/local)
#    - CMAKE_BUILD_TYPE Is this Debug, Release, or Profiling
#    - BUILD_SHARED_LIBS If ON Psi4 libraries will be shared
#    - PYTHON_EXECUTABLE Path to the python executable you want to use
#    - PYTHON_LIBRARY Path to the python library that goes with the exe
#    - PYTHON_INCLUDE_DIR Path to the python include files    
#
#

#Need 3.2 for build_always
cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(psi4 CXX C)
include(ExternalProject)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake
                              ${CMAKE_SOURCE_DIR}/cmake/compilers
)
include(Psi4Macros)
no_in_source()

############################## Options #########################################
option(ENABLE_PCMSOLVER "Enbales the PCM solver module" OFF)
option(ENABLE_GTFOCK "Enables GTFock (requires MPI)" OFF)
option(ENABLE_CHEMPS2 "Enables chemps2 for DMRG" OFF)
option(ENABLE_DKH "Enables DKH integrals" OFF)
option(ENABLE_ERD "Enables ERD instead of libint" OFF)
option(ENABLE_GDMA "Enables Stone's GDMA multipole code" OFF)
option(ENABLE_AMBIT "Enables the Ambit tensor library" OFF)
option(BUILD_SHARED_LIBS "Psi4 libraries will be shared objects" OFF)
option(ENABLE_XHOST "Enables processor-specific optimization")
option(ENABLE_CODE_COVERAGE "Enables details on code coverage")
option(ENABLE_ASAN "Enables address sanitizer")
option(ENABLE_TSAN "Enables thread sanitizer")
option(ENABLE_UBSAN "Enables undefined behavior sanitizer" OFF)
option_with_default(LIBINT_OPT_AM "Maximum angular momentum for integrals" 5)
option_with_default(CMAKE_BUILD_TYPE "The type of our build" Debug)#Change to release

############################## Compiler Features ###############################
set(ENABLE_CXX11_SUPPORT ON)
include(ConfigArchitecture)
include(ConfigCompilerFlags)
find_package(OpenMP)
optional_add_flags(OpenMP_FOUND ${OpenMP_CXX_FLAGS})

############################## Dependencies ####################################
find_package(PythonInterp)
find_package(PythonLibs)
add_subdirectory(external)

set(PSI4_STAGE_EX ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/external)
list(APPEND CMAKE_PREFIX_PATH  ${PSI4_STAGE_EX})

############################## Main Project ####################################

ExternalProject_Add(psi4_core
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/psi4
   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
              -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
              -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
              -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
              -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
              -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
              -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
              -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
              -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
              -DPYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR}
              -DPYTHON_LIBRARY=${PYTHON_LIBRARY}
              -DPSI4_ROOT=${CMAKE_CURRENT_SOURCE_DIR}
              -DBOOST_ROOT=${PSI4_STAGE_EX}/boost
              -DRESTRICT_KEYWORD=${RESTRICT_KEYWORD}
   BUILD_ALWAYS 1
   INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install DESTDIR=${CMAKE_BINARY_DIR}/stage
)
add_dependencies(psi4_core boost_external 
                           libint_external
                           libderiv_external 
                           libefp_external
)
install(DIRECTORY ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}
        DESTINATION ${CMAKE_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS)
include(ConfigInfo)
