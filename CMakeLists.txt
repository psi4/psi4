#       This is the top-level CMakeLists.txt for Psi4.  
# 
# This is version 2.0.0 of the build system.  Version 1 was a swell first
# pass, but quickly became too complicated.  For projects of this size, with
# lots of little (soon to be somewhat) independent pieces, the general CMake
# strategy is known as a superbuild.  In such a build all pieces of a package
# are treated as external projects, including the main project itself.  What
# this does is allow us to use the same mechanism (namely find_package) to
# find dependencies and configure each part.
#
# Big changes:
#   - Debug build is warning free with GNU compilers.  Keep it that way...
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
#       directory and were thus accessible like: <header.h> vs. <lib/header.h>
#     - The reason for this is we want to transition to a modular nature and
#       it is easiest for the users of these modules if they only have to
#       specify the root dir as an include and not each subdir.  Note that this
#       is how nearly every professional C++ library I am aware of does this 
#       (Boost, Eigen, TBB, ...)
#     - The second (and arguably better reason) is doing it any other way makes
#       the superbuild a lot harder when something is staged.
#   - I have attempted to purge convenience headers (i.e. mints.h)
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
#        - Splitting the header is a nightmare if code gets factored, which
#          really should happen to mints...
#        - I realize that getting rid of mints.h means you need to include a lot
#          of header files, but I'll argue this is a surefire sign that APIs are
#          not designed well and need reexamined (e.g. do you really need the
#          whole object, or would just a matrix suffice?).
#   - I have pulled most of the forward declarations out of libmints/tyepdefs.h.
#     I saw many cases where people were forward declaring them again.  The ones
#     that remain are those required for the actual typedefs.  There really
#     isn't anything bad about this code wise, but it makes us look like we 
#     don't understand our own code.
#   - psiconfig.h is gone.  The macros in it are easily replace by CMake's
#     add_definitions function and doing so clears up the macro namespace
#
# Small changes:
#    - Static vs. shared libraries can be toggled at compile time
#    - Time being: cut out setup.py it just wrapped the CMake call and did so
#      in a manner inconsistent with usual CMake procedures.  Now veteran CMake
#      users can just go right for the throat.
#    - The vast majority of includes that use angled-brackets have been changed
#      to quotes.  Quotes indicate that a header file is part of the current
#      project and angled-brackets indicate it is an external dependency.  This
#      matters to the compiler which uses this to determine the look-up order.
#    - exception.h has been moved into libpsi4util to be with its source file
#    - process.h   has been moved into libparallel to be with its source file
#    - boost::shared_ptr is forwarded declared even when actual instances
#      are made.  A forward declaration only gets you something if you have
#      pointers or references to that type.  Your code has been compiling
#      because somewhere up the include tree someone included boost::shared_ptr
#
#    General Comments:
#    - A lot of files have an artifical dependence on boost::python through
#      liboptions.h, via process.h, via psi4-dec.h.  Eliminating process.h will
#      fix this and then dependencies on boost::python should be reevaluated
#    - The header files in libmints need at least another pass to solidify the
#      actual dependencies
#    - Conventions need established.  *.h is usually used for header files in
#      Psi4, but *.hpp appears as well (*.h is by convention for C header files
#      so I recommend switching to *.hpp)
#    - Admittedly I am one of the worst culprits, but case semantics need to be
#      established as well.  Are we doing snake_case, CamelCase, etc.




# Useful CMake variables you should consider setting:
#    - CMAKE_C_COMPILER The C compiler
#    - CMAKE_CXX_COMPILER The C++ compiler
#    - CMAKE_C_FLAGS Any extra C flags you want to pass
#    - CMAKE_CXX_FLAGS Any extra C++ flags you want to pass
#    - CMAKE_BINARY_DIR The build directory (default: build)
#    - CMAKE_INSTALL_PREFIX The install directory (default: /usr/local)
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
include(cmake/Psi4Macros.cmake)
no_in_source()

############################## Options #########################################
option(ENABLE_MPI "Enables modules that depend on MPI (currently GTFock)" OFF)
option(ENABLE_PCMSOLVER "Enbales the PCM solver module" OFF)
option(ENABLE_GTFOCK "Enables GTFock" OFF)
option(BUILD_SHARED_LIBS "Psi4 libraries will be shared objects" OFF)
option_with_default(LIBINT_OPT_AM "Maximum angular momentum for integrals" 5)
option_with_default(CMAKE_BUILD_TYPE "The type of our build" Debug)#Change to release

############################## Dependencies ####################################
find_package(PythonInterp)
find_package(PythonLibs)
add_subdirectory(external)

list(APPEND CMAKE_PREFIX_PATH 
     ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/external)

############################## Features ########################################
list(APPEND CMAKE_CXX_FLAGS "-std=c++11")
test_restrict(RESTRICT_KEYWORD)

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
