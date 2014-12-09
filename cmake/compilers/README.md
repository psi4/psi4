
Supported compilers
===================

The supported compilers are GNU, Intel and Clang. A `FATAL_ERROR` will be issued if any compiler vendor
other than these three is provided. No checks on compiler version are performed, apart from the checks
on Intel compiler version in the C++11 support detection (these might one day go in favour of a better
approach)
Doing checks on compiler versions is not really portable, even using CMake. The best way to do it
would be to issue `try_compile` commands on some test source code and enable/disable compiler flags
and/or compiler versions based on the outcome of the tests.
This would raise the complexity of the compiler settings significantly and it's thus currently
avoided.

The compiler vendor `FATAL_ERROR` can be bypassed by passing custom flags for C++, C and Fortran (_vide infra_)
There is no guarantee that a successful build will result from this or that a sane executable
will be produced.

Compiler flags setup
====================

There are three ways of setting compiler flags for the whole project and the plugins:

1. by specifying the build type `--type=[{release,debug,profile}]` or `-DCMAKE_BUILD_TYPE={Release,Debug,Profile}`. 
This sets some sensible defaults for the three cases. Only profiling with `gprof` is currently supported. 
A profiling build type is just a release build type with the addition of the flags for linking against `gprof`. 
The default flags are defined in the files `CXXFlags.cmake`, `CFlags.cmake` and `FortranFlags.cmake` in `cmake/compilers`
2. by specifying the build type and some extra flags using `--extra-cc-flags`, `--extra-cxx-flags` and `--extra-fortran-flags` (or the corresponding -D CMake option) 
These flags are appended to the default ones. **No checks** are performed on the compatibility with the default flags.
3. completely bypass defaults and define your own compiler flags. Use the options `--custom-cc-flags`, `--custom-cxx-flags` and `--custom-fortran-flags`. 
**No checks** are performed on the suitability of these flags.

Globally available compiler flags sets
======================================

Given a language `LANG`, the compiler flags sets available **globally** within the project are:

1. `CMAKE_<LANG>_FLAGS` the basic flags. These are always set.
2. `CMAKE_<LANG>_FLAGS_DEBUG` the debug flags;
3. `CMAKE_<LANG>_FLAGS_RELEASE` the release flags;
4. `CMAKE_<LANG>_FLAGS_PROFILE` the profiling flags;

Additional flags are appended to the `CMAKE_<LANG>_FLAGS` if a coverage analysis (using `gcov`) is to be
performed.
Additional flags are appended to the `CMAKE_<LANG>_FLAGS_DEBUG` if sanitizers are used for dynamic analysis.
**WARNING** The use of sanitizers is currently only available with Clang 3.6 Neither the `setup` script nor CMake
perform tests to check that the sanitizers options are relevant for the chosen compiler!
