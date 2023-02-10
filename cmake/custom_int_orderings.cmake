# from https://github.com/loriab/libint/blob/new-cmake-2023-take2-b/cmake/modules/int_orderings.cmake
# hopefully someday from https://github.com/evaleev/libint/blob/master/cmake/modules/int_orderings.cmake

# <<<  solid harmonic Gaussian orderings  >>>

if (psi4_SHGAUSS_ORDERING STREQUAL "standard")
    set(psi4_SHGSHELL_ORDERING 1)
elseif (psi4_SHGAUSS_ORDERING STREQUAL "gaussian")
    set(psi4_SHGSHELL_ORDERING 2)
else()
    message(FATAL_ERROR "Invalid value for psi4_SHGAUSS_ORDERING (${psi4_SHGAUSS_ORDERING})")
endif()

