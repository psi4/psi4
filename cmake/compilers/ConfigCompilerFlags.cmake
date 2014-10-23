include(SaveCompilerFlags)
include(FindBrokenIntrinsics)
include(OptimizeForArchitecture)
include(TestRestrict)

if(ENABLE_VECTORIZATION)
	set(DEFINITIONS)
	set(CXX_ARCHITECTURE_FLAGS)
	set(C_ARCHITECTURE_FLAGS)
	set(Fortran_ARCHITECTURE_FLAGS)
	find_broken_intrinsics()
	OptimizeForArchitecture()
	message(STATUS "Vectorization flags used:")
	message(STATUS "  CXX:     ${CXX_ARCHITECTURE_FLAGS}")
	message(STATUS "  C:       ${C_ARCHITECTURE_FLAGS}")
	message(STATUS "  Fortran: ${Fortran_ARCHITECTURE_FLAGS}")
endif()	

# Test for restrict keyword
set(restrict_keyword "")
test_restrict(restrict_keyword)

if(CMAKE_C_COMPILER_WORKS)
    include(CFlags)
endif()

if(CMAKE_CXX_COMPILER_WORKS)
    include(CXXFlags)
endif()

if(CMAKE_Fortran_COMPILER_WORKS)
    include(FortranFlags)
endif()
