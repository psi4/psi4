#Maybe this tells it to do dynamic integral balancing?
add_definitions(-DDYNAMIC=1)

#Tells GTFock not build the SCF components
add_definitions(-DDO_SCF=0)

#Tells GTFock not to print SCF timings
add_definitions(-DSCF_TIMING=0)

#Tells GTFock (fuck if I know)
add_definitions(-DSCF_OUT=0)

#Tells GTFock the alignlen
add_definitions(-D__ALIGNLEN__=64)

#The additional flags that are wanted by GTFock
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")

#If we have an Intel compiler also add this one:
if(${CMAKE_C_COMPILER_ID} STREQUAL "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xHost -fasm-blocks -qno-offload")
endif() 

add_subdirectory(pfock)
add_subdirectory(libcint)
include_directories(${JK_GA_INC} ../libmints ../libqt pfock libcint)
add_library(GTFock STATIC ${PFOCKSRC} ${CINTSRC})
add_dependencies(GTFock GlobalArrays GTFockCInt)