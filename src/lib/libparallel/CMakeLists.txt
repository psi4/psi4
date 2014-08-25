set(SRC StreamBase.h
        StreamBase.cc
        process.cc 
        ParallelPrinter.h
        ParallelPrinter.cc
        TableSpecs.h
)
if(${MPI_FOUND})
   set(SRC ${SRC} mpi_wrapper.cc)
endif()
add_library(parallel ${SRC})
add_dependencies(parallel mints)
