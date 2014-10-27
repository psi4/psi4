set(SRC StreamBase.h
        PsiInStream.h
        PsiInStream.cc
        PsiOutStream.h
        PsiOutStream.cc
        process.cc 
        ParallelPrinter.h
        ParallelPrinter.cc
        BasesBase.h
        BasesBase.cc
        ParallelScanner.h
        ParallelScanner.cc
        PsiFileImpl.h 
        TableSpecs.h
)
if(${MPI_FOUND})
   set(SRC ${SRC} mpi_wrapper.cc)
endif()
add_library(parallel ${SRC})
add_dependencies(parallel mints)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(parallel custom_boost)
endif()   
