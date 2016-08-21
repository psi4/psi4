set(sources_list ParallelScanner.cc 
                 p_fprintf.cc 
                 PsiOutStream.cc 
                 BasesBase.cc 
                 process.cc 
                 ParallelPrinter.cc 
                 PsiInStream.cc 
)
message(STATUS "${CMAKE_INSTALL_PREFIX}")
add_definitions(-DINSTALLEDPSIDATADIR=${CMAKE_INSTALL_PREFIX}/share/psi4)
psi4_add_library(parallel sources_list options)
