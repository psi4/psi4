set(sources_list ParallelScanner.cc 
                 p_fprintf.cc 
                 PsiOutStream.cc 
                 BasesBase.cc 
                 process.cc 
                 ParallelPrinter.cc 
                 PsiInStream.cc 
)
add_definitions(-DINSTALLEDPSIDATADIR="${CMAKE_INSTALL_PREFIX}")
psi4_add_library(parallel sources_list)
target_include_directories(parallel PUBLIC options)
target_link_libraries(parallel options)