set(sources_list ParallelScanner.cc 
                 p_fprintf.cc 
                 PsiOutStream.cc 
                 BasesBase.cc 
                 process.cc 
                 ParallelPrinter.cc 
                 PsiInStream.cc 
)
add_definitions(-DINSTALLEDPSIDATADIR=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4)
psi4_add_module(lib parallel sources_list options)
