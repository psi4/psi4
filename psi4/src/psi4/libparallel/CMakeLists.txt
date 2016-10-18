set(sources_list ParallelScanner.cc 
                 p_fprintf.cc 
                 PsiOutStream.cc 
                 BasesBase.cc 
                 process.cc 
                 ParallelPrinter.cc 
                 PsiInStream.cc 
)
psi4_add_module(lib parallel sources_list options)
