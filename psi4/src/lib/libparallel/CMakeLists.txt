set(sources_list ParallelScanner.cc 
                 p_fprintf.cc 
                 PsiOutStream.cc 
                 BasesBase.cc 
                 process.cc 
                 ParallelPrinter.cc 
                 PsiInStream.cc 
)

add_library(parallel ${sources_list})
add_definitions(-DINSTALLEDPSIDATADIR="${CMAKE_INSTALL_PREFIX}")
target_include_directories(parallel PUBLIC options)
target_link_libraries(parallel options)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        DESTINATION psi4/src/lib/libparallel
        FILES_MATCHING PATTERN "*.h")
