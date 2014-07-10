FILE(GLOB ERDSrc *.f *.F) 
add_library(erd STATIC ${ERDSrc})
install(TARGETS erd DESTINATION ${CMAKE_BINARY_DIR}/lib)
