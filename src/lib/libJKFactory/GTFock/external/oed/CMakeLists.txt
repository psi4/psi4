FILE(GLOB OEDSrc *.f *.F)
add_library(oed STATIC ${OEDSrc})
install(TARGETS oed DESTINATION ${CMAKE_BINARY_DIR}/lib)