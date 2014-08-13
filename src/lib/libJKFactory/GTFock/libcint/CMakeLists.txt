set(CINT_SRCS basisset.h
              cint_def.h
              cint_type.h
              config.h
              erd_integral.h
              #oed_integral.h
              basisset.c
              erd_integral.c
              oed_integral.c)
include_directories(../include)
add_library(cint STATIC ${CINT_SRCS})
install(TARGETS cint DESTINATION ${CMAKE_BINARY_DIR}/lib)
