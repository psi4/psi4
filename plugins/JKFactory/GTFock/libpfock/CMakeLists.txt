set(LIBPFOCK_SRCS config.h
                  fock_task.h
                  one_electron.h
                  pfock_def.h
                  pfock_type.h
                  screening.h
                  taskq.h
                  utils.h
                  fock_task.c
                  one_electron.c
                  pfock.c
                  screening.c
                  taskq.c
                  utils.c)
include_directories(../include)
add_library(pfock STATIC ${LIBPFOCK_SRCS})
install(TARGETS pfock DESTINATION ${CMAKE_BINARY_DIR}/lib)
