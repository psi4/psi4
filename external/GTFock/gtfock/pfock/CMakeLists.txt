set(GTFockHeaders ${CMAKE_CURRENT_SOURCE_DIR}/config.h 
                  ${CMAKE_CURRENT_SOURCE_DIR}/fock_buf.h 
                  ${CMAKE_CURRENT_SOURCE_DIR}/fock_task.h 
                  ${CMAKE_CURRENT_SOURCE_DIR}/pfock.h
                  ${CMAKE_CURRENT_SOURCE_DIR}/screening.h 
                  ${CMAKE_CURRENT_SOURCE_DIR}/taskq.h)
set(GTFockSRCs ${CMAKE_CURRENT_SOURCE_DIR}/fock_buf.c 
               ${CMAKE_CURRENT_SOURCE_DIR}/fock_task.c 
               ${CMAKE_CURRENT_SOURCE_DIR}/screening.c 
               ${CMAKE_CURRENT_SOURCE_DIR}/taskq.c 
               ${CMAKE_CURRENT_SOURCE_DIR}/pfock.c
)



set(PFOCKSRC ${GTFockHeaders} ${GTFockSRCs} PARENT_SCOPE)