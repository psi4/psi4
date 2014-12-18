set(CINTHeaders ${CMAKE_CURRENT_SOURCE_DIR}/basisset.h 
                ${CMAKE_CURRENT_SOURCE_DIR}/CInt.h)
set(CINTSrcs ${CMAKE_CURRENT_SOURCE_DIR}/basisset.c)
set(CINTSRC ${CINTHeaders} ${CINTSrcs} PARENT_SCOPE)