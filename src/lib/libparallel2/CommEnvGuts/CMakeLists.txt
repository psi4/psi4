set(CEGSRC 
    ${CMAKE_CURRENT_SOURCE_DIR}/CommBase.h
    ${CMAKE_CURRENT_SOURCE_DIR}/CommGuts.h
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalComm.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ParallelEnvironmentGuts.h
    ${CMAKE_CURRENT_SOURCE_DIR}/CommGuts.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalComm.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ParallelEnvironmentGuts.cc
)
if(${MPI_FOUND})
   set(CEGSRC ${CEGSRC}
      ${CMAKE_CURRENT_SOURCE_DIR}/ParallelComm.h
      ${CMAKE_CURRENT_SOURCE_DIR}/ParallelComm.cc
   )
endif()
set(CEGSRC ${CEGSRC} PARENT_SCOPE)
