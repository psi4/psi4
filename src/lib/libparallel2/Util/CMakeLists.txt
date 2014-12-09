set(UTILSRC ${CMAKE_CURRENT_SOURCE_DIR}/MasterTools.h
            ${CMAKE_CURRENT_SOURCE_DIR}/MasterTools.cc
            ${CMAKE_CURRENT_SOURCE_DIR}/MPISetterUpper.h
            ${CMAKE_CURRENT_SOURCE_DIR}/MPISetterUpper.cc
            ${CMAKE_CURRENT_SOURCE_DIR}/MPITaskQueue.h
            ${CMAKE_CURRENT_SOURCE_DIR}/MPITaskQueue.cc
            ${CMAKE_CURRENT_SOURCE_DIR}/SlaveTools.h
            ${CMAKE_CURRENT_SOURCE_DIR}/SlaveTools.cc
            ${CMAKE_CURRENT_SOURCE_DIR}/TaskMap.h
            ${CMAKE_CURRENT_SOURCE_DIR}/TaskMap.cc
            ${CMAKE_CURRENT_SOURCE_DIR}/TaskStatistics.h
            ${CMAKE_CURRENT_SOURCE_DIR}/TaskStatistics.cc
            ${CMAKE_CURRENT_SOURCE_DIR}/ToolSet.h
            ${CMAKE_CURRENT_SOURCE_DIR}/ToolSet.cc
            ${CMAKE_CURRENT_SOURCE_DIR}/UpCastWrapper.h
)
set(UTILSRC ${UTILSRC} PARENT_SCOPE)