set(SSRC ${CMAKE_CURRENT_SOURCE_DIR}/DynamicRoundRobin.h
         ${CMAKE_CURRENT_SOURCE_DIR}/DynamicRoundRobin.cc
         ${CMAKE_CURRENT_SOURCE_DIR}/DynamicScheduler.h
         ${CMAKE_CURRENT_SOURCE_DIR}/DynamicScheduler.cc
         ${CMAKE_CURRENT_SOURCE_DIR}/MPIScheduler.h
         ${CMAKE_CURRENT_SOURCE_DIR}/MPIScheduler.cc
         ${CMAKE_CURRENT_SOURCE_DIR}/MPIStaticScheduler.h
         ${CMAKE_CURRENT_SOURCE_DIR}/MPIStaticScheduler.cc
         ${CMAKE_CURRENT_SOURCE_DIR}/SimpleScheduler.h
         ${CMAKE_CURRENT_SOURCE_DIR}/SimpleScheduler.cc
         ${CMAKE_CURRENT_SOURCE_DIR}/StaticRoundRobin.h
         ${CMAKE_CURRENT_SOURCE_DIR}/StaticRoundRobin.cc
)
set(SSRC ${SSRC} PARENT_SCOPE)