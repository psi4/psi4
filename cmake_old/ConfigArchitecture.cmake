set(LINUX_UBUNTU FALSE)
set(LINUX_FEDORA FALSE)

if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    add_definitions(-DSYS_LINUX)
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "i686")
        add_definitions(-DARCH32BIT)
    elseif(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64")
    else()
    endif()
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    add_definitions(-DSYS_DARWIN)
    add_definitions(-DVAR_MFDS)
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
    add_definitions(-DSYS_AIX)
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_definitions(-DSYS_WINDOWS)
endif()
