if(ENABLE_ASAN)
    if((CMAKE_CXX_COMPILER_ID MATCHES Clang) OR
        (CMAKE_CXX_COMPILER_ID MATCHES AppleClang) OR
        (CMAKE_CXX_COMPILER_ID MATCHES GNU))
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    else()
        message(WARNING "ASAN flags are not known for your compiler ${CMAKE_CXX_COMPILER_ID}")
    endif()
endif()

if(ENABLE_TSAN)
    if((CMAKE_CXX_COMPILER_ID MATCHES Clang) OR
        (CMAKE_CXX_COMPILER_ID MATCHES AppleClang) OR
        (CMAKE_CXX_COMPILER_ID MATCHES GNU))
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread -fno-omit-frame-pointer -pie")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=thread -fno-omit-frame-pointer -pie")
    else()
        message(WARNING "TSAN flags are not known for your compiler ${CMAKE_CXX_COMPILER_ID}")
    endif()
endif()

if(ENABLE_UBSAN)
    if((CMAKE_CXX_COMPILER_ID MATCHES Clang) OR
        (CMAKE_CXX_COMPILER_ID MATCHES AppleClang) OR
        (CMAKE_CXX_COMPILER_ID MATCHES GNU))
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined -fno-omit-frame-pointer")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=undefined -fno-omit-frame-pointer")
    else()
        message(WARNING "UBSAN flags are not known for your compiler ${CMAKE_CXX_COMPILER_ID}")
    endif()
endif()

