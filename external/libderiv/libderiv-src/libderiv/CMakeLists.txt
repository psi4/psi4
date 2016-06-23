include(ExternalProject)

ExternalProject_Add(libderiv_compiler
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libderiv_compiler
   CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
              -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
              -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
              -DLIBINT_OPT_AM=${LIBINT_OPT_AM}
              -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
              -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
   INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install DESTDIR=${CMAKE_BINARY_DIR}/stage
)

#We generalize download to "getting the source", which means generate it here
ExternalProject_Add(libderiv_library
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libderiv_library
    CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
               -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
               -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
    DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/src
    DOWNLOAD_COMMAND ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/libderiv_compiler
    LOG_DOWNLOAD 1
    BINARY_DIR ${CMAKE_BINARY_DIR}/src
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install DESTDIR=${CMAKE_BINARY_DIR}/stage
)

add_dependencies(libderiv_library libderiv_compiler)
configure_file("libderivConfig.cmake.in" "libderivConfig.cmake" @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/libderiv/libderivConfig.cmake 
        DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/lib
        DESTINATION ${CMAKE_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS)
install(DIRECTORY ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/include
        DESTINATION ${CMAKE_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS)
