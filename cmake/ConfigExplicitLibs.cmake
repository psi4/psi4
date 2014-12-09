set(EXPLICIT_LIBS
    "${EXPLICIT_LIBS}"
    CACHE STRING
    "User set math libraries"
    FORCE
    )

if(EXPLICIT_LIBS)
    set(EXTERNAL_LIBS
        ${EXTERNAL_LIBS}
        ${EXPLICIT_LIBS}
        )
    message("-- User set explicit libraries: ${EXPLICIT_LIBS}")
endif()

if(ENABLE_ACCELERATE)
   message(STATUS "Using Mac OS X Accelerate Framework")
   set(EXTERNAL_LIBS ${EXTERNAL_LIBS} "-framework accelerate")
   # We now set BLAS_FOUND and LAPACK_FOUND to TRUE
   set(BLAS_FOUND TRUE)
   set(LAPACK_FOUND TRUE)
endif()
