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
