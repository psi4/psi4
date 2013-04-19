add_custom_command(
    COMMENT Unpacking Boost
    COMMAND tar xzf ${PROJECT_SOURCE_DIR}/boost/${BOOSTVER}.tar.gz 
    OUTPUT ${BOOSTVER}
)

add_custom_command(
    COMMENT Configuring and building Boost 
    COMMAND ./bootstrap.sh --with-libraries=filesystem,python,regex,serialization,system,thread --with-python=${PYTHON}
    COMMAND ./b2
    COMMAND touch ../boost.compiled
    WORKING_DIRECTORY ${BOOSTVER}
    DEPENDS ${BOOSTVER}
    OUTPUT boost.compiled
)

# This target is always executed, but checks for the existence of boost.compiled
add_custom_target(
    boostlibs ALL
    DEPENDS boost.compiled
)
