add_custom_command(
    COMMAND tar xzf ${PROJECT_SOURCE_DIR}/boost/${BOOSTVER}.tar.gz 
    OUTPUT ${BOOSTVER}
    COMMENT Unpacking Boost
)

add_custom_command(
    COMMAND ./bootstrap.sh --with-libraries=filesystem,python,regex,serialization,system,thread --with-python=${PYTHON}
    COMMAND ./b2
    COMMAND touch ../boost.compiled
    WORKING_DIRECTORY ${BOOSTVER}
    DEPENDS ${BOOSTVER}
    OUTPUT boost.compiled
    COMMENT Configuring and building Boost 
)

# This target is always executed, but checks for the existence of boost.compiled
add_custom_target(
    boostlibs ALL
    DEPENDS boost.compiled
)
