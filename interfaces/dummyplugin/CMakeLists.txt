include(ExternalProject)

set(PLUG "dummyplugin")
set(PLUG_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/${PLUG}")

ExternalProject_Add(${PLUG}
    PREFIX ${PLUG_PREFIX}
    GIT_REPOSITORY https://github.com/loriab/dummyplugin.git
    CONFIGURE_COMMAND cat ${PLUG_PREFIX}/src/${PLUG}/Makefile.template | sed -e "s|PLUGIN_CXX|${CMAKE_CXX_COMPILER}|g" -e "s|PLUGIN_DEFINES|${PLUGIN_DEFINES}|g" -e "s|PLUGIN_FLAGS|${CMAKE_CXX_FLAGS}|g" -e "s|PLUGIN_INCLUDES|${PLUGIN_INCLUDE_DIRECTORIES}|g" -e "s|PLUGIN_OBJDIR|${CMAKE_BINARY_DIR}|g" -e "s|PLUGIN_INSTALLDIR|${PLUG_PREFIX}|g" > ${PLUG_PREFIX}/src/${PLUG}/Makefile
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${MAKE}
)

