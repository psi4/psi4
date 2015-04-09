set(OBABELVER openbabel_2_3_9)

include(ExternalProject)
#After extracting openbabel.tar.gz a directory called OpenBabel is made
#that directory is where the top-level CMakeLists.txt lives.
#ExternalProject puts us in a directory at the same level as this OpenBabel
#directory for build, thus the configure command is:
#     cmake ../OpenBabel
#to build openbabel

ExternalProject_Add(OpenBabel
    URL ${PROJECT_SOURCE_DIR}/openbabel/${OBABELVER}.tar.gz
    CONFIGURE_COMMAND ${CMAKE_COMMAND} ../OpenBabel/ -DRUN_SWIG=ON -DPYTHON_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} 
    BUILD_COMMAND ${make} 
    INSTALL_DIR ${PROJECT_BINARY_DIR}
    INSTALL_COMMAND make install
)
