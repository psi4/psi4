include(ExternalProject)

set(PLUG "gpu_dfcc")
set(PLUG_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/${PLUG}")

find_package(CUDA REQUIRED)
if(${CUDA_FOUND})
    message("-- NVCC=${CUDA_NVCC} version ${CUDA_VERSION_STRING} rootdir ${CUDA_TOOLKIT_ROOT_DIR} lib ${CUDA_LIBRARIES}")
else()
    message("nocuda")
endif()
set(MY_CUDAFLAGS "-arch sm_13 -Xcompiler -fopenmp --compiler-options '-fPIC' -O2")
set(MY_CUDALIBS "-lcublas -L/usr/local/cuda/lib64/")
set(MY_NVCC "nvcc")

ExternalProject_Add(gpu_dfcc
    PREFIX ${PLUG_PREFIX}
    GIT_REPOSITORY https://github.com/edeprince3/gpu_dfcc.git
    #CONFIGURE_COMMAND cat ${PLUG_PREFIX}/src/${PLUG}/Makefile.template | sed -e "s|MY_TOP_SRCDIR|${CMAKE_SOURCE_DIR}|g" -e "s|MY_TOP_OBJDIR|${CMAKE_BINARY_DIR}|g" -e "s|MY_NVCC|${CUDA_NVCC}|g" -e "s|MY_CUDAFLAGS|${MY_CUDAFLAGS}|g" -e "s|MY_CUDALIBS|${MY_CUDALIBS}|g" > ${PLUG_PREFIX}/src/${PLUG}/Makefile
    CONFIGURE_COMMAND cat ${PLUG_PREFIX}/src/${PLUG}/Makefile.template | sed -e "s|PLUGIN_CXX|${CMAKE_CXX_COMPILER}|g" -e "s|PLUGIN_DEFINES|${PLUGIN_DEFINES}|g" -e "s|PLUGIN_FLAGS|${CMAKE_CXX_FLAGS}|g" -e "s|PLUGIN_INCLUDES|${PLUGIN_INCLUDE_DIRECTORIES}|g" -e "s|PLUGIN_OBJDIR|${CMAKE_BINARY_DIR}|g" -e "s|PLUGIN_INSTALLDIR|${PLUG_PREFIX}|g" > ${PLUG_PREFIX}/src/${PLUG}/Makefile
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${MAKE}
)

