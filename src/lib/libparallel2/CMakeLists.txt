set(LPSRC Algorithms.h
        Communicator.h
        Communicator.cc
        LibParallelBase.h
        LibParallelBase.cc
        LibParallelHelper.h
        LibParallelHelper.cc
        MPIJob.h
        MPITask.h
        ParallelEnvironment.h
        ParallelEnvironment.cc
)
add_subdirectory(CommEnvGuts)
add_subdirectory(Schedulers)
add_subdirectory(TaskJobGuts)
add_subdirectory(Util)
add_library(parallel2 STATIC ${LPSRC} ${CEGSRC} ${SSRC} ${TJGSRC} ${UTILSRC})
set_property(GLOBAL APPEND PROPERTY LIBLIST parallel2)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(parallel2 custom_boost)
endif()
