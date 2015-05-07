set(SRC ProgressBar.h
        PsiMap.h
        PsiRdmNumGen.h
        PythonFxn.h
        Table.h
        Timer.h
        Units.h
        AtomicData.h
        Exception2.h
        Implementations/UnitsGuts.hh
        Implementations/PrefixConverter.cc
        Implementations/Timer.cc
        Implementations/BaseUnitConverter.cc
        Implementations/PythonFxnGuts.cc
        Implementations/AtomicData.cc
        Implementations/ProgressBar.cc
        Implementations/Exception2.cc
)
add_library(PsiUtil ${SRC})
add_dependencies(PsiUtil parallel2)
include_directories(libparallel2)
set_property(GLOBAL APPEND PROPERTY LIBLIST PsiUtil)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(PsiUtil custom_boost)
endif()
