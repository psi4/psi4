set(SRC PsiMap.h
        PsiRdmNumGen.h
        Units.h
        Implementations/UnitsGuts.hh
        Implementations/PrefixConverter.cc
        Implementations/BaseUnitConverter.cc
)
add_library(PsiUtil ${SRC})
if(BUILD_CUSTOM_BOOST)
   add_dependencies(PsiUtil custom_boost)
endif()
