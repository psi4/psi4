set(libfragsrc GMBE.h
               LibFragDriver.h
			   LibFragDriver.cc
			   MBE.h
			   MBEProp.h
)
add_subdirectory(libmolecule)
add_library(frag STATIC ${libfragsrc} ${MOLSRC})  
include_directories(../libPsiUtil ../libqt libmolecule)
add_dependencies(frag parallel2 PsiUtil)
target_link_libraries(frag parallel parallel2 PsiUtil) 
set_property(GLOBAL APPEND PROPERTY LIBLIST frag)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(frag custom_boost)
endif()            
