set(libfragsrc GMBE.h
               LibFragDriver.h
			   LibFragDriver.cc
			   MBE.h
			   MBE.cc
			   MBEProp.h
)
add_library(frag STATIC ${libfragsrc})  
include_directories(../libPsiUtil ../libqt ../libmolecule)
add_dependencies(frag parallel2 PsiUtil psiapi)
target_link_libraries(frag parallel parallel2 PsiUtil molecule psiapi) 
set_property(GLOBAL APPEND PROPERTY LIBLIST frag)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(frag custom_boost)
endif()            
