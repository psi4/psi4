set(CWD ${CMAKE_CURRENT_SOURCE_DIR})
set(Srcs ${CWD}/Srcs)
set(Utils ${CWD}/Utils)
set(PsiAPIHeaders ${CWD}/PsiMethod.h
                  ${CWD}/PsiAPIBase.h
                  ${CWD}/Molecule.h
                  ${CWD}/Atom.h
)
set(PsiAPISrcs ${Srcs}/Atom.cc
               ${Srcs}/Molecule.cc
               ${Srcs}/BasisSet.cc
               ${Utils}/Exception.cc
)

set(PSIAPISRC ${PsiAPIHeaders} ${PsiAPISrcs})
add_library(psiapi STATIC ${PSIAPISRC})
add_dependencies(psiapi PsiUtil)
target_link_libraries(psiapi PsiUtil)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(molecule custom_boost)
endif()  

