set(libfragsrc 
#               libmolecule/Atom.h
#               libmolecule/AtomTypes.h
#               libmolecule/LibMoleculeBase.h
#               libmolecule/Fragment.h
#               libmolecule/Molecule.h            
#               libmolecule/Atom.cc
#               libmolecule/AtomTypes.cc
#               libmolecule/Fragment.cc
#               libmolecule/Molecule.cc
               Fragmenter.h
               FragOptions.h
               BSSEer.h
               Embedder.h
               Capper.h
               GMBE.h
               Connections.h
               libfrag.h
               MBE.h
               Grouper.h
               Grouper.cc
               MBEFrag.h
               MBEFrag.cc
               MBEFragSet.h
               MBEFragSet.cc
               Set.h
               AtomSet.h
               Fragmenter.cc
               Embedder.cc
               GMBE.cc
               Connections.cc
               BSSEer.cc
               libfrag.cc
               FragOptions.cc
               MBE.cc
               Capper.cc
               AtomSet.cc
               Symmetrizer.h
               Symmetrizer.cc
               DirectSum.cc
)
add_library(frag STATIC ${libfragsrc})  
include_directories(../libPsiUtil)
target_link_libraries(frag parallel PsiUtil) 
if(BUILD_CUSTOM_BOOST)
   add_dependencies(frag custom_boost)
endif()            
