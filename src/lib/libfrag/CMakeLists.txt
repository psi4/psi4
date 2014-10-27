set(libfragsrc Fragmenter.h
               FragOptions.h
               BSSEer.h
               Embedder.h
               Capper.h
               GMBE.h
               Connections.h
               libfrag.h
               MBE.h
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
)
add_library(frag STATIC ${libfragsrc})  
target_link_libraries(frag parallel)             
if(BUILD_CUSTOM_BOOST)
   add_dependencies(frag custom_boost)
endif()   
