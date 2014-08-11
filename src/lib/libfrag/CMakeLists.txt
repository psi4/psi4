set(libfragsrc LibFragTypes.h
		       Fragmenter.h
               FragOptions.h
               BSSEer.h
               Embedder.h
               Capper.h
               GMBE.h
               Connections.h
               libfrag.h
               LibFragBase.h
               LibFragBase.cc
               MBE.h
               MBEFrag.h
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
               Set.cc)
add_library(frag STATIC ${libfragsrc})
               