set(libfragsrc LibFragTypes.h
		       Fragmenter.h
               FragOptions.h
               BSSEer.h
               Capper.h
               GMBE.h
               Connections.h
               libfrag.h
               LibFragPrinter.h
               MBE.h
               MBEFrag.h
               Set.h
               AtomSet.h
               Fragmenter.cc
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
               