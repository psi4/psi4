set(libfragsrc Fragmenter.h
               FragOptions.h
               BSSEer.h
               GMBE.h
               libfrag.h
               LibFragPrinter.h
               MBE.h
               MBEFrag.h
               Set.h
               Fragmenter.cc
               GMBE.cc
               BSSEer.cc
               libfrag.cc
               FragOptions.cc
               MBE.cc
               Set.cc)
add_library(frag STATIC ${libfragsrc})
               