set(libfragsrc Fragmenter.h
               FragOptions.h
               GMBE.h
               InputManip.h
               libfrag.h
               LibFragPrinter.h
               MBE.h
               MBEFrag.h
               Set.h
               Fragmenter.cc
               GMBE.cc
               InputManip.cc
               libfrag.cc
               MBE.cc
               Set.cc)
add_library(frag STATIC ${libfragsrc})
               