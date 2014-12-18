set(CINTHeaders basisset.h CInt.h)
set(CINTSrcs basisset.c)
add_library(CINT OBJECT ${CINTHeaders} ${CINTSrcs} )