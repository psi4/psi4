include_directories(../GTFock/include)
include_directories(../GTFock/libcint)
include_directories(../GTFock/libpfock)
set(SrcSource JKFactory.hpp
	          JKFactory.cpp  
               GTFockInter.hpp  
               GTFockInter.cpp  
               MPIManager.hpp
               Molecule.hpp
               Molecule.cpp
               BasisSet.hpp  
               BasisSet.cpp
               BasisReorderer.hpp
               BasisReorderer.cpp  
               Interface.hpp
               Interface.cpp
               Matrix.hpp
               Matrix.cpp)
add_library(SRCOBJS OBJECT ${SrcSource})
