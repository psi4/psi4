set(CWD ${CMAKE_CURRENT_SOURCE_DIR})
set(CWDSRCs ${CWD}/Srcs)
set(Impl ${CWD}/Implementations)
set(Utils ${CWD}/Utils)
set(AutoParam ${CWD}/AutoFxnalGroup)
set(MoleculeHeaders ${CWD}/Atom.h
                 ${CWD}/AtomTypes.h
                 ${CWD}/LibFragFragment.h
                 ${CWD}/FragmentedSys.h
                 ${CWD}/LibMolecule.h
                 ${CWD}/LibMoleculeBase.h
                 ${CWD}/LibFragMolecule.h
                 ${CWD}/MoleculeTypes.h
                 ${CWD}/MolItr.h
                 ${CWD}/SuperCell.h
                 ${CWD}/UnitCell.h
                 ${AutoParam}/Node.h
                 ${AutoParam}/MMAtomTypes.h
                 ${Impl}/AtomGuts.h
                 ${Impl}/FragItrGuts.h
                 ${Impl}/FragSysGuts.h
                 ${Impl}/IteratorGuts.h
                 ${Impl}/MoleculeGuts.h   
                 ${Impl}/MolItrGuts.h
                 ${Impl}/UnitCellGuts.h
                 ${Utils}/BSSEer.h
                 ${Utils}/BSSEFactory.h
                 ${Utils}/BioGeom.h
                 ${Utils}/Fragmenter.h
                 ${Utils}/GeomManipulator.h
                 ${Utils}/Geometry.h
                 ${Utils}/Orientation.h       
                 ${Utils}/OrganicGeom.h         
)
set(MoleculeSrcs ${CWDSRCs}/Atom.cc
                 ${CWDSRCs}/AtomTypes.cc
                 ${CWDSRCs}/Fragment.cc
                 ${CWDSRCs}/FragmentedSys.cc
                 ${CWDSRCs}/Molecule.cc
                 ${CWDSRCs}/MoleculeTypes.cc
                 ${CWDSRCs}/MolItr.cc
                 ${CWDSRCs}/SuperCell.cc
                 ${CWDSRCs}/UnitCell.cc
                 ${AutoParam}/Graph.cc
                 ${AutoParam}/MMAtomTypes.cc
                 ${AutoParam}/Node.cc
                 ${Impl}/AtomGuts.cc
                 ${Impl}/FragItrGuts.cc
                 ${Impl}/FragSysGuts.cc
                 ${Impl}/MoleculeGuts.cc
                 ${Impl}/MolItrGuts.cc
                 ${Impl}/UnitCellGuts.cc
                 ${Utils}/BSSEer.cc
                 ${Utils}/BSSEFactory.cc
                 ${Utils}/BioGeom.cc
                 ${Utils}/Fragmenter.cc 
                 ${Utils}/GeomManipulator.cc
                 ${Utils}/Geometry.cc
                 ${Utils}/Orientation.cc
                 ${Utils}/OrganicGeom.cc
)

set(MOLSRC ${MoleculeHeaders} ${MoleculeSrcs})
add_library(molecule STATIC ${MOLSRC})
include_directories(../libPsiUtil ../libqt ../libfrag)
set_property(GLOBAL APPEND PROPERTY LIBLIST molecule)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(molecule custom_boost)
endif()  

