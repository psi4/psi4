set(CWD ${CMAKE_CURRENT_SOURCE_DIR})
set(CWDSRCs ${CWD}/Srcs)
set(Impl ${CWD}/Implementations)
set(Utils ${CWD}/Utils)
set(MoleculeHeaders ${CWD}/Atom.h
                 ${CWD}/AtomTypes.h
                 ${CWD}/Fragment.h
                 ${CWD}/FragmentedSys.h
                 ${CWD}/LibMolecule.h
                 ${CWD}/LibMoleculeBase.h
                 ${CWD}/Molecule.h
                 ${CWD}/MoleculeTypes.h
                 ${CWD}/MolItr.h
                 ${CWD}/SuperCell.h
                 ${CWD}/UnitCell.h
                 ${Impl}/AtomGuts.h
                 ${Impl}/FragItrGuts.h
                 ${Impl}/FragSysGuts.h
                 ${Impl}/IteratorGuts.h
                 ${Impl}/MoleculeGuts.h   
                 ${Impl}/MolItrGuts.h
                 ${Impl}/UnitCellGuts.h
                 ${Utils}/AddDerivedGrpFunctor.h
                 ${Utils}/AddGroupFunctor.h
                 ${Utils}/AromaRingFinder.h
                 ${Utils}/Fragmenter.h
                 ${Utils}/FxnalGroup.h
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
                 ${Impl}/AtomGuts.cc
                 ${Impl}/FragItrGuts.cc
                 ${Impl}/FragSysGuts.cc
                 ${Impl}/MoleculeGuts.cc
                 ${Impl}/MolItrGuts.cc
                 ${Impl}/UnitCellGuts.cc
                 ${Utils}/AromaRingFinder.cc
                 ${Utils}/Fragmenter.cc
                 ${Utils}/FxnalGroup.cc
                 ${Utils}/GeomManipulator.cc
                 ${Utils}/Geometry.cc
                 ${Utils}/Orientation.cc
                 ${Utils}/OrganicGeom.cc
)

set(MOLSRC ${MoleculeHeaders} ${MoleculeSrcs} PARENT_SCOPE)