#include "mospace.h"
#include "spaceinfo.h"

namespace psi{ namespace libtrans{

/// Keeps track of which labels have been assigned, for safety
std::map<char, int> MOSpace::labelsUsed;
shared_ptr<MOSpace> MOSpace::occ(new MOSpace(MOSPACE_OCC));
shared_ptr<MOSpace> MOSpace::vir(new MOSpace(MOSPACE_VIR));
shared_ptr<MOSpace> MOSpace::all(new MOSpace(MOSPACE_ALL));
shared_ptr<MOSpace> MOSpace::nil(new MOSpace(MOSPACE_NIL));

/**
 * This creates an empty MOSpace with just a label.  This is solely for the
 * construction of the pre-defined spaces
 */
MOSpace::MOSpace(char label):
        _label(label)
{
    // Register this label as "taken"
    ++labelsUsed[_label];
}

/**
 * Defines a custom orbital space with different alpha and beta spaces
 * @param label   - a single character to label this space.  This must be unique to this
 *                  space, so see the MOSpace static member variables for a list of the
 *                  labels already used.  The uniqueness is checked internally.
 * @param aForbPI - an array containing the first orbital in this space, per irrep for the
 *                  alpha orbitals.  These numbers are the absolute numbers in Pitzer ordering.
 * @param bForbPI - the corresponding array for the beta orbitals
 * @param aOrbsPI - an array containing the the number of alpha orbitals in this space, per irrep
 * @param bOrbsPI - the corresponding array for the beta orbitals
 */
MOSpace::MOSpace(char label, int *aFOrbPI, int *bFOrbPI, int *aOrbsPI, int *bOrbsPI):
        _label(label),
        _spaceInfo(new SpaceInfo())
{
    _spaceInfo->aFOrbPI = aFOrbPI;
    _spaceInfo->bFOrbPI = bFOrbPI;
    _spaceInfo->aOrbsPI = aOrbsPI;
    _spaceInfo->bOrbsPI = bOrbsPI;
    
    // TODO check that it doesn't exist already..
    ++labelsUsed[label];
}


MOSpace::~MOSpace()
{
    --labelsUsed[_label];
}


}} // End namespaces

