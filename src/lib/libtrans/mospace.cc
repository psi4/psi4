#include "mospace.h"
#include "libciomr/libciomr.h"

using namespace psi;
using namespace boost;

/// Keeps track of which labels have been assigned, for safety
std::map<char, int> MOSpace::labelsUsed;
shared_ptr<MOSpace> MOSpace::frc(new MOSpace(MOSPACE_FRC));
shared_ptr<MOSpace> MOSpace::occ(new MOSpace(MOSPACE_OCC));
shared_ptr<MOSpace> MOSpace::vir(new MOSpace(MOSPACE_VIR));
shared_ptr<MOSpace> MOSpace::frv(new MOSpace(MOSPACE_FRV));
shared_ptr<MOSpace> MOSpace::all(new MOSpace(MOSPACE_ALL));
shared_ptr<MOSpace> MOSpace::nil(new MOSpace(MOSPACE_NIL));

/**
 * This creates an empty MOSpace with just a label.  This is solely for the
 * construction of the pre-defined spaces; use the longer form of the constructor
 * for custom spaces.
 */
MOSpace::MOSpace(char label):
        _label(label),
        _aOrbSym(0),
        _bOrbSym(0),
        _aOrbsPI(0),
        _bOrbsPI(0),
        _aEvecs(0),
        _bEvecs(0),
        _aIndex(0),
        _bIndex(0)
{
    // Register this label as "taken"
    ++labelsUsed[_label];
}

/**
 * Defines a custom orbital space with different alpha and beta spaces
 * @param label   - a single character to label this space.  This must be unique to this
 *                  space, so see the MOSpace static member variables for a list of the
 *                  labels already used.  The uniqueness is checked internally.
 * @param nirreps - the number of irreducible representations in the Abelian point group of the molecule.
 * @param aOrbsPI - an array of dimension #irreps, describing the number of alpha orbitals per irrep
 * @param bOrbsPI - an array of dimension #irreps, describing the number of beta orbitals per irrep.
 *                  This is assumed to be the same as a aOrbsPI for restricted transformations, so
 *                  in this case it can be passed in a NULL.
 * @param aEvecs  - an array of matrices of dimension #irreps.  For each irrep h, aEvecs[h] should be an
 *                  NSO[h] X aOrbsPI[h] block matrix containing the alpha MO coefficients for this space.
 * @param bEvecs  - an array of matrices of dimension #irreps.  For each irrep h, bEvecs[h] should be an
 *                  NSO[h] X bOrbsPI[h] block matrix containing the beta MO coefficients for this space.
 *                  Pass in as NULL for restricted transformations; the alpha coefficients will be used
 *                  in this case.
 * @param aIndex  - an array of dimension #orbitals describing the number of each alpha orbital in
 *                  the space. This is only for the purposes of IWL output, so it can be passed as
 *                  NULL for DPD output.
 * @param bIndex  - an array of dimension #orbitals describing the number of each beta orbital in
 *                  the space. For restricted transformations or for DPD output only, this can be
 *                  passed as NULL.
 */
MOSpace::MOSpace(const char label, const int nirreps, const int *aOrbsPI,
                 const int *bOrbsPI, const double ***aEvecs, const double ***bEvecs,
                 const int *aIndex, const int *bIndex):
        _label(label),
        _aOrbSym(0),
        _bOrbSym(0),
        _aOrbsPI(aOrbsPI),
        _bOrbsPI(bOrbsPI),
        _aEvecs(aEvecs),
        _bEvecs(bEvecs),
        _aIndex(aIndex),
        _bIndex(bIndex)
{
    if(labelsUsed.count(label)){
        std::string error("Space ");
        error += label;
        error += " is already in use.  Choose a unique name for the custom MOSpace.";
        throw SanityCheckError(error, __FILE__, __LINE__);
    }
    ++labelsUsed[label];

    // Count the number of alpha orbitals
    int nAOrbs = 0;
    for(int h = 0; h < nirreps; ++h) nAOrbs += aOrbsPI[h];
    // Define the alpha orbital symmetries
    _aOrbSym = init_int_array(nAOrbs);
    for(int h = 0, count = 0; h < nirreps; ++h) _aOrbSym[count++] = h;
    if(bOrbsPI != NULL){
        // Count the number of beta orbitals
        int nBOrbs = 0;
        for(int h = 0; h < nirreps; ++h) nBOrbs += bOrbsPI[h];
        // Define the beta orbital symmetries
        _bOrbSym = init_int_array(nBOrbs);
        for(int h = 0, count = 0; h < nirreps; ++h) _bOrbSym[count++] = h;
    }
}


MOSpace::~MOSpace()
{
    --labelsUsed[_label];
    if(_aOrbSym != NULL) free(_aOrbSym);
    if(_bOrbSym != NULL) free(_bOrbSym);
}
