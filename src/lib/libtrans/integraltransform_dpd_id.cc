#include "integraltransform.h"
#include "mospace.h"

namespace psi{ namespace libtrans{
/**
 * Computes the DPD number that gives the most packing for a given pair of spaces
 * @param s1    - the first MoSpace
 * @param spin1 - the spin of the first MoSpace; either Alpha or Beta can be specified
 *              - for SO spaces or for restricted transformations
 * @param s2    - the second MoSpace
 * @param spin2 - the spin of the second MoSpace; either Alpha or Beta can be specified
 *              - for SO spaces or for restricted transformations
 * @param pack  - if true, compute the pair number with maximum packing (for disk storage),
 *              - else assume no packing (for in-core handling)
 * @return  - the DPD number to use for disk storage
 */
int
IntegralTransform::DPD_ID(shared_ptr<MOSpace> s1, shared_ptr<MOSpace> s2,
                            SpinType spin, bool pack)
{
    int id;
    if(_aSpaceNum.count(s1->label()) == 0){
        std::string str = "This IntegralTransform object was not initialized for MO space " + s1->label();
        throw SanityCheckError(str, __FILE__, __LINE__);
    }
    if(_aSpaceNum.count(s2->label()) == 0){
        std::string str = "This IntegralTransform object was not initialized for MO space " + s2->label();
        throw SanityCheckError(str, __FILE__, __LINE__);
    }

    // These are the positions of the spaces in the list
    int pos1 = spin == Alpha ? _aSpaceNum[s1->label()] : _bSpaceNum[s1->label()];
    int pos2 = spin == Alpha ? _aSpaceNum[s2->label()] : _bSpaceNum[s2->label()];

    // Now we know where these things are in the list of spaces, go ahead and compute the pair number
    int numSpaces = _spacesUsed.size();
    if(s1 == s2){
        if(pack){
            // Pack the indices for p>=q
            id = 5 * pos1 + 3;
        }else{
            // No packing, but the same type of space: p,q
            id = 5 * pos1 + 0;
        }
    }else{
        // Offset by the number of same-space pairings created
        id = 5 * numSpaces;
        /* This looks a bit weird, but it comes from the closed form of an
         * arithmetic progression; the "mixed" pair indices are ordered
         * 0 1
         * 0 2
         * 0 3
         * . .
         * . .
         * 0 n-1
         * 1 2
         * 1 3
         * . .
         * . .
         */
        if(pos1 < pos2){
            id += pos1 * (2*(numSpaces-1)-(pos1-1)) + 2*(pos2 - pos1 - 1);
        }else{
            id += pos2 * (2*(numSpaces-1)-(pos2-1)) + 2*(pos1 - pos2 - 1) + 1;
        }
    }
    return id;
}

}} // Namespaces