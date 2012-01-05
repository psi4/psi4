#include "integraltransform.h"
#include "mospace.h"
#include <exception.h>

using namespace boost;
using namespace psi;

/**
 * Computes the DPD number that gives the most packing for a given pair of spaces
 * This version is specifically for chemists' notation integrals and is useful when
 * the spaces are not known at call time, because it can decide how to pack the
 * integrals based on whether s1 and s2 are equivalent or not.  For all other
 * purposes, use the string based routines of the same name.
 *
 * @param s1    - the first MoSpace
 * @param s2    - the second MoSpace
 * @param spin  - the spin of the first MoSpace; either Alpha or Beta can be specified
 *              - for SO spaces or for restricted transformations
 * @param pack  - if true, compute the pair number with maximum packing (for disk storage),
 *              - else assume no packing (for in-core handling)
 * @return the DPD number to use for disk storage
 */
int
IntegralTransform::DPD_ID(const shared_ptr<MOSpace> s1, const shared_ptr<MOSpace> s2,
                            SpinType spin, bool pack)
{

    std::string label = "[";
    if(s1->label() != MOSPACE_NIL && spin == Alpha){
        label += toupper(s1->label());
    }else{
        label += tolower(s1->label());
    }

    if(pack && s1->label() == s2->label()){
        label += ">=";
    }else{
        label += ",";
    }

    if(s2->label() != MOSPACE_NIL && spin == Alpha){
        label += toupper(s2->label());
    }else{
        label += tolower(s2->label());
    }

    if(pack && s1->label() == s2->label()){
        label += "]+";
    }else{
        label += "]";
    }
    if(print_>5)
        fprintf(outfile, "s1: %c s2: %c %s, label = %s, id = %d\n",
                 s1->label(), s2->label(), pack ? "packed" : "unpacked", label.c_str(), DPD_ID(label));
    return DPD_ID(label);
}

/**
 * Computes the DPD number of the space corresponding to a given MO space label
 *
 * @param c - the label of the MO space
 * @returns the number associated with the MO space in the transformation object's DPD instance
 */
int
IntegralTransform::DPD_ID(const char c)
{
    for(int i = 0; i < spacesUsed_.size(); ++i){
        if(spacesUsed_[i] == c) return(i);
    }
    std::string str("MOSpace ");
    str += c;
    str += " is not known to this transformation object";
    throw SanityCheckError(str, __FILE__, __LINE__);
}


/**
 * Computes the DPD number that gives the most packing for a given pair of spaces
 *
 * @param str - a string describing the pair packing scheme
 *
 * The following entries are possible for str
 * "[a,a]"    - Both indices belong to space 'a' and are not packed
 * "[a>a]+"   - Both indices belong to space 'a' and are packed without the diagonal elements.
 *              The tensor is symmetric w.r.t. permutation of these indices.
 * "[a>a]-"   - Both indices belong to space 'a' and are packed without the diagonal elements.
 *              The tensor is antisymmetric w.r.t. permutation of these indices.
 * "[a>=a]+"  - Both indices belong to space 'a' and are packed with the diagonal elements.
 *              The tensor is symmetric w.r.t. permutation of these indices.
 * "[a>=a]-"  - Both indices belong to space 'a' and are packed with the diagonal elements.
 *              The tensor is antisymmetric w.r.t. permutation of these indices.
 *  [a,b]     - One index belongs to space 'a', the other to 'b'.  No packing is possible.
 */
int
IntegralTransform::DPD_ID(const std::string &str)
{
    if(dpdLookup_.count(str)==0){
        std::string s = "Pair ";
        s+= str;
        s+= " has not been created.  Check the spaces passed into the IntegralTransform constructor";
        throw SanityCheckError(s, __FILE__, __LINE__);
    }
    return dpdLookup_[str];
}


/**
 * Just a wrapper to the string version, provided for those const-safe role models out there
 */
int
IntegralTransform::DPD_ID(const char *str)
{
    std::string s(str);
    return DPD_ID(s);
}


/**
 * Just a wrapper to the string version.
 */
int
IntegralTransform::DPD_ID(char *str)
{
    std::string s(str);
    return DPD_ID(s);
}
