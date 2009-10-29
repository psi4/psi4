#ifndef _PSI_SRC_LIB_LIBTRANS_SPACEINFO_H_
#define _PSI_SRC_LIB_LIBTRANS_SPACEINFO_H_

using namespace psi;

namespace psi{ namespace libtrans{

/**
* SpaceInfo is just a structure that contains all of the essential information
* about an MOSpace needed for the transformation.  It provides a convenient way to give
* default values to the components and allows multiple copies of the information for a
* given space.  There might be times where mulitple instances of the occupied space
* are required, for example, where one freezes the core and the other does not.
* For this purpose, it's better to include the SpaceInfo (which can be replicated)
* as a member of MOSpace (which cannot be replicated, because of the unique labelling).
*/
class SpaceInfo {
  public:
    SpaceInfo():
        aOrbsPI(NULL),
        bOrbsPI(NULL),
        aFOrbPI(NULL),
        bFOrbPI(NULL),
        Ca(NULL),
        Cb(NULL)
    {}

    ~SpaceInfo(){
    }

    int *aOrbsPI;
    int *bOrbsPI;
    int *aFOrbPI;
    int *bFOrbPI;
    double ***Ca;
    double ***Cb;
};


}} // End namespaces

// This is here so that files including this have clean(er) syntax
using namespace psi::libtrans;

#endif // Header guard

