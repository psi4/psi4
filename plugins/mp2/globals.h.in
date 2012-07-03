#include <liboptions/liboptions.h>
#include <libtrans/integraltransform.h>
#include <libchkpt/chkpt.hpp>

#ifdef EXTERN
    #undef EXTERN
    #define EXTERN extern
#else
    #define EXTERN
#endif

//#define INDEX(i,j) ((i>j) ? ((i*(i+1)/2)+j) : ((j*(j+1)/2)+i))
#define ID(x) ints.DPD_ID(x)

using namespace boost;

namespace psi{ namespace plugin_mp2{
    // Nasty, nasty global variables.
    EXTERN int nirreps, nmo;
    EXTERN shared_ptr<PSIO> psio;
    EXTERN int *aOccOrbsPI, *bOccOrbsPI, *aVirOrbsPI, *bVirOrbsPI;
    EXTERN int *mopi, *clsdpi, *openpi, *frzcpi, *frzvpi;
    EXTERN int numAOcc, numBOcc, numAVir, numBVir;
    EXTERN double eSCF;
    EXTERN double *aOccEvals, *bOccEvals, *aVirEvals, *bVirEvals;
    
    EXTERN double plugin_mp2_unrestricted(Options &options, shared_ptr<Chkpt> chkpt);
    EXTERN double plugin_mp2_restricted(Options &options, shared_ptr<Chkpt> chkpt);
}} // Namespaces
