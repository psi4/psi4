#include <libplugin/plugin.h>
#include"psi4-dec.h"
#include<libdpd/dpd.h>

#include<boost/shared_ptr.hpp>
#include<liboptions/liboptions.h>
#include<libtrans/integraltransform.h>
#include<libtrans/mospace.h>
#include<libmints/matrix.h>
#include<libmints/vector.h>
#include<libchkpt/chkpt.h>
#include<libiwl/iwl.h>
#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
#include <ccfiles.h>
#include <cstdio>
#include <cstdlib>


#include"globals.h"
#include"gpuhelper.h"
#include"blas.h"
#include"ccsd.h"

using namespace psi;
using namespace boost;

// Forward declaration to call cctriples
namespace psi { namespace cctriples {
PsiReturnType cctriples(Options &options);
}}
// Forward declaration to call ccsort
namespace psi { namespace ccsort {
PsiReturnType ccsort(Options &options);
}}
// Forward declaration to call ccenergy
namespace psi { namespace ccenergy {
  PsiReturnType ccenergy(Options &options);
  struct dpd_file4_cache_entry *priority_list(void);
  int **cacheprep_rhf(int level, int *cachefiles);
}}

// wrapper to crawford's triples code
namespace psi{
  PsiReturnType triples(boost::shared_ptr<psi::CoupledCluster>ccsd,Options&options);
}



namespace psi{

void RunCoupledCluster(Options &options){

  boost::shared_ptr<CoupledCluster> ccsd(new CoupledCluster);
  PsiReturnType status;

  tstart();
  ccsd->WriteBanner();
  ccsd->Initialize(options);
  ccsd->AllocateMemory(options);
  status = ccsd->CCSDIterations();
  tstop();

  if (options.get_str("WFN")=="CCSD_T"){
     free(ccsd->tempt);
     free(ccsd->tempv);
     free(ccsd->integrals);
     status = psi::triples(ccsd,options);
     if (status == Failure){
        throw PsiException( 
           "Whoops, the (T) correction died.",__FILE__,__LINE__);
     }
  }

  ccsd.reset();

}

}


