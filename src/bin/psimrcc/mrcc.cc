#include <liboptions/liboptions.h>
#include <libmoinfo/libmoinfo.h>
#include <libutil/libutil.h>

#include "algebra_interface.h"
#include "mrcc.h"
#include "matrix.h"
#include "blas.h"
#include "debugging.h"

namespace psi{ namespace psimrcc{

using namespace std;

CCMRCC::CCMRCC(Options &options):
        CCManyBody(options),
        options_(options)
{
  triples_type = ccsd;
  triples_coupling_type = cubic;
  ap_correction   = false; // Set tu true when computing the a posteriori correction
  current_energy  =  0.0;
  old_energy      = 10.0;

  // Parse the CORR_WFN parameter
  vector<string> theory_levels = split("PT2 CCSD CCSD_T CCSDT-1A CCSDT-1B CCSDT-2 CCSDT-3 CCSDT");
  for(size_t i=0;i<theory_levels.size();++i){
    if(options.get_str("CORR_WFN")==theory_levels[i])
      triples_type = TriplesType(i);
  }

  // Parse the COUPLING parameter
  vector<string> coupling_levels = split("NONE LINEAR QUADRATIC CUBIC");
  for(size_t i=0;i<coupling_levels.size();++i){
    if(options.get_str("COUPLING")==coupling_levels[i]){
      triples_coupling_type = TriplesCouplingType(i);
    }
  }

  // Parse the PERT_CBS parameter
  pert_cbs = options.get_bool("PERTURB_CBS");
  pert_cbs_coupling = options.get_bool("PERTURB_CBS_COUPLING");

  // Add the matrices that will store the intermediates
  add_matrices();

  // Generate the Fock matrices, Integrals and Denominators
  generate_integrals();
  generate_denominators();

  if(triples_type>ccsd)
    generate_triples_denominators();

  compute_reference_energy();

  DEBUGGING(1,
    blas->print_memory();
  )
}

CCMRCC::~CCMRCC()
{
}

}} /* End Namespaces */
