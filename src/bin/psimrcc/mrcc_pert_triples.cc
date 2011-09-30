/**
 *  @file ccmrcc_pert_triples.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Computes the (T) correction
*/

#include <libmoinfo/libmoinfo.h>
#include <libutil/libutil.h>
#include <liboptions/liboptions.h>
#include <boost/shared_ptr.hpp>
#include <libchkpt/chkpt.hpp>

#include "mrcc.h"
#include "mrccsd_t.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;

void CCMRCC::compute_perturbative_triples()
{
  Timer timer;

  h_eff.set_eigenvalue(current_energy);
  h_eff.set_matrix(Heff,moinfo->get_nrefs());
  h_eff.set_right_eigenvector(right_eigenvector,moinfo->get_nrefs());
  h_eff.set_left_eigenvector(left_eigenvector,moinfo->get_nrefs());
  h_eff.set_zeroth_order_eigenvector(zeroth_order_eigenvector,moinfo->get_nrefs());

  MRCCSD_T mrccsd_t(options_,&h_eff);

  if(options_.get_bool("DIAGONALIZE_HEFF")){
    fprintf(outfile,"\n\n  Diagonalizing Heff");
    current_energy = h_eff.diagonalize();
  }else{
    fprintf(outfile,"\n\n  Computing the expectation value of Heff");
    current_energy = h_eff.expectation_value();
  }
  _default_chkpt_lib_->wt_etot(current_energy);

  fprintf(outfile,"\n\n%6c* Mk-MRCCSD(T) total energy   =    %20.12f",' ',current_energy);
  fprintf(outfile,"\n\n  Timing for triples:             %20.6f s",timer.get());
  fflush(outfile);
}


}}  /* End Namespaces */
