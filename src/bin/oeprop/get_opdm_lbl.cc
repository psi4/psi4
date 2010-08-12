/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "includes.h"
#include "globals.h"
#include "prototypes.h"
#include <libqt/qt.h>
#include <psifiles.h>
#include <ccfiles.h>

namespace psi { namespace oeprop {

void get_opdm_lbl(Options & options) { 
  int i, errcod;
  int root;

  nrho = 1;
  if (options["NUM_ROOTS"].has_changed()) {
    nrho = options.get_int("NUM_ROOTS");
  }
  else if (cc_wfn(wfn)) {
    psio_open(CC_INFO,1);
    psio_read_entry(CC_INFO, "Num. of CC densities", (char *) &(nrho), 
      sizeof(int));
    psio_close(CC_INFO,1);
  }
  else if (ci_wfn(wfn)) {
    psio_open(opdm_file,1);
    if (transdens) 
      psio_read_entry(opdm_file, "Num MO-basis TDM", (char *) &(nrho),
        sizeof(int));
    else
      psio_read_entry(opdm_file, "Num MO-basis OPDM", (char *) &(nrho),
        sizeof(int));
    psio_close(opdm_file,1);
  }

  if (nrho < 1) {
    fprintf(outfile,"(oeprop): error - got nrho = %d\n", nrho);
    throw PsiException("error - got nrho", __FILE__, __LINE__);
  }

  /* if the ROOT keyword is specified, let's just analyze the ROOT given 
     and not however many there may be on disk */
  if (options["ROOT"].has_changed()) {
    root = options.get_int("ROOT");
    root -= 1;
    nrho = 1;
    opdm_lbl = (char **) malloc(sizeof(char *) * nrho);
    opdm_lbl[0] = (char *) malloc(32*sizeof(char));

    if (transdens)
      sprintf(opdm_lbl[0],"MO-basis TDM Root %d",root);
    else
      sprintf(opdm_lbl[0],"MO-basis OPDM Root %d",root);
  } /* end ROOT exists */

  /* if ROOT is not given and only one density specified, 
     then let's analyze "THE" density */
  else if (nrho == 1) {
    if ( ref == "RHF" || ref == "ROHF" ) {
      opdm_lbl = (char **) malloc(sizeof(char *) * nrho);
      opdm_lbl[0] = (char *) malloc(32*sizeof(char));
      sprintf(opdm_lbl[0], "MO-basis %s", transdens ? "TDM" : "OPDM");
    }
    else {
      opdm_a_lbl = (char **) malloc(sizeof(char *) * nrho);
      opdm_a_lbl[0] = (char *) malloc(32*sizeof(char));
      opdm_b_lbl = (char **) malloc(sizeof(char *) * nrho);
      opdm_b_lbl[0] = (char *) malloc(32*sizeof(char));
      sprintf(opdm_a_lbl[0], "MO-basis Alpha %s", transdens ? "TDM" : "OPDM");
      sprintf(opdm_b_lbl[0], "MO-basis Beta %s", transdens ? "TDM" : "OPDM");
    } 
  } /* end "THE" density */

  if (transdens) 
    fprintf(outfile, "\tTransition densities available up to root %d\n", nrho);
  else
    fprintf(outfile, "\tDensities available up to root %d\n", nrho);

  if (nrho == 1) return;

  if ( ref == "RHF" || ref == "ROHF" ) {
    opdm_lbl = (char **) malloc(sizeof(char *) * nrho);
    for (i=0;i<nrho;i++) {
      opdm_lbl[i] = (char *) malloc(32*sizeof(char));
      if (i==0 && cc_wfn(wfn)) sprintf(opdm_lbl[i], "MO-basis %s",
        transdens ? "TDM" : "OPDM");
      else sprintf(opdm_lbl[i], "MO-basis %s Root %d", 
        transdens ? "TDM" : "OPDM", i);
    }
  }
  else {
    opdm_a_lbl = (char **) malloc(sizeof(char *) * nrho);
    opdm_b_lbl = (char **) malloc(sizeof(char *) * nrho);
    for (i=0;i<nrho;i++) {
      opdm_a_lbl[i] = (char *) malloc(32*sizeof(char));
      opdm_b_lbl[i] = (char *) malloc(32*sizeof(char));
      if (i==0 && cc_wfn(wfn)) {
        sprintf(opdm_a_lbl[i], "MO-basis Alpha %s", 
          transdens ? "TDM" : "OPDM");
        sprintf(opdm_b_lbl[i], "MO-basis Beta %s",
          transdens ? "TDM" : "OPDM");
      }
      else {
        sprintf(opdm_a_lbl[i], "MO-basis Alpha %s Root %d", 
          transdens ? "TDM" : "OPDM", i);
        sprintf(opdm_b_lbl[i], "MO-basis Beta %s Root %d", 
          transdens ? "TDM" : "OPDM", i);
      }
    }
  }

  return;

}

}} // namespace psi::oeprop
