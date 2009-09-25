/*! 
** \file
** \ingroup MP2R12
** \brief MOInfo structure header file
*/

/* Struct for PSIF_CHKPT molecular orbital information */

#ifndef _psi_src_bin_mp2r12_moinfo_h_
#define _psi_src_bin_mp2r12_moinfo_h_

namespace psi{
  namespace mp2r12{

  struct MOInfo {
    int nmo;
    int nirreps;
    int iopen;
    int *orbspi;
    int *clsdpi;
    int *openpi;
    int *virtpi;
    int *frdocc;
    int *fruocc;
    int *first;
    int *last;
    int noeints;
    int nteints;
    char **labels;
    double enuc;
    double escf;
    double *evals;
    double *te_ints;
  };

  } /* Namespace mp2r12 */
} /* Namespace psi */
 
# endif /* Header guard */   
