/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cmath>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<memory.h>
#include<pthread.h>
#include<libipv1/ip_lib.h>
#include<libciomr/libciomr.h>
#include<libpsio/psio.h>
#include<libint/libint.h>
#include<libqt/qt.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

#include"read_scf_opdm.h"
#include"read_scf_occ_evec.h"
#include"shell_block_matrix.h"
#include"hf_fock.h"
#include"xc_fock.h"
#include"xc_fock_u.h"
#include"xc_grad_fock.h"

pthread_mutex_t fock_mutex;            /* Lock on the global AO matrix */

namespace psi {
  namespace CINTS {
    
    void fock()
    {
      int dum;
      int n, num;
      int i, j, k, l;

      int nstri;
      double temp;
      double **tmpmat1;
      double *Gtri, *Gtri_o;  /* Total and open-shell G matrices 
				 and lower triagonal form in SO basis */
  
      /*----------------------------------------
	Read in the difference HF/DFT densities
	----------------------------------------*/
      timer_on("HF_FOCK");
      read_scf_opdm();
      
      /*-------------------------------------------
	Compute HF contribution to the Fock matrix
	-------------------------------------------*/
      
      hf_fock();
      timer_off("HF_FOCK");
      /*-----------------------------------
	Do numerical interation for KS DFT
	-----------------------------------*/
      if(UserOptions.make_dft){
	/*--- Read in the SCF eigenvector density ---*/
	read_scf_occ_evec();
	/*-- Compute exch+corr contribution to the Fock matrix ---*/
	if(UserOptions.reftype == rhf)
	  xc_grad_fock();
    else if(UserOptions.reftype == uhf)
      xc_fock_u();
    else
	throw std::domain_error("\nUnrecognized Kohn-Sham DFT reference");
  }

  return;
}
}}
