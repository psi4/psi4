/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libiwl/iwl.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

/*
** This is the ugliest piece of code I have written in my life.
** T. Daniel Crawford, September, 1996
*/

void classify(int p, int q, int r, int s, double value,
	      struct iwlbuf *ABuf, struct iwlbuf *BBuf,
	      struct iwlbuf *CBuf, struct iwlbuf *DBuf,
	      struct iwlbuf *EBuf, struct iwlbuf *F1Buf,
	      struct iwlbuf *F2Buf)
{
  int *occ, *vir, *socc;
  int *cc_occ, *cc_vir;
  int dirac=1;
  int soccs;
  int nfzc, *frozen;

  nfzc = moinfo.nfzc;
  frozen = moinfo.frozen;

  if(params.dertype==1) { /* Skip the frozen orbitals in the list */
    if(frozen[p] || frozen[q] || frozen[r] || frozen[s]) return;
    else { /* Adjust the indices to account for the frozen core */
      p -= nfzc; q -= nfzc;
      r -= nfzc; s -= nfzc;
    }
  }

  occ = moinfo.occ; 
  vir = moinfo.vir; 
  socc = moinfo.socc;
  cc_occ = moinfo.cc_occ; 
  cc_vir = moinfo.cc_vir;

  soccs = socc[p] + socc[q] + socc[r] + socc[s];
 
  /* A (oo|oo) integrals */
  if((occ[p] && occ[q] && occ[r] && occ[s]))
    iwl_buf_wrt_val(ABuf, cc_occ[p], cc_occ[q], cc_occ[r], cc_occ[s],
		    value, 0, outfile, dirac);

  /* B (vv|vv) integrals */
  if(params.make_abcd) {
    if((vir[p] && vir[q] && vir[r] && vir[s]))
      iwl_buf_wrt_val(BBuf, cc_vir[p], cc_vir[q], cc_vir[r], cc_vir[s],
		      value, 0, outfile, dirac);
  }

  /* C (oo|vv) integrals */
  if(soccs > 1) {
    if((occ[p] && occ[q] && vir[r] && vir[s]))
      iwl_buf_wrt_val(CBuf, cc_occ[p], cc_occ[q], cc_vir[r], cc_vir[s],
		      value, 0, outfile, dirac);
    if((occ[r] && occ[s] && vir[p] && vir[q]))
      iwl_buf_wrt_val(CBuf, cc_occ[r], cc_occ[s], cc_vir[p], cc_vir[q],
		      value, 0, outfile, dirac);
  }
  else if((occ[p] && occ[q] && vir[r] && vir[s]))
    iwl_buf_wrt_val(CBuf, cc_occ[p], cc_occ[q], cc_vir[r], cc_vir[s],
		    value, 0, outfile, dirac);
  else if((occ[r] && occ[s] && vir[p] && vir[q]))
    iwl_buf_wrt_val(CBuf, cc_occ[r], cc_occ[s], cc_vir[p], cc_vir[q],
		    value, 0, outfile, dirac);

  /* D (ov|ov) integrals */
  if(soccs > 1) {
    if((occ[p] && vir[q] && occ[r] && vir[s]))
      iwl_buf_wrt_val(DBuf, cc_occ[p], cc_vir[q], cc_occ[r], cc_vir[s],
		      value, 0, outfile, dirac);
    if((occ[q] && vir[p] && occ[r] && vir[s]))
      iwl_buf_wrt_val(DBuf, cc_occ[q], cc_vir[p], cc_occ[r], cc_vir[s],
		      value, 0, outfile, dirac);
    if((occ[p] && vir[q] && occ[s] && vir[r]))
      iwl_buf_wrt_val(DBuf, cc_occ[p], cc_vir[q], cc_occ[s], cc_vir[r],
		      value, 0, outfile, dirac);
    if((occ[q] && vir[p] && occ[s] && vir[r]))
      iwl_buf_wrt_val(DBuf, cc_occ[q], cc_vir[p], cc_occ[s], cc_vir[r],
		      value, 0, outfile, dirac);
  }
  else if((occ[p] && vir[q] && occ[r] && vir[s]))
    iwl_buf_wrt_val(DBuf, cc_occ[p], cc_vir[q], cc_occ[r], cc_vir[s],
		    value, 0, outfile, dirac);
  else if((occ[q] && vir[p] && occ[r] && vir[s]))
    iwl_buf_wrt_val(DBuf, cc_occ[q], cc_vir[p], cc_occ[r], cc_vir[s],
		    value, 0, outfile, dirac);
  else if((occ[p] && vir[q] && occ[s] && vir[r]))
    iwl_buf_wrt_val(DBuf, cc_occ[p], cc_vir[q], cc_occ[s], cc_vir[r],
		    value, 0, outfile, dirac);
  else if((occ[q] && vir[p] && occ[s] && vir[r]))
    iwl_buf_wrt_val(DBuf, cc_occ[q], cc_vir[p], cc_occ[s], cc_vir[r],
		    value, 0, outfile, dirac);

  /* E (vo|oo) integrals */
  if(soccs > 1) {
    if((vir[p] && occ[q] && occ[r] && occ[s]))
      iwl_buf_wrt_val(EBuf, cc_vir[p], cc_occ[q], cc_occ[r], cc_occ[s],
		      value, 0, outfile, dirac);
    if((vir[q] && occ[p] && occ[r] && occ[s]))
      iwl_buf_wrt_val(EBuf, cc_vir[q], cc_occ[p], cc_occ[r], cc_occ[s],
		      value, 0, outfile, dirac);
    if((vir[r] && occ[s] && occ[p] && occ[q]))
      iwl_buf_wrt_val(EBuf, cc_vir[r], cc_occ[s], cc_occ[p], cc_occ[q],
		      value, 0, outfile, dirac);
    if((vir[s] && occ[r] && occ[p] && occ[q]))
      iwl_buf_wrt_val(EBuf, cc_vir[s], cc_occ[r], cc_occ[p], cc_occ[q],
		      value, 0, outfile, dirac);
  } 
  else if((vir[p] && occ[q] && occ[r] && occ[s]))
    iwl_buf_wrt_val(EBuf, cc_vir[p], cc_occ[q], cc_occ[r], cc_occ[s],
		    value, 0, outfile, dirac);
  else if((vir[q] && occ[p] && occ[r] && occ[s]))
    iwl_buf_wrt_val(EBuf, cc_vir[q], cc_occ[p], cc_occ[r], cc_occ[s],
		    value, 0, outfile, dirac);
  else if((vir[r] && occ[s] && occ[p] && occ[q]))
    iwl_buf_wrt_val(EBuf, cc_vir[r], cc_occ[s], cc_occ[p], cc_occ[q],
		    value, 0, outfile, dirac);
  else if((vir[s] && occ[r] && occ[p] && occ[q]))
    iwl_buf_wrt_val(EBuf, cc_vir[s], cc_occ[r], cc_occ[p], cc_occ[q], 
		    value, 0, outfile, dirac);

  /* F (ov|vv) integrals */
  if(soccs > 1) {
    if((occ[p] && vir[q] && vir[r] && vir[s])) {
      iwl_buf_wrt_val(F1Buf, cc_occ[p], cc_vir[q], cc_vir[r], cc_vir[s],
		      value, 0, outfile, dirac);
      if(params.make_aibc)
	iwl_buf_wrt_val(F2Buf, cc_vir[r], cc_vir[s], cc_occ[p], cc_vir[q],
		      value, 0, outfile, dirac);
    }
    if((occ[q] && vir[p] && vir[r] && vir[s])) {
      iwl_buf_wrt_val(F1Buf, cc_occ[q], cc_vir[p], cc_vir[r], cc_vir[s],
		      value, 0, outfile, dirac);
      if(params.make_aibc)
	iwl_buf_wrt_val(F2Buf, cc_vir[r], cc_vir[s], cc_occ[q], cc_vir[p],
		      value, 0, outfile, dirac);
    }
    if((occ[r] && vir[s] && vir[p] && vir[q])) {
      iwl_buf_wrt_val(F1Buf, cc_occ[r], cc_vir[s], cc_vir[p], cc_vir[q],
		      value, 0, outfile, dirac);
      if(params.make_aibc)
	iwl_buf_wrt_val(F2Buf, cc_vir[p], cc_vir[q], cc_occ[r], cc_vir[s],
		      value, 0, outfile, dirac);
    }
    if((occ[s] && vir[r] && vir[p] && vir[q])) {
      iwl_buf_wrt_val(F1Buf, cc_occ[s], cc_vir[r], cc_vir[p], cc_vir[q],
		      value, 0, outfile, dirac);
      if(params.make_aibc)
	iwl_buf_wrt_val(F2Buf, cc_vir[p], cc_vir[q], cc_occ[s], cc_vir[r],
		      value, 0, outfile, dirac);
    }
  }
  else if((occ[p] && vir[q] && vir[r] && vir[s])) {
    iwl_buf_wrt_val(F1Buf, cc_occ[p], cc_vir[q], cc_vir[r], cc_vir[s],
		    value, 0, outfile, dirac);
    if(params.make_aibc)
      iwl_buf_wrt_val(F2Buf, cc_vir[r], cc_vir[s], cc_occ[p], cc_vir[q],
		    value, 0, outfile, dirac);
  }
  else if((occ[q] && vir[p] && vir[r] && vir[s])) {
    iwl_buf_wrt_val(F1Buf, cc_occ[q], cc_vir[p], cc_vir[r], cc_vir[s],
		    value, 0, outfile, dirac);
    if(params.make_aibc)
      iwl_buf_wrt_val(F2Buf, cc_vir[r], cc_vir[s], cc_occ[q], cc_vir[p],
		    value, 0, outfile, dirac);
  }
  else if((occ[r] && vir[s] && vir[p] && vir[q])) {
    iwl_buf_wrt_val(F1Buf, cc_occ[r], cc_vir[s], cc_vir[p], cc_vir[q],
		    value, 0, outfile, dirac);
    if(params.make_aibc)
      iwl_buf_wrt_val(F2Buf, cc_vir[p], cc_vir[q], cc_occ[r], cc_vir[s],
		    value, 0, outfile, dirac);
  }
  else if((occ[s] && vir[r] && vir[p] && vir[q])) {
    iwl_buf_wrt_val(F1Buf, cc_occ[s], cc_vir[r], cc_vir[p], cc_vir[q],
		    value, 0, outfile, dirac);
    if(params.make_aibc)
      iwl_buf_wrt_val(F2Buf, cc_vir[p], cc_vir[q], cc_occ[s], cc_vir[r],
		    value, 0, outfile, dirac);
  }
}


void classify_uhf(int p, int q, int r, int s, double value, const char *spin,
	          struct iwlbuf *ABuf1, struct iwlbuf *BBuf1,
		  struct iwlbuf *CBuf1, struct iwlbuf *CBuf2, 
		  struct iwlbuf *DBuf1, struct iwlbuf *EBuf1,
		  struct iwlbuf *EBuf2, struct iwlbuf *FBuf1, 
		  struct iwlbuf *FBuf2, struct iwlbuf *FBuf3,
		  struct iwlbuf *FBuf4)
{
  int *occ1, *occ2, *vir1, *vir2;
  int *cc_occ1, *cc_occ2, *cc_vir1, *cc_vir2;
  int dirac=1;
  int *frozen, nfzc;

  nfzc = moinfo.nfzc;
  frozen = moinfo.frozen;

  if(params.dertype==1) { /* Skip the frozen orbitals in the list */
    if(frozen[p] || frozen[q] || frozen[r] || frozen[s]) return;
    else { /* Adjust the indices to account for the frozen core */
      p -= nfzc; q -= nfzc;
      r -= nfzc; s -= nfzc;
    }
  }

  if(!strcmp(spin,"AA")) {
    occ1 = moinfo.aocc; 
    occ2 = moinfo.aocc;
    vir1 = moinfo.avir;
    vir2 = moinfo.avir;
    cc_occ1 = moinfo.cc_aocc;
    cc_occ2 = moinfo.cc_aocc;
    cc_vir1 = moinfo.cc_avir;
    cc_vir2 = moinfo.cc_avir;
  }
  else if(!strcmp(spin,"BB")) {
    occ1 = moinfo.bocc; 
    occ2 = moinfo.bocc;
    vir1 = moinfo.bvir;
    vir2 = moinfo.bvir;
    cc_occ1 = moinfo.cc_bocc;
    cc_occ2 = moinfo.cc_bocc;
    cc_vir1 = moinfo.cc_bvir;
    cc_vir2 = moinfo.cc_bvir;
  }
  else if(!strcmp(spin,"AB")) {
    occ1 = moinfo.aocc; 
    occ2 = moinfo.bocc;
    vir1 = moinfo.avir;
    vir2 = moinfo.bvir;
    cc_occ1 = moinfo.cc_aocc;
    cc_occ2 = moinfo.cc_bocc;
    cc_vir1 = moinfo.cc_avir;
    cc_vir2 = moinfo.cc_bvir;
  }

  /* B (vv|vv) integrals */
  if(params.make_abcd) {
    if(vir1[p] && vir1[q] && vir2[r] && vir2[s])
      iwl_buf_wrt_val(BBuf1, cc_vir1[p], cc_vir1[q], cc_vir2[r], cc_vir2[s],
		      value, 0, outfile, dirac);
  }

  /* F (ov|vv) and (vv|ov) integrals */
  if(occ1[p] && vir1[q] && vir2[r] && vir2[s]) {
    iwl_buf_wrt_val(FBuf1, cc_occ1[p], cc_vir1[q], cc_vir2[r], cc_vir2[s],
		    value, 0, outfile, dirac);
    iwl_buf_wrt_val(FBuf2, cc_vir2[r], cc_vir2[s], cc_occ1[p], cc_vir1[q],
		    value, 0, outfile, dirac);
  }
  else if(occ1[q] && vir1[p] && vir2[r] && vir2[s]) {
    iwl_buf_wrt_val(FBuf1, cc_occ1[q], cc_vir1[p], cc_vir2[r], cc_vir2[s],
		    value, 0, outfile, dirac);
    iwl_buf_wrt_val(FBuf2, cc_vir2[r], cc_vir2[s], cc_occ1[q], cc_vir1[p],
		    value, 0, outfile, dirac);
  }
  else if(occ1[r] && vir1[s] && vir2[p] && vir2[q] && strcmp(spin,"AB")) {
    iwl_buf_wrt_val(FBuf1, cc_occ1[r], cc_vir1[s], cc_vir2[p], cc_vir2[q],
		    value, 0, outfile, dirac);
    iwl_buf_wrt_val(FBuf2, cc_vir2[p], cc_vir2[q], cc_occ1[r], cc_vir1[s],
		    value, 0, outfile, dirac);
  }
  else if(occ1[s] && vir1[r] && vir2[p] && vir2[q] && strcmp(spin,"AB")) {
    iwl_buf_wrt_val(FBuf1, cc_occ1[s], cc_vir1[r], cc_vir2[p], cc_vir2[q],
		    value, 0, outfile, dirac);
    iwl_buf_wrt_val(FBuf2, cc_vir2[p], cc_vir2[q], cc_occ1[s], cc_vir1[r],
		    value, 0, outfile, dirac);
  }


  if(!strcmp(spin,"AB")) {
    /* F (vv|vo) and (vv|ov) integrals */
    if(vir1[p] && vir1[q] && vir2[r] && occ2[s]) {
      iwl_buf_wrt_val(FBuf3, cc_vir1[p], cc_vir1[q], cc_vir2[r], cc_occ2[s],
		      value, 0, outfile, dirac);
      iwl_buf_wrt_val(FBuf4, cc_vir1[p], cc_vir1[q], cc_occ2[s], cc_vir2[r],
		      value, 0, outfile, dirac);
    }
    else if(vir1[p] && vir1[q] && vir2[s] && occ2[r]) {
      iwl_buf_wrt_val(FBuf3, cc_vir1[p], cc_vir1[q], cc_vir2[s], cc_occ2[r],
		      value, 0, outfile, dirac); 
      iwl_buf_wrt_val(FBuf4, cc_vir1[p], cc_vir1[q], cc_occ2[r], cc_vir2[s],
		      value, 0, outfile, dirac); 
   }
  }

  /* C (oo|vv) integrals */
  if(occ1[p] && occ1[q] && vir2[r] && vir2[s])
    iwl_buf_wrt_val(CBuf1, cc_occ1[p], cc_occ1[q], cc_vir2[r], cc_vir2[s],
		    value, 0, outfile, dirac);
  else if(occ1[r] && occ1[s] && vir2[p] && vir2[q] && strcmp(spin,"AB"))
    iwl_buf_wrt_val(CBuf1, cc_occ1[r], cc_occ1[s], cc_vir2[p], cc_vir2[q],
		    value, 0, outfile, dirac);

  if(!strcmp(spin,"AB")) {
    /* C (vv|oo) integrals */
    if(vir1[p] && vir1[q] && occ2[r] && occ2[s])
      iwl_buf_wrt_val(CBuf2, cc_vir1[p], cc_vir1[q], cc_occ2[r], cc_occ2[s],
		      value, 0, outfile, dirac);
  }

  /* D (ov|ov) integrals */
  if(occ1[p] && vir1[q] && occ2[r] && vir2[s])
    iwl_buf_wrt_val(DBuf1, cc_occ1[p], cc_vir1[q], cc_occ2[r], cc_vir2[s],
		    value, 0, outfile, dirac);
  else if(occ1[q] && vir1[p] && occ2[r] && vir2[s])
    iwl_buf_wrt_val(DBuf1, cc_occ1[q], cc_vir1[p], cc_occ2[r], cc_vir2[s],
		    value, 0, outfile, dirac);
  else if(occ1[p] && vir1[q] && occ2[s] && vir2[r])
    iwl_buf_wrt_val(DBuf1, cc_occ1[p], cc_vir1[q], cc_occ2[s], cc_vir2[r],
		    value, 0, outfile, dirac);
  else if(occ1[q] && vir1[p] && occ2[s] && vir2[r])
    iwl_buf_wrt_val(DBuf1, cc_occ1[q], cc_vir1[p], cc_occ2[s], cc_vir2[r],
		    value, 0, outfile, dirac);

  /* E (vo|oo) integrals */
  if(vir1[p] && occ1[q] && occ2[r] && occ2[s])
    iwl_buf_wrt_val(EBuf1, cc_vir1[p], cc_occ1[q], cc_occ2[r], cc_occ2[s],
		    value, 0, outfile, dirac);
  else if(vir1[q] && occ1[p] && occ2[r] && occ2[s])
    iwl_buf_wrt_val(EBuf1, cc_vir1[q], cc_occ1[p], cc_occ2[r], cc_occ2[s],
		    value, 0, outfile, dirac);
  else if(vir1[r] && occ1[s] && occ2[p] && occ2[q] && strcmp(spin,"AB"))
    iwl_buf_wrt_val(EBuf1, cc_vir1[r], cc_occ1[s], cc_occ2[p], cc_occ2[q],
		    value, 0, outfile, dirac);
  else if(vir1[s] && occ1[r] && occ2[p] && occ2[q] && strcmp(spin,"AB"))
    iwl_buf_wrt_val(EBuf1, cc_vir1[s], cc_occ1[r], cc_occ2[p], cc_occ2[q], 
		    value, 0, outfile, dirac);

  if(!strcmp(spin,"AB")) {
    /* E (oo|ov) integrals */
    if(occ1[p] && occ1[q] && occ2[r] && vir2[s])
      iwl_buf_wrt_val(EBuf2, cc_occ1[p], cc_occ1[q], cc_occ2[r], cc_vir2[s],
		      value, 0, outfile, dirac);
    else if(occ1[p] && occ1[q] && occ2[s] && vir2[r])
      iwl_buf_wrt_val(EBuf2, cc_occ1[p], cc_occ1[q], cc_occ2[s], cc_vir2[r],
		      value, 0, outfile, dirac);
  }

  /* A (oo|oo) integrals */
  if(occ1[p] && occ1[q] && occ2[r] && occ2[s])
    iwl_buf_wrt_val(ABuf1, cc_occ1[p], cc_occ1[q], cc_occ2[r], cc_occ2[s],
		    value, 0, outfile, dirac);

}

}} // namespace psi::ccsort
