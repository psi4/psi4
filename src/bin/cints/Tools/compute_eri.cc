/*! \file compute_eri.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/

#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <libint/libint.h>
#include <libqt/qt.h>
#include "defines.h"
#define EXTERN
#include "global.h"
#include "norm_quartet.h"
#include "quartet_permutations.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"int_fjt.h"
  #include"fjt.h"
#endif
#include "quartet_data.h"
#include "symmetrize.h"
#include "small_fns.h"

namespace psi { 
  namespace CINTS {
    
    static void inline switch_ij(int& a, int& b) {int dum = a; a = b; b = dum;};
    
    int compute_eri(double* target,
		    Libint_t* Libint, int& si, int& sj, int& sk, int& sl,
		    int& inc1, int& inc2, int& inc3, int& inc4, const bool do_not_permute)
    {
      
#ifndef USE_TAYLOR_FM
      static double_array_t fjt_table;
      init_fjt_table(&fjt_table);
#endif
      
      /* place in "ascending" angular mom-
	 my simple way of optimizing PHG recursion (VRR) */
      
      bool switch_bra_ket = false;
      /* this should be /good/ for the VRR */
      if ( BasisSet.shells[si].am + BasisSet.shells[sj].am + inc1 + inc2 >
	   BasisSet.shells[sk].am + BasisSet.shells[sl].am + inc3 + inc4 ) {
	switch_bra_ket = true;
	switch_ij(si,sk);
	switch_ij(sj,sl);
	switch_ij(inc1, inc3);
	switch_ij(inc2, inc4);
      }
      
      /* these two are good for the HRR */
      bool switch_bra = false;
      if(BasisSet.shells[si].am + inc1 < BasisSet.shells[sj].am + inc2){
	switch_bra = true;
	switch_ij(si,sj);
	switch_ij(inc1,inc2);
      }
      bool switch_ket = false;
      if(BasisSet.shells[sk].am + inc3 < BasisSet.shells[sl].am + inc4){
	switch_ket = true;
	switch_ij(sk,sl);
	switch_ij(inc3,inc4);
      }
      
      int am1 = BasisSet.shells[si].am - 1 + inc1;
      int am2 = BasisSet.shells[sj].am - 1 + inc2;
      int am3 = BasisSet.shells[sk].am - 1 + inc3;
      int am4 = BasisSet.shells[sl].am - 1 + inc4;
      int am = am1 + am2 + am3 + am4;
      
      int c1 = BasisSet.shells[si].center - 1;
      int c2 = BasisSet.shells[sj].center - 1;
      int c3 = BasisSet.shells[sk].center - 1;
      int c4 = BasisSet.shells[sl].center - 1;
      
      // If all centers are the same and angular momentum of functions is odd -- the ERI will be 0 by AM selection
      if ( (c1 == c2) && (c1 == c3) && (c1 == c4) && (am%2 == 1)) {
	return 0;
      }
      
      struct shell_pair* sp_ij = &(BasisSet.shell_pairs[si][sj]);
      struct shell_pair* sp_kl = &(BasisSet.shell_pairs[sk][sl]);
      
      Libint->AB[0] = sp_ij->AB[0];
      Libint->AB[1] = sp_ij->AB[1];
      Libint->AB[2] = sp_ij->AB[2];
      Libint->CD[0] = sp_kl->AB[0];
      Libint->CD[1] = sp_kl->AB[1];
      Libint->CD[2] = sp_kl->AB[2];
      
      double AB2 = 
	Libint->AB[0]*Libint->AB[0]+
	Libint->AB[1]*Libint->AB[1]+
	Libint->AB[2]*Libint->AB[2];
      double CD2 = 
	Libint->CD[0]*Libint->CD[0]+
	Libint->CD[1]*Libint->CD[1]+
	Libint->CD[2]*Libint->CD[2];
      
      /*--- Compute data for primitive quartets here ---*/
      int np1 = BasisSet.shells[si].n_prims;
      int np2 = BasisSet.shells[sj].n_prims;
      int np3 = BasisSet.shells[sk].n_prims;
      int np4 = BasisSet.shells[sl].n_prims;
      
      /*--- Compute data for primitive quartets here ---*/
      int num_prim_comb = 0;
      //  bool si_eq_sj = (si == sj && am1 == am2);
      //  bool sk_eq_sl = (sk == sl && am3 == am4);
      bool si_eq_sj = false;
      bool sk_eq_sl = false;
      for (int p1 = 0; p1 < np1; p1++) {
	int max_p2 = si_eq_sj ? p1+1 : np2;
	for (int p2 = 0; p2 < max_p2; p2++) {
	  int m = (1 + (si_eq_sj && p1 != p2));
	  for (int p3 = 0; p3 < np3; p3++) {
	    int max_p4 = sk_eq_sl ? p3+1 : np4;
	    for (int p4 = 0; p4 < max_p4; p4++){
	      int n = m * (1 + (sk_eq_sl && p3 != p4));
#ifdef USE_TAYLOR_FM
	      quartet_data(&(Libint->PrimQuartet[num_prim_comb++]), NULL, AB2, CD2,
			   sp_ij, sp_kl, am, p1, p2, p3, p4, n);
#else
	      quartet_data(&(Libint->PrimQuartet[num_prim_comb++]), &fjt_table, AB2, CD2,
			   sp_ij, sp_kl, am, p1, p2, p3, p4, n);
#endif
	    }
	  }
	}
      }
      
      int nao[4];
      nao[0] = ioff[am1+1];
      nao[1] = ioff[am2+1];
      nao[2] = ioff[am3+1];
      nao[3] = ioff[am4+1];
      int nbra = nao[0] * nao[1];
      int nket = nao[2] * nao[3];
      int size = nbra * nket;
      
      // In general we need a couple extra buffers to be able to copy and permute integrals
      static int max_size = ioff[BasisSet.max_am]*ioff[BasisSet.max_am]*ioff[BasisSet.max_am]*ioff[BasisSet.max_am];
      static double* perm_buf = new double[max_size];
      static double* copy_buf = new double[max_size];
      if (size > max_size) {
	delete[] perm_buf;
	delete[] copy_buf;
	perm_buf = new double[size];
	copy_buf = new double[size];
	max_size = size;
      }
#ifdef NONDOUBLE_INTS
      double* raw_data = copy_buf;
#endif
      
      if (am) {
#ifdef NONDOUBLE_INTS
	REALTYPE* target_ptr = build_eri[am1][am2][am3][am4](Libint,num_prim_comb);
	for(i=0;i<size;i++)
	  raw_data[i] = (double) target_ptr[i];
#else
	double* raw_data = build_eri[am1][am2][am3][am4](Libint,num_prim_comb);
#endif
	
	//
	// Permute shells back, if necessary
	if (!do_not_permute) {
	  int num_perm = (int) switch_bra_ket + (int) switch_bra + (int) switch_ket;
	  
	  double *sbuf1, *sbuf2, *sbuf3;
	  double *tbuf1, *tbuf2, *tbuf3;
	  // in general, permutations will require careful orchestration of
	  // the utilization of the buffers. If no permutations are necessary -- just copy the data
	  if (num_perm == 0) {
	    for(int i=0; i<size; i++)
	      target[i] = raw_data[i];
	  }
	  // otherwise -- do the dance
	  else if (num_perm == 1) {
	    sbuf1 = sbuf2 = sbuf3 = raw_data;
	    tbuf1 = tbuf2 = tbuf3 = target;
	  }
	  else if (num_perm == 2) {
	    if (!switch_bra_ket) {
	      sbuf1 = raw_data;
	      tbuf1 = sbuf2 = perm_buf;
	      tbuf2 = target;
	    }
	    else if (!switch_ket) {
	      sbuf1 = raw_data;
	      tbuf1 = sbuf3 = perm_buf;
	      tbuf3 = target;
	    }
	    else if (!switch_bra) {
	      sbuf2 = raw_data;
	      tbuf2 = sbuf3 = perm_buf;
	      tbuf3 = target;
	    }
	  }
	  else if (num_perm == 3) {
	    sbuf1 = raw_data;
	    tbuf1 = sbuf2 = perm_buf;
	    tbuf2 = sbuf3 = raw_data;
	    tbuf3 = target;
	  }
	  
	  // Note that the order of permutations is opposite to the order of shell index permutations
	  // at the beginning of this function
	  if (switch_bra) {
	    ijkl_to_jikl(sbuf1, tbuf1, nao[0], nao[1], nket);
	    switch_ij(si,sj);
	    switch_ij(inc1,inc2);
	  }
	  if (switch_ket) {
	    ijkl_to_ijlk(sbuf2, tbuf2, nbra, nao[2], nao[3]);
	    switch_ij(sk,sl);
	    switch_ij(inc3,inc4);
	  }
	  if (switch_bra_ket) {
	    ijkl_to_klij(sbuf3, tbuf3, nbra, nket);
	    switch_ij(si,sk);
	    switch_ij(sj,sl);
	    switch_ij(inc1, inc3);
	    switch_ij(inc2, inc4);
	  }
	}
      }
      else {
	REALTYPE temp = 0.0;
	prim_data* prim_quartet = Libint->PrimQuartet;
	for(int p=0;p<num_prim_comb;p++,prim_quartet++) {
	  temp += prim_quartet->F[0];
	}
	target[0] = (double) temp;
      }
      
      return size;
    }
  }
}
