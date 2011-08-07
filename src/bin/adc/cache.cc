
/*
 *  cache.cc
 *  
 *
 *  Created by M.Saitow on 10/07/28.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <ccfiles.h>
#include <psifiles.h>

namespace psi { namespace adc {

    void cache_abcd_rhf(int **cachelist);
    void cache_iabc_rhf(int **cachelist);
    void cache_ijab_rhf(int **cachelist);
    void cache_iajb_rhf(int **cachelist);
    void cache_ijka_rhf(int **cachelist);
    void cache_ijkl_rhf(int **cachelist);	
    
    int **cacheprep_rhf(int level, int *cachefiles)
    {
      int **cachelist;
      
      /* The listing of CC files whose entries may be cached */
      cachefiles[CC_AINTS] = 1;
      cachefiles[CC_CINTS] = 1;
      cachefiles[CC_DINTS] = 1;
      cachefiles[CC_EINTS] = 1;
      cachefiles[CC_DENOM] = 1;
      cachefiles[CC_TAMPS] = 1;
      cachefiles[CC_LAMPS] = 1;
      cachefiles[CC_HBAR] = 1;
      
      /* The listing of DPD patterns which may be cached */
      cachelist = init_int_matrix(12,12);
      
      if(level == 0) return cachelist;
      else if(level == 1) {
	
	/*** Cache oooo and ooov ***/
	cache_ijkl_rhf(cachelist);
	cache_ijka_rhf(cachelist);
	
	return cachelist;
      }
      else if(level == 2) {
	
	/*** Cache oooo, ooov, oovv, and ovov ***/
	cache_ijkl_rhf(cachelist);
	cache_ijka_rhf(cachelist);
	cache_ijab_rhf(cachelist);
	cache_iajb_rhf(cachelist);
	
	return cachelist;
      }
      else if(level == 3) {
	
	/*** Cache, oooo, oov, oovv, ovov, and ovvv ***/
	
	cache_ijkl_rhf(cachelist);
	cache_ijka_rhf(cachelist);
	cache_ijab_rhf(cachelist);
	cache_iajb_rhf(cachelist);
	cache_iabc_rhf(cachelist);
	
	return cachelist;
      }
      else if(level == 4) {
	
	/*** Cache everything ***/
	cache_ijkl_rhf(cachelist);
	cache_ijka_rhf(cachelist);
	cache_ijab_rhf(cachelist);
	cache_iajb_rhf(cachelist);
	cache_iabc_rhf(cachelist);
	cache_abcd_rhf(cachelist);
	
	return cachelist;
      }
      else { printf("Error: invalid cache level!\n"); exit(PSI_RETURN_FAILURE); }
    }
    
    void cache_abcd_rhf(int **cachelist)
    {
      /* <ab|cd> */
      cachelist[5][5] = 1;
      cachelist[5][6] = 1;
      cachelist[5][7] = 1;
      cachelist[5][8] = 1;
      cachelist[5][9] = 1;
      cachelist[6][5] = 1;
      cachelist[6][6] = 1;
      cachelist[6][7] = 1;
      cachelist[6][8] = 1;
      cachelist[6][9] = 1;
      cachelist[7][5] = 1;
      cachelist[7][6] = 1;
      cachelist[7][7] = 1;
      cachelist[7][8] = 1;
      cachelist[7][9] = 1;
      cachelist[8][5] = 1;
      cachelist[8][6] = 1;
      cachelist[8][7] = 1;
      cachelist[8][8] = 1;
      cachelist[8][9] = 1;
      cachelist[9][5] = 1;
      cachelist[9][6] = 1;
      cachelist[9][7] = 1;
      cachelist[9][8] = 1;
      cachelist[9][9] = 1;
    }
    
    void cache_iabc_rhf(int **cachelist)
    {
      /* <ia|bc> */
      cachelist[10][5] = 1;
      cachelist[10][6] = 1;
      cachelist[10][7] = 1;
      cachelist[10][8] = 1;
      cachelist[10][9] = 1;
      cachelist[11][5] = 1;
      cachelist[11][6] = 1;
      cachelist[11][7] = 1;
      cachelist[11][8] = 1;
      cachelist[11][9] = 1;
      /* <ab|ci> */
      cachelist[5][10] = 1;
      cachelist[5][11] = 1;
      cachelist[6][10] = 1;
      cachelist[6][11] = 1;
      cachelist[7][10] = 1;
      cachelist[7][11] = 1;
      cachelist[8][10] = 1;
      cachelist[8][11] = 1;
      cachelist[9][10] = 1;
      cachelist[9][11] = 1;
    }
    
    void cache_ijab_rhf(int **cachelist)
    {
      /* <ij|ab> */
      cachelist[0][5] = 1;
      cachelist[0][6] = 1;
      cachelist[0][7] = 1;
      cachelist[0][8] = 1;
      cachelist[0][9] = 1;
      cachelist[1][5] = 1;
      cachelist[1][6] = 1;
      cachelist[1][7] = 1;
      cachelist[1][8] = 1;
      cachelist[1][9] = 1;
      cachelist[2][5] = 1;
      cachelist[2][6] = 1;
      cachelist[2][7] = 1;
      cachelist[2][8] = 1;
      cachelist[2][9] = 1;
      cachelist[3][5] = 1;
      cachelist[3][6] = 1;
      cachelist[3][7] = 1;
      cachelist[3][8] = 1;
      cachelist[3][9] = 1;
      cachelist[4][5] = 1;
      cachelist[4][6] = 1;
      cachelist[4][7] = 1;
      cachelist[4][8] = 1;
      cachelist[4][9] = 1;
      /* <ab|ij> */
      cachelist[5][0] = 1;
      cachelist[5][1] = 1;
      cachelist[5][2] = 1;
      cachelist[5][3] = 1;
      cachelist[5][4] = 1;
      cachelist[6][0] = 1;
      cachelist[6][1] = 1;
      cachelist[6][2] = 1;
      cachelist[6][3] = 1;
      cachelist[6][4] = 1;
      cachelist[7][0] = 1;
      cachelist[7][1] = 1;
      cachelist[7][2] = 1;
      cachelist[7][3] = 1;
      cachelist[7][4] = 1;
      cachelist[8][0] = 1;
      cachelist[8][1] = 1;
      cachelist[8][2] = 1;
      cachelist[8][3] = 1;
      cachelist[8][4] = 1;
      cachelist[9][0] = 1;
      cachelist[9][1] = 1;
      cachelist[9][2] = 1;
      cachelist[9][3] = 1;
      cachelist[9][4] = 1;
    }
    
    void cache_iajb_rhf(int **cachelist)
    {
      /* <ia|jb> */
      cachelist[10][10] = 1;
      cachelist[10][11] = 1;
      cachelist[11][10] = 1;
      cachelist[11][11] = 1;
    }
    
    void cache_ijka_rhf(int **cachelist)
    {
      /* <ij|ka> */
      cachelist[0][10] = 1;
      cachelist[0][11] = 1;
      cachelist[1][10] = 1;
      cachelist[1][11] = 1;
      cachelist[2][10] = 1;
      cachelist[2][11] = 1;
      cachelist[3][10] = 1;
      cachelist[3][11] = 1;
      cachelist[4][10] = 1;
      cachelist[4][11] = 1;
      /* <ia|jk> */
      cachelist[10][0] = 1;
      cachelist[10][1] = 1;
      cachelist[10][2] = 1;
      cachelist[10][3] = 1;
      cachelist[10][4] = 1;
      cachelist[11][0] = 1;
      cachelist[11][1] = 1;
      cachelist[11][2] = 1;
      cachelist[11][3] = 1;
      cachelist[11][4] = 1;
    }
    
    void cache_ijkl_rhf(int **cachelist)
    {
      /* <ij|kl> */
      cachelist[0][0] = 1;
      cachelist[0][1] = 1;
      cachelist[0][2] = 1;
      cachelist[0][3] = 1;
      cachelist[0][4] = 1;
      cachelist[1][0] = 1;
      cachelist[1][1] = 1;
      cachelist[1][2] = 1;
      cachelist[1][3] = 1;
      cachelist[1][4] = 1;
      cachelist[2][0] = 1;
      cachelist[2][1] = 1;
      cachelist[2][2] = 1;
      cachelist[2][3] = 1;
      cachelist[2][4] = 1;
      cachelist[3][0] = 1;
      cachelist[3][1] = 1;
      cachelist[3][2] = 1;
      cachelist[3][3] = 1;
      cachelist[3][4] = 1;
      cachelist[4][0] = 1;
      cachelist[4][1] = 1;
      cachelist[4][2] = 1;
      cachelist[4][3] = 1;
      cachelist[4][4] = 1;
    }
    
    void cachedone_rhf(int **cachelist)
    {
      free_int_matrix(cachelist);
    }

}}
