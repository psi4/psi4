/*! \file hash.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<cstdlib>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"hash.h"

#define hashing_function(a) (a)%htable->size   /* division method */

namespace psi { namespace CINTS {
/*-------------------------
  Initialize hashing table
 -------------------------*/
void init_htable(htable_t *htable, int nirreps)
{
  int i;
  
  switch(nirreps) {
  case 1:
      htable->size = 7;
      break;
  case 2:
      htable->size = 97;
      break;
  case 4:
      htable->size = 1543;
      break;
  case 8:
  default:
      htable->size = 24593;
      break;
  }

  htable->table = (htable_entry *) malloc(htable->size*sizeof(htable_entry));
  for(i=0;i<htable->size;i++)
    htable->table[i].key = EMPTY_KEY;

  return;
}

      
/*-------------------------
  Deallocate hashing table
 -------------------------*/
void free_htable(htable_t *htable)
{
  free(htable->table);
  htable->table = NULL;
  return;
}


/*------------------------------------------
  Compute the key from four (shell) indices
 ------------------------------------------*/
PSI_INT_LEAST64 compute_key(int si, int sj, int sk, int sl)
{
  PSI_INT_LEAST64 ij, kl, ijkl;
  PSI_INT_LEAST64 key;

  ij = INDEX(si,sj);
  kl = INDEX(sk,sl);
  key = (ij > kl) ? ij*(ij+1)/2 + kl : kl*(kl+1)/2 + ij;
/*  key = INDEX(INDEX(si,sj),INDEX(sk,sl));*/

  return key;
}


/*----------------------------------------------------
  int put_entry() :  Put an entry into the hashing table.
                     Return the location if it's a new entry,
		     or -1 if it has been in the table.

  One must remember that given set of si, sj, sk, sl and
  the set found in the table may be a permutation of each
  other which belong to Q4 - group which consists of all
  permutations under which an ERI is unchanged
  (Q4 is a subgroup of T4 of order 8).
  Another complication
  is that if si == sj or sk == sl - only htable[].q4ikjl
  is relevant, yet due to possibility of permutations
  it is safer to just increment it by THE SUM of q4ikjl and
  q4ilkj (since one of them is zero);

  Let P be the permutation which transforms the given {si,sj,sk,sl}
  into htable[].{si,sj,sk,sl}. P must be a product of the following
  (commuting) operations P0 (trivial transposition), P12 (swaps
  si and sj), P34 (swaps sk and sl), and P12,34 = P13 P24 (swaps bra
  and ket). Let K be a transposition that describes the contribution
  of a given {si,sj,sk,sl} to the current P(si,sj|sk,sl). K must be
  either P0 (q4ijkl is non-zero), P23 (q4ikjl non-zero) or P24
  (q4ilkj non-zero). The following set of rules may be easily derived
  using multiplication table of T4:

  P0  {P} equiv(Q4) P0
  P23 P12 equiv(Q4) P24
  P24 P12 equiv(Q4) P23
  P23 P34 equiv(Q4) P24
  P24 P34 equiv(Q4) P23
  P23 P12,34 equiv(Q4) P23
  P24 P12,34 equiv(Q4) P24

  equiv(Q4) means equivalent up to any permutation from Q4.

  This means that P is either P12 or P34 we need to swap q4ikjl and q4ilkj.

 ----------------------------------------------------*/
int put_entry(htable_t *htable, PSI_INT_LEAST64 key, int si, int sj, int sk, int sl, double q4ijkl, double q4ikjl, double q4ilkj)
{
  int curr_ptr;
  int hvalue = hashing_function(key);
  int return_code = -1;
  int P_includes_P12, P_includes_P34;
  
  curr_ptr = hvalue;
  while(htable->table[curr_ptr].key != key && htable->table[curr_ptr].key != EMPTY_KEY) {
    curr_ptr++;
    if (htable->size == curr_ptr)
      curr_ptr = 0;
  }
  if (htable->table[curr_ptr].key == EMPTY_KEY) {
    return_code = curr_ptr;
    htable->table[curr_ptr].key = key;
    htable->table[curr_ptr].si = si;
    htable->table[curr_ptr].sj = sj;
    htable->table[curr_ptr].sk = sk;
    htable->table[curr_ptr].sl = sl;
    htable->table[curr_ptr].q4ijkl = q4ijkl;
    if (si == sj || sk == sl) {
      htable->table[curr_ptr].q4ikjl = q4ikjl + q4ilkj;
      htable->table[curr_ptr].q4ilkj = 0;
    }
    else {
      htable->table[curr_ptr].q4ikjl = q4ikjl;
      htable->table[curr_ptr].q4ilkj = q4ilkj;
    }
  }
  else {
    htable->table[curr_ptr].q4ijkl += q4ijkl;
    if (si == sj || sk == sl) {
      htable->table[curr_ptr].q4ikjl += q4ikjl + q4ilkj;
    }
    else {
      P_includes_P34 = ( (htable->table[curr_ptr].sk == sl && htable->table[curr_ptr].sl == sk) ||
			 (htable->table[curr_ptr].sk == sj && htable->table[curr_ptr].sl == si) );
      P_includes_P12 = ( (htable->table[curr_ptr].si == sj && htable->table[curr_ptr].sj == si) ||
			 (htable->table[curr_ptr].si == sl && htable->table[curr_ptr].sj == sk) );
      if (P_includes_P12 ^ P_includes_P34) {
	htable->table[curr_ptr].q4ikjl += q4ilkj;
	htable->table[curr_ptr].q4ilkj += q4ikjl;
      }
      else {
	htable->table[curr_ptr].q4ikjl += q4ikjl;
	htable->table[curr_ptr].q4ilkj += q4ilkj;
      }
    }
  }

  return return_code;
}

/*----------------------------------------------------
  Put an entry into the table. Return the location if
  it's a new entry, or -1 if it has been in the table
 ----------------------------------------------------*/
int find_entry(htable_t *htable, PSI_INT_LEAST64 key)
{
  int curr_ptr;
  int hvalue = hashing_function(key);
  int return_code = -1;
  
  curr_ptr = hvalue;
  while(htable->table[curr_ptr].key != key && htable->table[curr_ptr].key != EMPTY_KEY) {
    curr_ptr++;
    if (htable->size == curr_ptr)
      curr_ptr = 0;
  }
  if (htable->table[curr_ptr].key == EMPTY_KEY)
    return -1;
  else
    return curr_ptr;

}
};};
