#ifndef _psi_src_bin_cints_Fock_hash_h
#define _psi_src_bin_cints_Fock_hash_h

/*! \file
    \ingroup CINTS
-----------------------------------
  Declarations of htable_entry, etc.
 -----------------------------------*/
#include <psitypes.h>

namespace psi { namespace CINTS {

typedef struct {
    PSI_INT_LEAST64 key;
    int si, sj, sk, sl;
    double q4ijkl, q4ikjl, q4ilkj;
} htable_entry;

typedef struct {
 htable_entry *table;
 int size;
} htable_t;

#define EMPTY_KEY -1

PSI_INT_LEAST64 compute_key(int si, int sj, int sk, int sl);
void init_htable(htable_t *htable, int nirreps);
void free_htable(htable_t *htable);
int put_entry(htable_t *htable, PSI_INT_LEAST64 key,
              int si, int sj, int sk, int sl,
	      double q4ijkl, double q4ikjl, double q4iljk);
int find_entry(htable_t *htable, PSI_INT_LEAST64 key);
}}
#endif
