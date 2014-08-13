/* $Id: matmul.h,v 1.15.4.1 2006-12-22 13:05:22 manoj Exp $ */
#ifndef _MATMUL_H_
#define _MATMUL_H_

#include "ga.h"
#include "globalp.h"
#include "message.h"
#include "base.h"

#if HAVE_MATH_H
#   include <math.h>
#endif
#include "armci.h"

#include "galinalg.h"

/* min acceptable amount of memory (in elements) and default chunk size */
#  define MINMEM 64
#  define CHUNK_SIZE 256
#  define MAX_CHUNKS 1024
#  define BLOCK_SIZE 1024 /* temp buf size for pinning */
#  define GA_ASPECT_RATIO 3
#  define NUM_MATS 3 
#  define MINTASKS 10 /* increase this if there is high load imbalance */
#  define EXTRA 4

#define MIN_CHUNK_SIZE 256

#define SET   1
#define UNSET 0

extern void gai_matmul_patch_flag(int flag);

typedef struct {
  int lo[2]; /* 2 elements: ilo and klo */
  int hi[2];
  int dim[2];
  int chunkBId;
  short int do_put;
}task_list_t;

#define VECTORCHECK(rank,dims,dim1,dim2, ilo, ihi, jlo, jhi) \
  if(rank>2)  pnga_error("rank is greater than 2",rank); \
  else if(rank==2) {dim1=dims[0]; dim2=dims[1];} \
  else if(rank==1) {if((ihi-ilo)>0) { dim1=dims[0]; dim2=1;} \
                    else { dim1=1; dim2=dims[0];}} \
  else pnga_error("rank must be atleast 1",rank);

#define WAIT_GET_BLOCK(nbhdl) pnga_nbwait(nbhdl)

#endif /* _MATMUL_H_ */
