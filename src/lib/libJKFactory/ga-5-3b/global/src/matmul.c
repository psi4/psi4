#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: matmul.c,v 1.60.4.1 2006-12-22 13:05:22 manoj Exp $ */
/*===========================================================
 *
 *         GA_Dgemm(): Parallel Matrix Multiplication
 *              (i.e.  C = alpha*A*B + beta*C)
 *
 *===========================================================*/

#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#include "matmul.h"
#include "ga-papi.h"
#include "ga-wapi.h"

#define DEBUG_ 0 /*set 1, to verify the correctness of parallel matrix mult.*/

/* some optimization macros */
#define KCHUNK_OPTIMIZATION 0 /* This Opt performing well only for m=1000;n=1000'k=2000 kinda cases and not for the opposite*/

/* Optimization flags: Initialized everytime in pnga_matmul() */
static short int CYCLIC_DISTR_OPT_FLAG  = SET;
static short int CONTIG_CHUNKS_OPT_FLAG = SET;
static short int DIRECT_ACCESS_OPT_FLAG = SET;

Integer gNbhdlA[2], gNbhdlB[2], gNbhdlC[2];/* for A and B matrix */

static int _gai_matmul_patch_flag = 0;
void gai_matmul_patch_flag(int flag)
{
    _gai_matmul_patch_flag = flag;
}

static inline int max3(int ichunk, int jchunk, int kchunk) {
  if(ichunk>jchunk) return GA_MAX(ichunk,kchunk);
  else return GA_MAX(jchunk, kchunk);
}

static inline void init_task_list(task_list_t *thing)
{
    thing->lo[0] = 0;
    thing->lo[1] = 0;
    thing->hi[0] = 0;
    thing->hi[1] = 0;
    thing->dim[0] = 0;
    thing->dim[1] = 0;
    thing->chunkBId = 0;
    thing->do_put = 0;
}

static void GET_BLOCK(Integer g_x, task_list_t *chunk, void *buf, 
		      char *trans, Integer xilo, Integer xjlo, 
		      Integer *dim_next, Integer *nbhdl) {

    Integer i0, i1, j0, j1;
    Integer lo[2], hi[2];

    if(*trans == 'n' || *trans == 'N') {
       *dim_next = chunk->dim[0];
       i0= xilo+chunk->lo[0]; i1= xilo+chunk->hi[0];
       j0= xjlo+chunk->lo[1]; j1= xjlo+chunk->hi[1];
    }
    else {
       *dim_next = chunk->dim[1];
       i0= xjlo+chunk->lo[1]; i1= xjlo+chunk->hi[1];
       j0= xilo+chunk->lo[0]; j1= xilo+chunk->hi[0];
    }

    lo[0] = i0;
    lo[1] = j0;
    hi[0] = i1;
    hi[1] = j1;
    pnga_nbget(g_x, lo, hi, buf, dim_next, nbhdl);
}

static short int
gai_get_task_list(task_list_t *taskListA, task_list_t *taskListB, 
		  task_list_t *state, Integer istart, Integer jstart,
		  Integer kstart, Integer iend, Integer jend, Integer kend, 
		  Integer Ichunk, Integer Jchunk, Integer Kchunk, 
		  int *max_tasks, Integer g_a) {
    
    int ii, jj, nloops=0;
    short int do_put, more_chunks_left=0, recovery=0;
    Integer ilo, ihi, jlo, jhi, klo, khi, get_new_B;
    Integer jstart_=jstart, kstart_=kstart;
    
    if(state->lo[0] != -1) recovery = 1;

    nloops = (iend-istart+1)/Ichunk + ( ((iend-istart+1)%Ichunk)?1:0 );
    if(nloops>MAX_CHUNKS) pnga_error("Increase MAX_CHUNKS value in matmul.h",0L);

    if(recovery) jstart_ = state->lo[0]; /* recovering the previous state */
    for(ii=jj=0, jlo = jstart_; jlo <= jend; jlo += Jchunk) {
       jhi = GA_MIN(jend, jlo+Jchunk-1);

       if(recovery) {
	  do_put = state->do_put;
	  kstart_ =  state->lo[1];
       }
       else do_put = SET; /* for 1st shot we can put, instead of accumulate */
       
       for(klo = kstart_; klo <= kend; klo += Kchunk) {
	  khi = GA_MIN(kend, klo+Kchunk-1); 
	  get_new_B = TRUE;
	  
	  /* set it back after the first loop */
	  recovery = 0;
	  jstart_ = jstart;
	  kstart_ = kstart;
	  
	  /* save CURRENT STATE. Saving state before "i" loop helps to avoid 
	     tracking get_new_B, which is hassle in ga_matmul_regular() */
	  if(ii+nloops >= MAX_CHUNKS) {
	     more_chunks_left = 1;
	     state->lo[0]  = jlo;
	     state->lo[1]  = klo;
	     state->do_put   = do_put;
	     break;
	  }
	  
	  for(ilo = istart; ilo <= iend; ilo += Ichunk){ 	     
	     ihi = GA_MIN(iend, ilo+Ichunk-1);
	     taskListA[ii].dim[0] = ihi - ilo + 1; 
	     taskListA[ii].dim[1] = khi - klo + 1;
	     taskListA[ii].lo[0]  = ilo; taskListA[ii].hi[0] = ihi;
	     taskListA[ii].lo[1]  = klo; taskListA[ii].hi[1] = khi;
	     taskListA[ii].do_put = do_put;
	     if(get_new_B) { /* B matrix */
		ihi = GA_MIN(iend, ilo+Ichunk-1);
		taskListB[jj].dim[0] = khi - klo + 1; 
		taskListB[jj].dim[1] = jhi - jlo + 1;
		taskListB[jj].lo[0]  = klo; taskListB[jj].hi[0] = khi;
		taskListB[jj].lo[1]  = jlo; taskListB[jj].hi[1] = jhi;
		get_new_B = FALSE; /* Until J or K change again */
		taskListA[ii].chunkBId = jj;
		++jj;
	     }
	     else taskListA[ii].chunkBId = taskListA[ii-1].chunkBId;
	     ++ii;
	  }
	  if (more_chunks_left) break;
	  do_put = UNSET;
       }
       if (more_chunks_left) break;
    }

    *max_tasks = ii;

    /* Optimization disabled if chunks exceeds buffer space */
    if(more_chunks_left) CYCLIC_DISTR_OPT_FLAG = UNSET;

    if(CYCLIC_DISTR_OPT_FLAG) { /* should not be called for irregular matmul */
       int prow, pcol, offset, grp_me;
       Integer a_grp = pnga_get_pgroup(g_a);
       grp_me = (int)pnga_pgroup_nodeid(a_grp);
       prow = GA[GA_OFFSET + g_a].nblock[0];
       pcol = GA[GA_OFFSET + g_a].nblock[1];
       offset = (grp_me/prow + grp_me%prow) % pcol;
       for(jj=0, ilo = istart; ilo <= iend; jj++, ilo += Ichunk)
	  taskListA[jj].do_put = UNSET;
       for(jj=0, ilo = istart; ilo <= iend; jj++, ilo += Ichunk)
	  taskListA[jj+offset].do_put = SET;
    }

    return more_chunks_left;
}

static void gai_get_chunk_size(int irregular,Integer *Ichunk,Integer *Jchunk,
			       Integer *Kchunk,Integer *elems,Integer atype, 
			       Integer m,Integer n,Integer k, short int nbuf,
			       short int use_armci_memory, Integer a_grp) {
    double temp;
    Integer min_tasks = MINTASKS; /* Increase tasks if there is load imbalance.
				     This controls the granularity of chunks */
    Integer  max_chunk, nproc=pnga_nnodes(), tmpa, tmpb;
    Integer avail = pnga_memory_avail_type(atype);

    tmpa = *Ichunk;
    tmpb = *Jchunk;
    
    if(irregular) {
       temp = (k*(double)(m*(double)n)) / (min_tasks * nproc);
       max_chunk = (Integer)pow(temp, (1.0/3.0) );
       if (max_chunk < MIN_CHUNK_SIZE) max_chunk = MIN_CHUNK_SIZE;  
    }
    else
       max_chunk = (Integer) max3(*Ichunk, *Jchunk, *Kchunk);

    pnga_pgroup_gop(a_grp, pnga_type_f2c(MT_F_INT), &avail, (Integer)1, "min");
    
    if ( max_chunk > CHUNK_SIZE/nbuf) {
       /*if memory if very limited, performance degrades for large matrices
	 as chunk size is very small, which leads to communication overhead)*/
      if(avail<MINMEM && pnga_pgroup_nodeid(a_grp)==0) pnga_error("NotEnough memory",avail);
      *elems = (Integer)(avail*0.9); /* Donot use every last drop */
      
      /* MAX: get the maximum chunk (or, block) size i.e  */
      max_chunk=GA_MIN(max_chunk, (Integer)(sqrt( (double)((*elems-nbuf*NUM_MATS)/(nbuf*NUM_MATS)))));

      if(!irregular && use_armci_memory==SET) 
	 max_chunk = *Ichunk = *Jchunk = *Kchunk = BLOCK_SIZE;
    
      if(irregular) {
	 /* NOTE:enable this part for regular cases, if later 
	    part of the code is buggy or inefficient */
	 *Ichunk = GA_MIN(m,max_chunk);
	 *Jchunk = GA_MIN(n,max_chunk);
	 *Kchunk = GA_MIN(k,max_chunk);      
      }
      else { /* This part of the code takes care of rectangular chunks and
		most probably gives optimum rectangular chunk size */
	 temp = max_chunk*max_chunk;
	 if(*Ichunk < max_chunk && *Kchunk > max_chunk) {
	    *Kchunk = GA_MIN(*Kchunk,(Integer)(temp/(*Ichunk)));
	    *Jchunk = GA_MIN(*Jchunk,(Integer)(temp/(*Kchunk)));
	 }
	 else if(*Kchunk < max_chunk && *Ichunk > max_chunk) {
	    temp *= 1.0/(*Kchunk);
	    *Ichunk = GA_MIN(*Ichunk,(Integer)temp);
	    *Jchunk = GA_MIN(*Jchunk,(Integer)temp);
	 }
	 else *Ichunk = *Jchunk = *Kchunk = max_chunk;
      }
    }
    else 
       *Ichunk = *Jchunk = *Kchunk = CHUNK_SIZE/nbuf;

    /* Try to use 1-d data transfer & take advantage of zero-copy protocol */
    if(CONTIG_CHUNKS_OPT_FLAG) { /* select a contiguous piece */
       if(!irregular) {
	  if(*Ichunk > tmpa && *Jchunk > tmpb) {
	     *Ichunk = tmpa;
	     *Jchunk = tmpb;
	     *Kchunk = GA_MIN(*Ichunk,*Jchunk);
	  }
	  else {
	     int i=1;/* i should be >=1 , to avoid divide by zero error */
	     temp = max_chunk*max_chunk;
	     if(temp > tmpa) {
		*Ichunk = tmpa;
		*Jchunk = (Integer)(temp/(*Ichunk));
		if(*Jchunk < tmpb) {
		   while(tmpb/i > *Jchunk) ++i;
		   *Jchunk = tmpb/i;
		}
		else *Jchunk = tmpb;
		*Kchunk = GA_MIN(*Ichunk, *Jchunk);
	     }
	  }
       }
    }

    if(*Ichunk<=0) *Ichunk = 1; /* should be atleast 1 */
    if(*Jchunk<=0) *Jchunk = 1;
    if(*Kchunk<=0) *Kchunk = 1;

    /* Total elements "NUM_MAT" extra elems for safety - just in case */
    *elems = ( nbuf*(*Ichunk)*(*Kchunk) + nbuf*(*Kchunk)*(*Jchunk) + 
	       (*Ichunk)*(*Jchunk) );
    *elems += nbuf*NUM_MATS*sizeof(DoubleComplex)/GAsizeofM(atype);
}

static DoubleComplex* 
gai_get_armci_memory(Integer Ichunk, Integer Jchunk, Integer Kchunk,
		     short int nbuf, Integer atype) {

    DoubleComplex *tmp = NULL;
    Integer elems;

    elems = (Integer) pow((double)BLOCK_SIZE,(double)2);
    elems = nbuf*elems + nbuf*elems + elems; /* A,B,C temporary buffers */
    
    /* add extra elements for safety */
    elems += nbuf*NUM_MATS*sizeof(DoubleComplex)/GAsizeofM(atype);

    /* allocate temporary storage using ARMCI_Malloc */
    if( (Integer) (((double)nbuf)*(Ichunk* Kchunk) + 
		   ((double)nbuf)*(Kchunk* Jchunk) + 
		   Ichunk* Jchunk ) < elems) {
       tmp=(DoubleComplex*)ARMCI_Malloc_local(elems*GAsizeofM(atype));
    }
    return tmp;
}
	  
/************************************
 * Sequential DGEMM 
 *      i.e. BLAS dgemm Routines
 ************************************/

static void GAI_DGEMM(Integer atype, char *transa, char *transb, 
        Integer idim, Integer jdim, Integer kdim, void *alpha, 
        DoubleComplex *a, Integer adim, DoubleComplex *b, 
        Integer bdim, DoubleComplex *c, Integer cdim) {

    BlasInt idim_t, jdim_t, kdim_t, adim_t, bdim_t, cdim_t;
    DoubleComplex ZERO;
    SingleComplex ZERO_CF;

    idim_t=idim; jdim_t=jdim; kdim_t=kdim;
    adim_t=adim; bdim_t=bdim; cdim_t=cdim;
    ZERO.real = 0.; ZERO.imag = 0.;
    ZERO_CF.real = 0.; ZERO_CF.imag = 0.;

    switch(atype) {
        case C_FLOAT:
            BLAS_SGEMM(transa, transb, &idim_t, &jdim_t, &kdim_t,
                    (Real *)alpha, (Real *)a, &adim_t,
                    (Real *)b, &bdim_t, (Real *)&ZERO_CF,
                    (Real *)c, &cdim_t);
            break;
        case C_DBL:
            BLAS_DGEMM(transa, transb, &idim_t, &jdim_t, &kdim_t,
                    (DoublePrecision *)alpha, (DoublePrecision *)a, &adim_t,
                    (DoublePrecision *)b, &bdim_t, (DoublePrecision *)&ZERO,
                    (DoublePrecision *)c, &cdim_t);
            break;
        case C_DCPL:
            BLAS_ZGEMM(transa, transb, &idim_t, &jdim_t, &kdim_t,
                    (DoubleComplex *)alpha, (DoubleComplex *)a, &adim_t,
                    (DoubleComplex *)b, &bdim_t, (DoubleComplex *)&ZERO,
                    (DoubleComplex *)c, &cdim_t);
            break;
        case C_SCPL:
            BLAS_CGEMM(transa, transb, &idim_t, &jdim_t, &kdim_t,
                    (SingleComplex *)alpha, (SingleComplex *)a, &adim_t,
                    (SingleComplex *)b, &bdim_t, (SingleComplex *)&ZERO_CF,
                    (SingleComplex *)c, &cdim_t);
            break;
        default:
            pnga_error("ga_matmul_patch: wrong data type", atype);
    }
}



static void gai_matmul_shmem(transa, transb, alpha, beta, atype,
			     g_a, ailo, aihi, ajlo, ajhi,
			     g_b, bilo, bihi, bjlo, bjhi,
			     g_c, cilo, cihi, cjlo, cjhi,
			     Ichunk, Kchunk, Jchunk, a,b,c, need_scaling)
     
     Integer g_a, ailo, aihi, ajlo, ajhi;    /* patch of g_a */
     Integer g_b, bilo, bihi, bjlo, bjhi;    /* patch of g_b */
     Integer g_c, cilo, cihi, cjlo, cjhi;    /* patch of g_c */
     Integer Ichunk, Kchunk, Jchunk, atype;
     void    *alpha, *beta;
     char    *transa, *transb;
     DoubleComplex *a, *b, *c;
     short int need_scaling;
{

  Integer me = pnga_nodeid();
  Integer get_new_B, loC[2]={0,0}, hiC[2]={0,0}, ld[2];
  Integer i0, i1, j0, j1;
  Integer ilo, ihi, idim, jlo, jhi, jdim, klo, khi, kdim, adim, bdim=0, cdim;
  int istart, jstart, kstart, iend, jend, kend;
  short int do_put=UNSET, single_task_flag=UNSET;
  DoubleComplex ONE = {1.,0.};
  SingleComplex ONE_CF = {1.,0.};
  Integer clo[2], chi[2];

  GA_PUSH_NAME("ga_matmul_shmem");

  /* to skip accumulate and exploit data locality:
     get chunks according to "C" matrix distribution*/
  pnga_distribution(g_c, me, loC, hiC);
  istart = loC[0]-1; iend = hiC[0]-1;
  jstart = loC[1]-1; jend = hiC[1]-1;
  kstart = 0       ; kend = ajhi-ajlo;

  if(DIRECT_ACCESS_OPT_FLAG) {
    /* check if there is only one task. If so, then it is contiguous */
    if( (iend-istart+1 <= Ichunk) && (jend-jstart+1 <= Jchunk) &&
        (kend-kstart+1 <= Kchunk) ) {
      single_task_flag = SET;
      pnga_access_ptr(g_c, loC, hiC, &c, ld);
    }
  }

  /* loop through columns of g_c patch */
  for(jlo = jstart; jlo <= jend; jlo += Jchunk) { 
    jhi  = GA_MIN(jend, jlo+Jchunk-1);
    jdim = jhi - jlo +1;

    /* if beta=0,then for first shot we can do put,instead of accumulate */
    if(need_scaling == UNSET) do_put = SET;

    /* loop cols of g_a patch : loop rows of g_b patch*/
    for(klo = kstart; klo <= kend; klo += Kchunk) { 
      khi = GA_MIN(kend, klo+Kchunk-1);
      kdim= khi - klo +1;
      get_new_B = TRUE; /* Each pass thru' outer 2 loops means we 
                           need a different patch of B.*/
      /*loop through rows of g_c patch */
      for(ilo = istart; ilo <= iend; ilo += Ichunk){ 
        ihi = GA_MIN(iend, ilo+Ichunk-1);
        idim= cdim = ihi - ilo +1;

        /* STEP1(a): get matrix "A" chunk */
        i0= ailo+ilo; i1= ailo+ihi;
        j0= ajlo+klo; j1= ajlo+khi;
        if (*transa == 'n' || *transa == 'N'){
          clo[0] = i0;
          clo[1] = j0;
          chi[0] = i1;
          chi[1] = j1;
          adim=idim; pnga_get(g_a, clo, chi, a, &idim);
        }else{
          clo[0] = j0;
          clo[1] = i0;
          chi[0] = j1;
          chi[1] = i1;
          adim=kdim; pnga_get(g_a, clo, chi, a, &kdim);
        }

        /* STEP1(b): get matrix "B" chunk*/
        if(get_new_B) {/*Avoid rereading B if same patch as last time*/
          i0= bilo+klo; i1= bilo+khi;
          j0= bjlo+jlo; j1= bjlo+jhi;
          if (*transb == 'n' || *transb == 'N'){ 
            clo[0] = i0;
            clo[1] = j0;
            chi[0] = i1;
            chi[1] = j1;
            bdim=kdim; pnga_get(g_b, clo, chi, b, &kdim);  
          }else {
            clo[0] = j0;
            clo[1] = i0;
            chi[0] = j1;
            chi[1] = i1;
            bdim=jdim; pnga_get(g_b, clo, chi, b, &jdim);
          }
          get_new_B = FALSE; /* Until J or K change again */
        }

        /* STEP2: Do the sequential matrix multiply - i.e.BLAS dgemm */
        GAI_DGEMM(atype, transa, transb, idim, jdim, kdim, alpha, 
            a, adim, b, bdim, c, cdim);

        /* STEP3: put/accumulate into "C" matrix */
        i0= cilo+ilo; i1= cilo+ihi;   
        j0= cjlo+jlo; j1= cjlo+jhi;	 
        /* if single_task_flag is SET (i.e =1), then there is no need to 
           update "C" matrix, as we use pointer directly in GAI_DGEMM */
        if(single_task_flag != SET) {
          switch(atype) {
            case C_FLOAT:
            case C_SCPL:
              clo[0] = i0;
              clo[1] = j0;
              chi[0] = i1;
              chi[1] = j1;
              if(do_put==SET) /* i.e.beta == 0.0 */
                pnga_put(g_c, clo, chi, (float *)c, &cdim);
              else {
                pnga_acc(g_c, clo, chi, (float*)c, &cdim, &ONE_CF);
              }
              break;
            default:
              clo[0] = i0;
              clo[1] = j0;
              chi[0] = i1;
              chi[1] = j1;
              if(do_put==SET) /* i.e.beta == 0.0 */
                pnga_put(g_c, clo, chi, (DoublePrecision*)c, &cdim);
              else {
                pnga_acc(g_c, clo, chi, (DoublePrecision*)c,
                    &cdim, (DoublePrecision*)&ONE);
              }
              break;
          }
        }
      }
      do_put = UNSET; /* In the second loop, accumulate should be done */
    }
  }
  GA_POP_NAME;
}


static
void init_block_info(Integer g_c, Integer *proc_index, Integer *index,
                     Integer *blocks, Integer *block_dims, Integer *topology,
                     Integer *iblock) 
{
    Integer me= pnga_nodeid();

    /* Uses simple block-cyclic data distribution */
    if(!pnga_uses_proc_grid(g_c))
    {
       *iblock = me;
    }
    else /* Uses scalapack block-cyclic data distribution */ 
    {   
       *iblock = 0;
       pnga_get_proc_index(g_c, me, proc_index);
       pnga_get_proc_index(g_c, me, index);
       pnga_get_block_info(g_c, blocks, block_dims);
       pnga_get_proc_grid(g_c, topology);
    }    
}

/**
 * get the lo/hi distribution info of the next block
 *   return 0 indicates no more blocks available
 *   return 1 indicates there is a block available
 */
static
int get_next_block_info(Integer g_c, Integer *proc_index, Integer *index,
                        Integer *blocks, Integer *block_dims, Integer*topology,
                        Integer *iblock, Integer *blo, Integer *bhi) 
{
    Integer dims[MAXDIM], ndim, type;
    int i;
    
    /* works only upto 2 dims - i.e vectors/matrices*/
    pnga_inquire(g_c,  &type, &ndim, dims);
    if(ndim>2) pnga_error("get_next_block_info() supports upto 2-d only ", 0L);
    
    /* Uses simple block-cyclic data distribution */
    if (!pnga_uses_proc_grid(g_c)) 
    {
       if(*iblock < pnga_total_blocks(g_c)) 
       {
          pnga_distribution(g_c, *iblock, blo, bhi);
          *iblock += pnga_nnodes();
          return 1;
       }
    }
    else /* Uses scalapack block-cyclic data distribution */
    {
       if (index[ndim-1] < blocks[ndim-1]) 
       {
          /* find bounding coordinates of block */
          for (i = 0; i < ndim; i++) 
          {
             blo[i] = index[i]*block_dims[i]+1;
             bhi[i] = (index[i] + 1)*block_dims[i];
             if (bhi[i] > dims[i]) bhi[i] = dims[i];
          }

          /* increment index to get next block on processor */
          index[0] += topology[0];
          for (i = 0; i < ndim; i++) 
          {
             if (index[i] >= blocks[i] && i<ndim-1) 
             {
                index[i]    = proc_index[i];
                index[i+1] += topology[i+1];
             }
          }

          return 1;
       } /* end of while */ 
    }

    return 0; 
}

                        

static void gai_matmul_regular(transa, transb, alpha, beta, atype,
			       g_a, ailo, aihi, ajlo, ajhi,
			       g_b, bilo, bihi, bjlo, bjhi,
			       g_c, cilo, cihi, cjlo, cjhi,
			       Ichunk, Kchunk, Jchunk, a_ar,b_ar,c_ar, 
			       need_scaling, irregular) 
     
     Integer g_a, ailo, aihi, ajlo, ajhi;    /* patch of g_a */
     Integer g_b, bilo, bihi, bjlo, bjhi;    /* patch of g_b */
     Integer g_c, cilo, cihi, cjlo, cjhi;    /* patch of g_c */
     Integer Ichunk, Kchunk, Jchunk, atype;
     void    *alpha, *beta;
     char    *transa, *transb;
     DoubleComplex **a_ar, **b_ar, **c_ar;
     short int need_scaling, irregular;
{

  Integer me= pnga_nodeid();
  Integer get_new_B=TRUE, i0, i1, j0, j1;
  Integer idim, jdim, kdim;
  Integer k, adim=0, bdim, cdim, adim_next, bdim_next;
  Integer clo[2], chi[2], loC[2]={1,1}, hiC[2]={1,1}, ld[2];
  int max_tasks=0, shiftA=0, shiftB=0;
  int currA, nextA, currB, nextB=0; /* "current" and "next" task Ids */
  task_list_t taskListA[MAX_CHUNKS], taskListB[MAX_CHUNKS], state; 
  short int do_put=UNSET, single_task_flag=UNSET, chunks_left=0;
  DoubleComplex ONE, *a, *b, *c;
  SingleComplex ONE_CF;
  int offset=0, gTaskId=0;
  int numblocks=0, has_more_blocks=1;
  Integer ctype, cndim, cdims[2];
  Integer iblock=0, proc_index[2], index[2];
  Integer blocks[2], block_dims[2], topology[2];

  GA_PUSH_NAME("ga_matmul_regular");
  if(irregular) pnga_error("irregular flag set", 0L);

  init_task_list(&state);

#if DEBUG_
  if(me==0) { printf("@@ga_matmul_regular:m,n,k=%ld %ld %ld\n",aihi-ailo+1,
      bjhi-bjlo+1,ajhi-ajlo+1);  fflush(stdout); }
#endif

  ONE.real =1.;    ONE.imag =0.;   
  ONE_CF.real =1.; ONE_CF.imag =0.;   
  clo[0] = cilo; chi[0] = cihi;
  clo[1] = cjlo; chi[1] = cjhi;
  k = ajhi - ajlo +1;

  numblocks = pnga_total_blocks(g_c);
  if(numblocks>=0) init_block_info(g_c, proc_index, index, blocks,
      block_dims, topology, &iblock);

  has_more_blocks = 1;
  while(has_more_blocks) 
  {
    /* if block cyclic distribution and loop accordigly. In case of simple
       block distribution we process this loop only once */
    if(numblocks<0)
    { /* simple block distribution */
      has_more_blocks = 0; 
      pnga_distribution(g_c, me, loC, hiC);
    }
    else
    { /* block cyclic */

      if(!get_next_block_info(g_c, proc_index, index, blocks, block_dims,
            topology, &iblock, loC, hiC))
        break;
    }

    /* If loC and hiC intersects with current patch region, then they will
     * be updated accordingly. Else it returns FALSE */
    pnga_inquire(g_c, &ctype, &cndim, cdims);
    if(!pnga_patch_intersect(clo,chi,loC,hiC,cndim)) continue;

#if DEBUG_
    printf("%d: Processing block #%d [%d,%d] - [%d,%d]\n", GAme, iblock,
        loC[0], loC[1], hiC[0], hiC[1]);
#endif

    state.lo[0] = -1; /* just for first do-while loop */
    do {

      /* Inital Settings */
      a = a_ar[0];
      b = b_ar[0];
      c = c_ar[0];
      do_put = single_task_flag = UNSET;
      offset = 0;

      /*****************************************************************
       * Task list: Collect information of all chunks. Matmul using 
       * Non-blocking call needs this list 
       *****************************************************************/
      gTaskId=0;

      /* to skip accumulate and exploit data locality:
         get chunks according to "C" matrix distribution*/
      /* pnga_distribution(g_c, me, loC, hiC); */
      chunks_left=gai_get_task_list(taskListA, taskListB, &state,loC[0]-1,
          loC[1]-1, 0, hiC[0]-1, hiC[1]-1, k-1,
          Ichunk,Jchunk,Kchunk, &max_tasks,g_a);
      currA = nextA = 0;

      if(chunks_left) { /* then turn OFF this optimization */
        if(DIRECT_ACCESS_OPT_FLAG) {
          /* check if there is only one task.If so,then it is contiguous */
          if(max_tasks == 1) {
            if( !((hiC[0]-loC[0]+1 <= Ichunk) &&(hiC[1]-loC[1]+1 <=Jchunk)
                  && (k <= Kchunk))) 
              pnga_error("Invalid task list", 0L);
            single_task_flag = SET;
            pnga_access_ptr(g_c, loC, hiC, &c, ld);
          }
        }
      }

      if(CYCLIC_DISTR_OPT_FLAG) {
        int prow,pcol,grp_me;
        Integer a_grp=pnga_get_pgroup(g_a);
        grp_me = (int)pnga_pgroup_nodeid(a_grp);
        prow = GA[GA_OFFSET + g_a].nblock[0];
        pcol = GA[GA_OFFSET + g_a].nblock[1];
        offset = (grp_me/prow + grp_me%prow) % pcol;
        currA = nextA = nextA + offset;
      }

      /*************************************************
       * Do the setup & issue non-blocking calls to get 
       * the first block/chunk I'm gonna work 
       *************************************************/
      shiftA=0; shiftB=0;
      if(nextA < max_tasks) {
        currB = nextB = taskListA[currA].chunkBId;

        GET_BLOCK(g_a, &taskListA[nextA], a_ar[shiftA], transa, 
            ailo, ajlo, &adim_next, &gNbhdlA[shiftA]);

        GET_BLOCK(g_b, &taskListB[nextB], b_ar[shiftB], transb,
            bilo, bjlo, &bdim_next, &gNbhdlB[shiftB]);

        adim=adim_next; bdim=bdim_next;
        get_new_B = TRUE;
      }

      /*************************************************************
       * Main Parallel DGEMM Loop.
       *************************************************************/
      while(nextA < max_tasks) {
        currA = nextA;
        currB = taskListA[currA].chunkBId;

        idim = cdim = taskListA[currA].dim[0];
        jdim = taskListB[currB].dim[1];
        kdim = taskListA[currA].dim[1];
        bdim=bdim_next;

        /* if beta=0.0 (i.e.if need_scaling=UNSET), then for first shot,
           we can do put, instead of accumulate */
        if(need_scaling == UNSET) do_put = taskListA[currA].do_put; 

        nextA = ++gTaskId; /* get the next task id */

        if(CYCLIC_DISTR_OPT_FLAG && nextA < max_tasks) 
          nextA = (offset+nextA) % max_tasks;


        /* ---- WAIT till we get the current A & B block ---- */
        a = a_ar[shiftA];
        WAIT_GET_BLOCK(&gNbhdlA[shiftA]);
        if(get_new_B){/*Avoid rereading B if it is same patch as last time*/
          get_new_B = FALSE;
          b = b_ar[shiftB];
          WAIT_GET_BLOCK(&gNbhdlB[shiftB]);
        }

        /* ---- GET the next A & B block ---- */
        if(nextA < max_tasks) {
          GET_BLOCK(g_a, &taskListA[nextA], a_ar[(shiftA+1)%2], transa, 
              ailo, ajlo, &adim_next, &gNbhdlA[(shiftA+1)%2]);

          nextB = taskListA[nextA].chunkBId;
          if(currB != nextB) {
            shiftB=((shiftB+1)%2);

            GET_BLOCK(g_b, &taskListB[nextB], b_ar[shiftB], transb, 
                bilo, bjlo, &bdim_next, &gNbhdlB[shiftB]);
          }
        }
        if(currB != nextB) get_new_B = TRUE;

        /* Do the sequential matrix multiply - i.e.BLAS dgemm */
        GAI_DGEMM(atype, transa, transb, idim, jdim, kdim, alpha, 
            a, adim, b, bdim, c, cdim);

        /* Non-blocking Accumulate Operation. Note: skip wait in 1st loop*/
        i0 = cilo + taskListA[currA].lo[0];
        i1 = cilo + taskListA[currA].hi[0];
        j0 = cjlo + taskListB[currB].lo[1];
        j1 = cjlo + taskListB[currB].hi[1];

        if(currA < max_tasks) {
          if (single_task_flag != SET) {
            switch(atype) {
              case C_FLOAT:
              case C_SCPL:
                clo[0] = i0;
                clo[1] = j0;
                chi[0] = i1;
                chi[1] = j1;
                if(do_put==SET) /* Note:do_put is UNSET, if beta!=0.0*/
                  pnga_put(g_c, clo, chi, (float *)c, &cdim);
                else {
                  pnga_acc(g_c, clo, chi, (float *)c, &cdim, &ONE_CF);
                }
                break;
              default:
                clo[0] = i0;
                clo[1] = j0;
                chi[0] = i1;
                chi[1] = j1;
                if(do_put==SET) /* i.e.beta ==0.0 */
                  pnga_put(g_c, clo, chi, (DoublePrecision*)c, &cdim);
                else {
                  pnga_acc(g_c, clo, chi, (DoublePrecision*)c, &cdim,(DoublePrecision*)&ONE);
                }
                break;
            }
          }
        }

        /* shift next buffer..toggles between 0 and 1: as we use 2 buffers, 
           one for computation and the other for communication (overlap) */
        shiftA = ((shiftA+1)%2); 
        adim = adim_next;
      }
    } while(chunks_left);
  } /* while(has_more_blocks) */

  GA_POP_NAME;
}



static void gai_matmul_irreg(transa, transb, alpha, beta, atype,
			     g_a, ailo, aihi, ajlo, ajhi,
			     g_b, bilo, bihi, bjlo, bjhi,
			     g_c, cilo, cihi, cjlo, cjhi,
			     Ichunk, Kchunk, Jchunk, a_ar,b_ar,c_ar, 
			     need_scaling, irregular) 
     
     Integer g_a, ailo, aihi, ajlo, ajhi;    /* patch of g_a */
     Integer g_b, bilo, bihi, bjlo, bjhi;    /* patch of g_b */
     Integer g_c, cilo, cihi, cjlo, cjhi;    /* patch of g_c */
     Integer Ichunk, Kchunk, Jchunk, atype;
     void    *alpha, *beta;
     char    *transa, *transb;
     DoubleComplex **a_ar, **b_ar, **c_ar;
     short int need_scaling, irregular;
{

#if DEBUG_
  Integer me= pnga_nodeid();
#endif
  Integer nproc=pnga_nnodes();
  Integer get_new_B, i, i0, i1, j0, j1;
  Integer ilo, ihi, idim, jlo, jhi, jdim, klo, khi, kdim, ijk=0;
  Integer n, m, k, adim, bdim=0, cdim;
  Integer idim_prev=0, jdim_prev=0, kdim_prev=0;
  Integer adim_prev=0, bdim_prev=0, cdim_prev=0;
  task_list_t taskListC; 
  short int compute_flag=0, shiftA=0, shiftB=0;
  DoubleComplex ONE, *a, *b, *c;
  SingleComplex ONE_CF; 
  Integer grp_me, a_grp = pnga_get_pgroup(g_a);
  Integer clo[2], chi[2];

  GA_PUSH_NAME("ga_matmul_irreg");
  init_task_list(&taskListC);
  ONE.real =1.; ONE.imag =0.;
  ONE_CF.real =1.; ONE_CF.imag =0.;
#if DEBUG_
  if(me==0) { printf("@@ga_matmul_irreg:m,n,k=%ld %ld %ld\n", aihi-ailo+1,
      bjhi-bjlo+1,ajhi-ajlo+1); fflush(stdout); }
#endif

  m = aihi - ailo +1;
  n = bjhi - bjlo +1;
  k = ajhi - ajlo +1;
  a = a_ar[0];
  b = b_ar[0];
  c = c_ar[0];

  grp_me = pnga_pgroup_nodeid(a_grp);
  clo[0] = cilo; clo[1] = cjlo;
  chi[0] = cihi; chi[1] = cjhi;
  if(!need_scaling) pnga_fill_patch(g_c, clo, chi, beta);

  compute_flag=0;     /* take care of the last chunk */

  for(jlo = 0; jlo < n; jlo += Jchunk){ /* loop thru columns of g_c patch */
    jhi = GA_MIN(n-1, jlo+Jchunk-1);
    jdim= jhi - jlo +1;

    for(klo = 0; klo < k; klo += Kchunk){    /* loop cols of g_a patch */
      khi = GA_MIN(k-1, klo+Kchunk-1);          /* loop rows of g_b patch */
      kdim= khi - klo +1;                                     

      /** Each pass through the outer two loops means we need a
        different patch of B.*/
      get_new_B = TRUE;

      for(ilo = 0; ilo < m; ilo+=Ichunk){ /* loop thru rows of g_c patch */

        if(ijk%nproc == grp_me){

          ihi = GA_MIN(m-1, ilo+Ichunk-1);
          idim= cdim = ihi - ilo +1;


          if (*transa == 'n' || *transa == 'N'){ 
            adim = idim;
            i0= ailo+ilo; i1= ailo+ihi;   
            j0= ajlo+klo; j1= ajlo+khi;
            clo[0] = i0;
            clo[1] = j0;
            chi[0] = i1;
            chi[1] = j1;
            pnga_nbget(g_a, clo, chi, a_ar[shiftA], 
                &idim, &gNbhdlA[shiftA]);
          }else{
            adim = kdim;
            i0= ajlo+klo; i1= ajlo+khi;   
            j0= ailo+ilo; j1= ailo+ihi;
            clo[0] = i0;
            clo[1] = j0;
            chi[0] = i1;
            chi[1] = j1;
            pnga_nbget(g_a, clo, chi, a_ar[shiftA],
                &kdim, &gNbhdlA[shiftA]);
          }

          /* Avoid rereading B if it is same patch as last time. */
          if(get_new_B) { 
            if (*transb == 'n' || *transb == 'N'){ 
              bdim = kdim;
              i0= bilo+klo; i1= bilo+khi;
              j0= bjlo+jlo; j1= bjlo+jhi;
              clo[0] = i0;
              clo[1] = j0;
              chi[0] = i1;
              chi[1] = j1;
              pnga_nbget(g_b, clo, chi, b_ar[shiftB], 
                  &kdim, &gNbhdlB[shiftB]);
            }else{
              bdim = jdim;
              i0= bjlo+jlo; i1= bjlo+jhi;   
              j0= bilo+klo; j1= bilo+khi;
              clo[0] = i0;
              clo[1] = j0;
              chi[0] = i1;
              chi[1] = j1;
              pnga_nbget(g_b, clo, chi, b_ar[shiftB], 
                  &jdim, &gNbhdlB[shiftB]);
            }
          }

          if(compute_flag) { /* compute loop */

            if(atype == C_FLOAT) 
              for(i=0;i<idim_prev*jdim_prev;i++) *(((float*)c)+i)=0;
            else if(atype ==  C_DBL)
              for(i=0;i<idim_prev*jdim_prev;i++) *(((double*)c)+i)=0;
            else if(atype ==  C_SCPL)
              for(i=0;i<idim_prev*jdim_prev;i++) {
                ((SingleComplex*)c)[i].real=0;
                ((SingleComplex*)c)[i].imag=0;
              }
            else for(i=0;i<idim_prev*jdim_prev;i++) {
              c[i].real=0;c[i].imag=0; }

            /* wait till we get the previous block */
            a = a_ar[shiftA^1];
            WAIT_GET_BLOCK(&gNbhdlA[shiftA^1]);
            if(taskListC.chunkBId) {
              b = b_ar[shiftB^1];
              WAIT_GET_BLOCK(&gNbhdlB[shiftB^1]);
            }

            /* Do the sequential matrix multiply - i.e.BLAS dgemm */
            GAI_DGEMM(atype, transa, transb, idim_prev, jdim_prev, 
                kdim_prev, alpha, a, adim_prev, b, bdim_prev, 
                c, cdim_prev);

            i0= cilo + taskListC.lo[0];
            i1= cilo + taskListC.hi[0];
            j0= cjlo + taskListC.lo[1];
            j1= cjlo + taskListC.hi[1];

            if(atype == C_FLOAT || atype == C_SCPL) {
              clo[0] = i0;
              clo[1] = j0;
              chi[0] = i1;
              chi[1] = j1;
              pnga_acc(g_c, clo, chi, (float *)c, &cdim_prev, &ONE_CF);
            } else {
              clo[0] = i0;
              clo[1] = j0;
              chi[0] = i1;
              chi[1] = j1;
              pnga_acc(g_c, clo, chi, (DoublePrecision*)c, &cdim_prev, (DoublePrecision*)&ONE);
            }
          }
          compute_flag=1;

          /* meta-data of current block for next compute loop */
          taskListC.lo[0] = ilo; taskListC.hi[0] = ihi;
          taskListC.lo[1] = jlo; taskListC.hi[1] = jhi;
          taskListC.chunkBId = get_new_B;
          idim_prev = idim;   adim_prev = adim;
          jdim_prev = jdim;   bdim_prev = bdim;
          kdim_prev = kdim;   cdim_prev = cdim;

          /* shift bext buffer */
          shiftA ^= 1;
          if(get_new_B) shiftB ^= 1;

          get_new_B = FALSE; /* Until J or K change again */
        }
        ++ijk;
      }
    }
  }

  /* -------- compute the last chunk --------- */
  if(compute_flag) {
    if(atype == C_FLOAT) 
      for(i=0;i<idim_prev*jdim_prev;i++) *(((float*)c)+i)=0;
    else if(atype ==  C_DBL)
      for(i=0;i<idim_prev*jdim_prev;i++) *(((double*)c)+i)=0;
    else if(atype ==  C_SCPL)
      for(i=0;i<idim_prev*jdim_prev;i++) {
        ((SingleComplex*)c)[i].real=0;
        ((SingleComplex*)c)[i].imag=0;
      }
    else for(i=0;i<idim_prev*jdim_prev;i++) {
      c[i].real=0;c[i].imag=0; }

    /* wait till we get the previous block */
    a = a_ar[shiftA^1];
    WAIT_GET_BLOCK(&gNbhdlA[shiftA^1]);
    if(taskListC.chunkBId) {
      b = b_ar[shiftB^1];
      WAIT_GET_BLOCK(&gNbhdlB[shiftB^1]);
    }

    /* Do the sequential matrix multiply - i.e.BLAS dgemm */
    GAI_DGEMM(atype, transa, transb, idim_prev, jdim_prev, 
        kdim_prev, alpha, a, adim_prev, b, bdim_prev, 
        c, cdim_prev);

    i0= cilo + taskListC.lo[0];
    i1= cilo + taskListC.hi[0];
    j0= cjlo + taskListC.lo[1];
    j1= cjlo + taskListC.hi[1];

    if(atype == C_FLOAT || atype == C_SCPL) {
      clo[0] = i0;
      clo[1] = j0;
      chi[0] = i1;
      chi[1] = j1;
      pnga_acc(g_c, clo, chi, (float *)c, &cdim_prev, &ONE_CF);
    } else {
      clo[0] = i0;
      clo[1] = j0;
      chi[0] = i1;
      chi[1] = j1;
      pnga_acc(g_c, clo, chi, (DoublePrecision*)c, &cdim_prev, (DoublePrecision*)&ONE);
    }
  }
  /* ----------------------------------------- */
  GA_POP_NAME;
}

#if DEBUG_

static void check_result(cond, transa, transb, alpha, beta, atype,
                         g_a, ailo, aihi, ajlo, ajhi,
                         g_b, bilo, bihi, bjlo, bjhi,
                         g_c, cilo, cihi, cjlo, cjhi)
 
     Integer g_a, ailo, aihi, ajlo, ajhi;    /* patch of g_a */
     Integer g_b, bilo, bihi, bjlo, bjhi;    /* patch of g_b */
     Integer g_c, cilo, cihi, cjlo, cjhi;    /* patch of g_c */
     Integer atype, cond;
     void    *alpha, *beta;
     char    *transa, *transb;
{
 
    DoubleComplex *tmpa=NULL, *tmpb=NULL, *tmpc2=NULL;
    static DoubleComplex *tmpc_orig = NULL;
    Integer i,j,m,n,k,adim,bdim,cdim;
    Integer factor=sizeof(DoubleComplex)/GAsizeofM(atype);
    BlasInt m_t, n_t, k_t, adim_t, bdim_t, cdim_t;
    Integer alo[2], ahi[2], blo[2], bhi[2], clo[2], chi[2];
    
    m = aihi - ailo +1;
    n = bjhi - bjlo +1;
    k = ajhi - ajlo +1;
    cdim = m;
 
    if(cond==0) { /* store the original matrix C before matmul starts, as 
		     matrix C is subject to change during pnga_matmul */
       tmpc_orig= (DoubleComplex*)malloc(sizeof(DoubleComplex)*(m*n/factor+1));
       if(tmpc_orig==NULL) pnga_error("check_result: malloc failed", 0);
       
       /* get matrix C */

       clo[0] = cilo;
       clo[1] = cjlo;
       chi[0] = cihi;
       chi[1] = cjhi;
       pnga_get(g_c, clo, chi, tmpc_orig, &m);
    }
    else { /* check for CORRECTNESS */
       
       /* Memory Allocation */
       tmpa = (DoubleComplex*)malloc( sizeof(DoubleComplex)* (m*k/factor+1));
       tmpb = (DoubleComplex*)malloc( sizeof(DoubleComplex)* (k*n/factor+1));
       if(tmpa==NULL || tmpb==NULL) pnga_error("check_result: malloc failed", 0);
       
       switch(atype) {
	  case C_FLOAT:
	     for(i=0; i<m*k; i++) ((float*)tmpa)[i] = -1.0;
	     for(i=0; i<k*n; i++) ((float*)tmpb)[i] = -1.0;
	     break;
	  case C_DBL:
	     for(i=0; i<m*k; i++) ((double*)tmpa)[i] = -1.0;
	     for(i=0; i<k*n; i++) ((double*)tmpb)[i] = -1.0;
	     break;
	  case C_DCPL:
	     for(i=0; i<m*k; i++) {tmpa[i].real=-1.0; tmpa[i].imag=0.0;}
	     for(i=0; i<k*n; i++) {tmpb[i].real=-1.0; tmpb[i].imag=0.0;}
	     break;
	  case C_SCPL:
	     for(i=0; i<m*k; i++) {((SingleComplex*)tmpa)[i].real=-1.0; ((SingleComplex*)tmpa)[i].imag=0.0;}
	     for(i=0; i<k*n; i++) {((SingleComplex*)tmpb)[i].real=-1.0; ((SingleComplex*)tmpb)[i].imag=0.0;}
	     break;
          default: 
            pnga_error("ga_matmul_patch: wrong data type", atype);
       }
       
       /* get matrix A */
       if (*transa == 'n' || *transa == 'N'){
         alo[0] = ailo;
         alo[1] = ajlo;
         ahi[0] = aihi;
         ahi[1] = ajhi;
         adim=m;
         pnga_get(g_a, alo, ahi, tmpa, &m);
       } else {
         alo[0] = ajlo;
         alo[1] = ailo;
         ahi[0] = ajhi;
         ahi[1] = aihi;
         adim=k;
         pnga_get(g_a, alo, ahi, tmpa, &k);
       }

       /* get matrix B */
       if (*transb == 'n' || *transb == 'N'){
         blo[0] = bilo;
         blo[1] = bjlo;
         bhi[0] = bihi;
         bhi[1] = bjhi;
         bdim=k;
         pnga_get(g_b, blo, bhi, tmpb, &k);
       } else { 
         blo[0] = bjlo;
         blo[1] = bilo;
         bhi[0] = bjhi;
         bhi[1] = bihi;
         bdim=n;
         pnga_get(g_b, blo, bhi, tmpb, &n);
       }
       

       m_t=m; n_t=n; k_t=k;
       adim_t=adim; bdim_t=bdim; cdim_t=cdim;
#if (defined(CRAY) || defined(WIN32)) && !NOFORT
       pnga_error("check_result: Serial dgemms not defined", 0L);
#else
       switch(atype) {
	  case C_DBL:
#   if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
	     dgemm_(transa, transb, &m_t, &n_t, &k_t, alpha, tmpa, &adim_t,
		    tmpb, &bdim_t, beta, tmpc_orig, &cdim_t, 1, 1);
#   else
	     dgemm_(transa, 1, transb, 1, &m_t, &n_t, &k_t, alpha, tmpa, &adim_t,
		    tmpb, &bdim_t, beta, tmpc_orig, &cdim_t);
#   endif
	     break;
	  case C_DCPL: 
#   if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
	     zgemm_(transa, transb, &m_t, &n_t, &k_t, (DoubleComplex*)alpha,
		    tmpa, &adim_t, tmpb, &bdim_t, beta, tmpc_orig, &cdim_t, 1, 1);
#   else
	     zgemm_(transa, 1, transb, 1, &m_t, &n_t, &k_t, (DoubleComplex*)alpha,
		    tmpa, &adim_t, tmpb, &bdim_t, beta, tmpc_orig, &cdim_t);
#   endif
	     break;
	  default:
	     pnga_error("check_result: data type not supported here", atype);
       }
#endif
       
       printf("CHK:%c%c : %ld %ld %ld %ld: %ld %ld %ld %ld: %ld %ld %ld %ld\n",
	      *transa, *transb, ailo, aihi, ajlo, ajhi, 
	      bilo, bihi, bjlo, bjhi, cilo, cihi, cjlo, cjhi);
       
       free(tmpa); free(tmpb);
       /* after computing c locally, verify it with the values in g_c */
       tmpc2 = (DoubleComplex*)malloc( sizeof(DoubleComplex)* (m*n/factor+1));
       if(tmpc2==NULL) pnga_error("check_result: malloc failed for tmpc2", 0);
       clo[0] = cilo;
       clo[1] = cjlo;
       chi[0] = cihi;
       chi[1] = cjhi;
       pnga_get(g_c, clo, chi, tmpc2, &m);
       
#define _GA_TOL_ 0.1 /* error tolerance */
       
       switch(atype) {
	  
	  case C_FLOAT:
	    {
	       float abs_value=0.0;
	       for(i=0; i<m*n; i++) {
		  abs_value = ((float*)tmpc_orig)[i] - ((float*)tmpc2)[i];
		  if(abs_value > _GA_TOL_ || abs_value < -(_GA_TOL_)) {
		     printf("Values are = %f : %f\n Alpha=%f Beta=%f\n", 
			    ((float*)tmpc_orig)[i], ((float*)tmpc2)[i], 
			    *((float*)alpha), *((float*)beta));
		     pnga_error("Matmul (type:float) check failed", 0);
		  }
	       }
	    }
	    break;
	    
	  case C_DBL:
	    {
	       double abs_value=0.0;
	       for(i=0; i<m*n; i++) {
		  abs_value = ((double*)tmpc_orig)[i] - ((double*)tmpc2)[i];
		  if(abs_value > _GA_TOL_ || abs_value < -(_GA_TOL_)) {
		     printf("Values are = %lf : %lf\n Alpha=%lf Beta=%lf\n", 
			    ((double*)tmpc_orig)[i+j],	((double*)tmpc2)[i+j],
			    *((double*)alpha),*((double*)beta));
		     pnga_error("Matmul (type:double) check failed", 0);
		  }
	       }
	    }
	    break;
	    
	  case C_DCPL:
	    {
	       DoubleComplex abs_value;
	       for(i=0; i<m*n; i++) {
		  abs_value.real = tmpc_orig[i].real - tmpc2[i].real;
		  abs_value.imag = tmpc_orig[i].imag - tmpc2[i].imag;
		  if(abs_value.real>_GA_TOL_ || abs_value.real<-(_GA_TOL_) ||
		     abs_value.imag>_GA_TOL_ || abs_value.imag<-(_GA_TOL_)) {
		     printf("Values= %lf, %lf : %lf, %lf\n", tmpc_orig[i].real,
			    tmpc_orig[i].imag,tmpc2[i].real,tmpc2[i].imag);
		     pnga_error("Matmul (DoubleComplex) check failed", 0);
		  }
	       }
	    }
	    break;

	  case C_SCPL:
	    {
	       SingleComplex abs_value;
	       for(i=0; i<m*n; i++) {
		  abs_value.real = ((SingleComplex*)tmpc_orig)[i].real - ((SingleComplex*)tmpc2)[i].real;
		  abs_value.imag = ((SingleComplex*)tmpc_orig)[i].imag - ((SingleComplex*)tmpc2)[i].imag;
		  if(abs_value.real>_GA_TOL_ || abs_value.real<-(_GA_TOL_) ||
		     abs_value.imag>_GA_TOL_ || abs_value.imag<-(_GA_TOL_)) {
		     printf("Values= %lf, %lf : %lf, %lf\n", ((SingleComplex*)tmpc_orig)[i].real,
			    ((SingleComplex*)tmpc_orig)[i].imag,((SingleComplex*)tmpc2)[i].real,((SingleComplex*)tmpc2)[i].imag);
		     pnga_error("Matmul (SingleComplex) check failed", 0);
		  }
	       }
	    }
	    break;
	    
	  default:
	     pnga_error("ga_matmul_patch: wrong data type", atype);
       }
       printf("Matrix Multiplication check (m,n,k=%ld %ld %ld)...O.K\n",m,n,k);
       fflush(stdout);
       free(tmpc_orig); free(tmpc2);
       tmpc_orig = NULL;
    }
}

#endif

/******************************************
 * PARALLEL DGEMM
 *     i.e.  C = alpha*A*B + beta*C
 ******************************************/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_matmul = pnga_matmul
#endif
void pnga_matmul(transa, transb, alpha, beta,
	       g_a, ailo, aihi, ajlo, ajhi,
	       g_b, bilo, bihi, bjlo, bjhi,
	       g_c, cilo, cihi, cjlo, cjhi)
     
     Integer g_a, ailo, aihi, ajlo, ajhi;    /* patch of g_a */
     Integer g_b, bilo, bihi, bjlo, bjhi;    /* patch of g_b */
     Integer g_c, cilo, cihi, cjlo, cjhi;    /* patch of g_c */
     void    *alpha, *beta;
     char    *transa, *transb;
{
    DoubleComplex *a=NULL, *b, *c, *a_ar[2], *b_ar[2], *c_ar[2];
    Integer adim1=0, adim2=0, bdim1=0, bdim2=0, cdim1=0, cdim2=0, dims[2];
    Integer atype, btype, ctype, rank, me= pnga_nodeid();
    Integer n, m, k, Ichunk, Kchunk, Jchunk;
    Integer loA[2]={0,0}, hiA[2]={0,0};
    Integer loB[2]={0,0}, hiB[2]={0,0};
    Integer loC[2]={0,0}, hiC[2]={0,0};
    int local_sync_begin,local_sync_end;
    short int need_scaling=SET,use_NB_matmul=SET;
    short int irregular=UNSET, use_armci_memory=UNSET;
    Integer a_grp=pnga_get_pgroup(g_a), b_grp=pnga_get_pgroup(g_b);
    Integer c_grp=pnga_get_pgroup(g_c);
    Integer numblocks;
    Integer clo[2], chi[2];

    /* OPTIMIZATIONS FLAGS. To unset an optimization, replace SET by UNSET) */
    CYCLIC_DISTR_OPT_FLAG  = UNSET;
    CONTIG_CHUNKS_OPT_FLAG = SET;
    DIRECT_ACCESS_OPT_FLAG = SET;

    local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
    _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
    if(local_sync_begin)pnga_pgroup_sync(a_grp);

    GA_PUSH_NAME("pnga_matmul");

    if (a_grp != b_grp || a_grp != c_grp)
       pnga_error("Arrays must be defined on same group",0L);
# if 0 /* disabled. should not fail if there are non-overlapping patches*/
    /* check if C is different from A and B */
    if (g_c == g_a || g_c == g_b)
       pnga_error("Global Array C should be different from A and B", 0);
#endif
    
    /**************************************************
     * Do All Sanity Checks 
     **************************************************/

    /* Check to make sure all global arrays are of the same type */
    if (!(pnga_is_mirrored(g_a) == pnga_is_mirrored(g_b) &&
	  pnga_is_mirrored(g_a) == pnga_is_mirrored(g_c))) {
       pnga_error("Processors do not match for all arrays",pnga_nnodes());
    }

    /* check if ranks are O.K. */
    pnga_inquire(g_a, &atype, &rank, dims); 
    VECTORCHECK(rank, dims, adim1, adim2, ailo, aihi, ajlo, ajhi);
    pnga_inquire(g_b, &btype, &rank, dims); 
    VECTORCHECK(rank, dims, bdim1, bdim2, bilo, bihi, bjlo, bjhi);
    pnga_inquire(g_c, &ctype, &rank, dims); 
    VECTORCHECK(rank, dims, cdim1, cdim2, cilo, cihi, cjlo, cjhi);

    /* check for data-types mismatch */
    if(atype != btype || atype != ctype ) pnga_error(" types mismatch ", 0L);
    if(atype != C_DCPL && atype != C_DBL && atype != C_FLOAT && atype!=C_SCPL)
       pnga_error(" type error",atype);
   
    /* check if patch indices and dims match */
    if (*transa == 'n' || *transa == 'N'){
       if (ailo <= 0 || aihi > adim1 || ajlo <= 0 || ajhi > adim2)
	  pnga_error("  g_a indices out of range ", g_a);
    }else
       if (ailo <= 0 || aihi > adim2 || ajlo <= 0 || ajhi > adim1)
	  pnga_error("  g_a indices out of range ", g_a);
   
    if (*transb == 'n' || *transb == 'N'){
       if (bilo <= 0 || bihi > bdim1 || bjlo <= 0 || bjhi > bdim2)
	  pnga_error("  g_b indices out of range ", g_b);
    }else
       if (bilo <= 0 || bihi > bdim2 || bjlo <= 0 || bjhi > bdim1)
	  pnga_error("  g_b indices out of range ", g_b);
   
    if (cilo <= 0 || cihi > cdim1 || cjlo <= 0 || cjhi > cdim2)
       pnga_error("  g_c indices out of range ", g_c);

    /* verify if patch dimensions are consistent */
    m = aihi - ailo +1;
    n = bjhi - bjlo +1;
    k = ajhi - ajlo +1;
    if( (cihi - cilo +1) != m) pnga_error(" a & c dims error",m);
    if( (cjhi - cjlo +1) != n) pnga_error(" b & c dims error",n);
    if( (bihi - bilo +1) != k) pnga_error(" a & b dims error",k);

#if DEBUG_
    if(me==0) check_result(0, transa, transb, alpha, beta, atype,
			   g_a, ailo, aihi, ajlo, ajhi,
			   g_b, bilo, bihi, bjlo, bjhi,
			   g_c, cilo, cihi, cjlo, cjhi);
    pnga_sync();
#endif

    /* switch to various matmul algorithms here. more to come */
    if( GA[GA_OFFSET + g_c].irreg == 1 ||
	GA[GA_OFFSET + g_b].irreg == 1 ||
	GA[GA_OFFSET + g_a].irreg == 1 ||
	_gai_matmul_patch_flag == SET) irregular = SET;

    /* even ga_dgemm is called, m,n & k might not match GA dimensions */
    pnga_inquire(g_c, &ctype, &rank, dims);
    if(dims[0] != m || dims[1] != n) irregular = SET; /* C matrix dims */

    if(!irregular) {
       if((adim1=GA_Cluster_nnodes()) > 1) use_NB_matmul = SET;
       else {
	  use_NB_matmul = UNSET;
	  CONTIG_CHUNKS_OPT_FLAG = UNSET;
	  DIRECT_ACCESS_OPT_FLAG = UNSET;
       }
#    if defined(__crayx1) || defined(NEC)
       use_NB_matmul = UNSET;
#    endif
    }

    /* if block cyclic, then use regular algorithm. This is turned on for now
     * to test block cyclic */ 
    numblocks = pnga_total_blocks(g_c);
    if(numblocks>=0) {
       irregular     = UNSET;
       use_NB_matmul = SET; 
    }
    
    /****************************************************************
     * Get the memory (i.e.static or dynamic) for temporary buffers 
     ****************************************************************/

    /* to skip accumulate and exploit data locality:
       get chunks according to "C" matrix distribution*/
    pnga_distribution(g_a, me, loA, hiA);
    pnga_distribution(g_b, me, loB, hiB);
    pnga_distribution(g_c, me, loC, hiC);

       {
	  Integer elems, factor=sizeof(DoubleComplex)/GAsizeofM(atype);
	  short int nbuf=1;
	  DoubleComplex *tmp = NULL;

	  Ichunk = GA_MIN( (hiC[0]-loC[0]+1), (hiA[0]-loA[0]+1) );
	  Jchunk = GA_MIN( (hiC[1]-loC[1]+1), (hiB[1]-loB[1]+1) );
	  Kchunk = GA_MIN( (hiA[1]-loA[1]+1), (hiB[0]-loB[0]+1) );

#if KCHUNK_OPTIMIZATION /*works great for m=1000,n=1000,k=4000 kinda cases*/
	  pnga_distribution(g_a, me, loC, hiC);
	  Kchunk = hiC[1]-loC[1]+1;
	  pnga_distribution(g_b, me, loC, hiC);
	  Kchunk = GA_MIN(Kchunk, (hiC[0]-loC[0]+1));
#endif

	  /* Just to avoid divide by zero error */
          if(Ichunk<=0) Ichunk = 1;
          if(Jchunk<=0) Jchunk = 1;
          if(Kchunk<=0) Kchunk = 1;

	  {
	     Integer irreg=0;
	     if(Ichunk/Kchunk > GA_ASPECT_RATIO || Kchunk/Ichunk > GA_ASPECT_RATIO || 
		Jchunk/Kchunk > GA_ASPECT_RATIO || Kchunk/Jchunk > GA_ASPECT_RATIO) {
                irreg = SET;
             }
	     pnga_pgroup_gop(a_grp, pnga_type_f2c(MT_F_INT), &irreg, (Integer)1, "max");   
	     if(irreg==SET) irregular = SET;
	  }
	  
	  /* If non-blocking, we need 2 temporary buffers for A and B matrix */
	  if(use_NB_matmul) nbuf = 2; 
	  
	  if(!irregular) {
	     tmp = a_ar[0] =a=gai_get_armci_memory(Ichunk,Jchunk,Kchunk,
						   nbuf, atype);
	     if(tmp != NULL) use_armci_memory = SET;
	  }
	  
	  /* get ChunkSize (i.e.BlockSize), that fits in temporary buffer */
	  gai_get_chunk_size(irregular, &Ichunk, &Jchunk, &Kchunk, &elems, 
			     atype, m, n, k, nbuf, use_armci_memory, a_grp);
	  
	  if(tmp == NULL) { /* try once again from armci for new chunk sizes */
	     tmp = a_ar[0] =a=gai_get_armci_memory(Ichunk,Jchunk,Kchunk,
						   nbuf, atype);
	     if(tmp != NULL) use_armci_memory = SET;
	  }

	  if(tmp == NULL) { /*if armci malloc fails again, then get from MA */
	     tmp = a_ar[0] = a =(DoubleComplex*) ga_malloc(elems,atype,
							   "GA mulmat bufs");
	  }

	  if(use_NB_matmul) tmp = a_ar[1] = a_ar[0] + (Ichunk*Kchunk)/factor+1;
	  
	  tmp = b_ar[0] = b = tmp + (Ichunk*Kchunk)/factor + 1;
	  if(use_NB_matmul) tmp = b_ar[1] = b_ar[0] + (Kchunk*Jchunk)/factor+1;
	  
	  c_ar[0] = c = tmp + (Kchunk*Jchunk)/factor + 1;
       }
       
       /** check if there is a need for scaling the data. 
	   Note: if beta=0, then need_scaling=0  */
       if(atype==C_DCPL){
	  if((((DoubleComplex*)beta)->real == 0) && 
	     (((DoubleComplex*)beta)->imag ==0)) need_scaling =0;} 
       else if(atype==C_SCPL){
	  if((((SingleComplex*)beta)->real == 0) && 
	     (((SingleComplex*)beta)->imag ==0)) need_scaling =0;} 
       else if((atype==C_DBL)){
	  if(*(DoublePrecision *)beta == 0) need_scaling =0;}
       else if( *(float*)beta ==0) need_scaling =0;

       clo[0] = cilo; clo[1] = cjlo;
       chi[0] = cihi; chi[1] = cjhi;
       if(need_scaling) pnga_scale_patch(g_c, clo, chi, beta);

       /********************************************************************
	* Parallel Matrix Multiplication Starts Here.
	* 3 Steps:
	*    1. Get a chunk of A and B matrix, and store it in local buffer.
	*    2. Do sequential dgemm.
	*    3. Put/accumulate the result into C matrix.
	*********************************************************************/

       /* if only one node, then enable the optimized shmem code */
       if(use_NB_matmul==UNSET) { 
	  gai_matmul_shmem(transa, transb, alpha, beta, atype,
			   g_a, ailo, aihi, ajlo, ajhi,
			   g_b, bilo, bihi, bjlo, bjhi,
			   g_c, cilo, cihi, cjlo, cjhi,
			   Ichunk, Kchunk, Jchunk, a,b,c, need_scaling);
       }
       else {
	  if(irregular)
	     gai_matmul_irreg(transa, transb, alpha, beta, atype,
			      g_a, ailo, aihi, ajlo, ajhi,
			      g_b, bilo, bihi, bjlo, bjhi,
			      g_c, cilo, cihi, cjlo, cjhi,
			      Ichunk, Kchunk, Jchunk, a_ar, b_ar, c_ar,
			      need_scaling, irregular);
	  else
	     gai_matmul_regular(transa, transb, alpha, beta, atype,
				g_a, ailo, aihi, ajlo, ajhi,
				g_b, bilo, bihi, bjlo, bjhi,
				g_c, cilo, cihi, cjlo, cjhi,
				Ichunk, Kchunk, Jchunk, a_ar, b_ar, c_ar, 
				need_scaling, irregular);
       }
	     
       a = a_ar[0];
       if(use_armci_memory == SET) ARMCI_Free_local(a);
       else ga_free(a);
       
#if DEBUG_
       Integer grp_me;
       grp_me = pnga_pgroup_nodeid(a_grp);
       pnga_pgroup_sync(a_grp);
       if(me==0) check_result(1, transa, transb, alpha, beta, atype,
			      g_a, ailo, aihi, ajlo, ajhi,
			      g_b, bilo, bihi, bjlo, bjhi,
			      g_c, cilo, cihi, cjlo, cjhi);
       pnga_pgroup_sync(a_grp);
#endif
       
       GA_POP_NAME;   
       if(local_sync_end)pnga_pgroup_sync(a_grp);
}

/* This is the old matmul code. It is enabled now for mirrored matrix multiply. 
   It also work for normal matrix/vector multiply with no changes */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_matmul_mirrored = pnga_matmul_mirrored
#endif
void pnga_matmul_mirrored(transa, transb, alpha, beta,
			g_a, ailo, aihi, ajlo, ajhi,
			g_b, bilo, bihi, bjlo, bjhi,
			g_c, cilo, cihi, cjlo, cjhi)

     Integer g_a, ailo, aihi, ajlo, ajhi;    /* patch of g_a */
     Integer g_b, bilo, bihi, bjlo, bjhi;    /* patch of g_b */
     Integer g_c, cilo, cihi, cjlo, cjhi;    /* patch of g_c */
     void    *alpha, *beta;
     char    *transa, *transb;
{

    #ifdef STATBUF
  /* approx. sqrt(2) ratio in chunk size to use the same buffer space */
   DoubleComplex a[ICHUNK*KCHUNK], b[KCHUNK*JCHUNK], c[ICHUNK*JCHUNK];
#else
   DoubleComplex *a, *b, *c;
#endif
Integer atype, btype, ctype, adim1=0, adim2=0, bdim1=0, bdim2=0, cdim1=0, cdim2=0, dims[2], rank;
Integer me= pnga_nodeid(), nproc;
Integer i, ijk = 0, i0, i1, j0, j1;
Integer ilo, ihi, idim, jlo, jhi, jdim, klo, khi, kdim;
Integer n, m, k, adim, bdim=0, cdim;
Integer Ichunk, Kchunk, Jchunk;
DoubleComplex ONE;
SingleComplex ONE_CF;

DoublePrecision chunk_cube;
Integer min_tasks = 10, max_chunk;
int need_scaling=1;
Integer ZERO_I = 0, inode, iproc;
Integer get_new_B;
int local_sync_begin,local_sync_end;
BlasInt idim_t, jdim_t, kdim_t, adim_t, bdim_t, cdim_t;
Integer clo[2], chi[2];

   ONE.real =1.; ONE.imag =0.;
   ONE_CF.real =1.; ONE_CF.imag =0.;

   local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
   _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
   if(local_sync_begin)pnga_sync();

   GA_PUSH_NAME("ga_matmul_patch");

   /* Check to make sure all global arrays are of the same type */
   if (!(pnga_is_mirrored(g_a) == pnga_is_mirrored(g_b) &&
        pnga_is_mirrored(g_a) == pnga_is_mirrored(g_c))) {
     pnga_error("Processors do not match for all arrays",pnga_nnodes());
   }
   if (pnga_is_mirrored(g_a)) {
     inode = pnga_cluster_nodeid();
     nproc = pnga_cluster_nprocs(inode);
     iproc = me - pnga_cluster_procid(inode, ZERO_I);
   } else {
     nproc = pnga_nnodes();
     iproc = me;
   }

   pnga_inquire(g_a, &atype, &rank, dims); 
   VECTORCHECK(rank, dims, adim1, adim2, ailo, aihi, ajlo, ajhi);
   pnga_inquire(g_b, &btype, &rank, dims); 
   VECTORCHECK(rank, dims, bdim1, bdim2, bilo, bihi, bjlo, bjhi);
   pnga_inquire(g_c, &ctype, &rank, dims); 
   VECTORCHECK(rank, dims, cdim1, cdim2, cilo, cihi, cjlo, cjhi);

   if(atype != btype || atype != ctype ) pnga_error(" types mismatch ", 0L);
   if(atype != C_DCPL && atype != C_DBL && atype != C_FLOAT && atype != C_SCPL)
     pnga_error(" type error",atype);
   
   
   
   /* check if patch indices and dims match */
   if (*transa == 'n' || *transa == 'N'){
     if (ailo <= 0 || aihi > adim1 || ajlo <= 0 || ajhi > adim2)
       pnga_error("  g_a indices out of range ", g_a);
   }else
     if (ailo <= 0 || aihi > adim2 || ajlo <= 0 || ajhi > adim1)
       pnga_error("  g_a indices out of range ", g_a);
   
   if (*transb == 'n' || *transb == 'N'){
     if (bilo <= 0 || bihi > bdim1 || bjlo <= 0 || bjhi > bdim2)
       pnga_error("  g_b indices out of range ", g_b);
   }else
     if (bilo <= 0 || bihi > bdim2 || bjlo <= 0 || bjhi > bdim1)
       pnga_error("  g_b indices out of range ", g_b);
   
   if (cilo <= 0 || cihi > cdim1 || cjlo <= 0 || cjhi > cdim2)
     pnga_error("  g_c indices out of range ", g_c);
   
   /* verify if patch dimensions are consistent */
   m = aihi - ailo +1;
   n = bjhi - bjlo +1;
   k = ajhi - ajlo +1;
   if( (cihi - cilo +1) != m) pnga_error(" a & c dims error",m);
   if( (cjhi - cjlo +1) != n) pnga_error(" b & c dims error",n);
   if( (bihi - bilo +1) != k) pnga_error(" a & b dims error",k);
   

   /* In 32-bit platforms, k*m*n might exceed the "long" range(2^31), 
      eg:k=m=n=1600. So casting the temporary value to "double" helps */
   chunk_cube = (k*(double)(m*n)) / (min_tasks * nproc);
   max_chunk = (Integer)pow(chunk_cube, (DoublePrecision)(1.0/3.0) );
   if (max_chunk < 32) max_chunk = 32;

#ifdef STATBUF
   if(atype ==  C_DBL || atype == C_FLOAT){
      Ichunk=D_CHUNK, Kchunk=D_CHUNK, Jchunk=D_CHUNK;
   }else{
      Ichunk=ICHUNK; Kchunk=KCHUNK; Jchunk=JCHUNK;
   }
#else
   {
     /**
      * Find out how much memory we can grab.  It will be used in
      * three chunks, and the result includes only the first one.
      */
     
     Integer elems, factor = sizeof(DoubleComplex)/GAsizeofM(atype);
     Ichunk = Jchunk = Kchunk = CHUNK_SIZE;
     
     if ( max_chunk > Ichunk) {       
       /*if memory if very limited, performance degrades for large matrices
	 as chunk size is very small, which leads to communication overhead)*/
       Integer avail = pnga_memory_avail_type(atype);
       if (pnga_is_mirrored(g_a)) {
         fflush(stdout);
         if (sizeof(Integer)/sizeof(int) > 1)
           armci_msg_gop_scope(SCOPE_NODE, &avail, 1, "min", ARMCI_LONG);
         else
           armci_msg_gop_scope(SCOPE_NODE, &avail, 1, "min", ARMCI_INT);
         fflush(stdout);
       } else {
         fflush(stdout);
         pnga_gop(pnga_type_f2c(MT_F_INT), &avail, (Integer)1, "min");
       }
       if(avail<MINMEM && pnga_nodeid()==0) pnga_error("NotEnough memory",avail);
       elems = (Integer)(avail*0.9); /* Donot use every last drop */
       
       max_chunk=GA_MIN(max_chunk, (Integer)(sqrt( (double)((elems-EXTRA)/3))));
       Ichunk = GA_MIN(m,max_chunk);
       Jchunk = GA_MIN(n,max_chunk);
       Kchunk = GA_MIN(k,max_chunk);
     }
     else /* "EXTRA" elems for safety - just in case */
       elems = 3*Ichunk*Jchunk + EXTRA*factor;
     
     a = (DoubleComplex*) ga_malloc(elems, atype, "GA mulmat bufs");
     b = a + (Ichunk*Kchunk)/factor + 1; 
     c = b + (Kchunk*Jchunk)/factor + 1;
   }
#endif

   if(atype==C_DCPL){if((((DoubleComplex*)beta)->real == 0) &&
	       (((DoubleComplex*)beta)->imag ==0)) need_scaling =0;} 
   else if(atype==C_SCPL){if((((SingleComplex*)beta)->real == 0) &&
	       (((SingleComplex*)beta)->imag ==0)) need_scaling =0;} 
   else if((atype==C_DBL)){if(*(DoublePrecision *)beta == 0) need_scaling =0;}
   else if( *(float*)beta ==0) need_scaling =0;

   pnga_mask_sync(ZERO_I, ZERO_I);
   clo[0] = cilo; clo[1] = cjlo;
   chi[0] = cihi; chi[1] = cjhi;
   if(need_scaling) pnga_scale_patch(g_c, clo, chi, beta);
   else  pnga_fill_patch(g_c, clo, chi, beta);

   for(jlo = 0; jlo < n; jlo += Jchunk){ /* loop through columns of g_c patch */
       jhi = GA_MIN(n-1, jlo+Jchunk-1);
       jdim= jhi - jlo +1;

       for(klo = 0; klo < k; klo += Kchunk){    /* loop cols of g_a patch */
	 khi = GA_MIN(k-1, klo+Kchunk-1);          /* loop rows of g_b patch */
	 kdim= khi - klo +1;                                     
	 
	 /** Each pass through the outer two loops means we need a
	     different patch of B.*/
	 get_new_B = TRUE;
	 
	 for(ilo = 0; ilo < m; ilo += Ichunk){ /*loop through rows of g_c patch */
	   
	   if(ijk%nproc == iproc){

	     ihi = GA_MIN(m-1, ilo+Ichunk-1);
	     idim= cdim = ihi - ilo +1;
	     
	     if(atype == C_FLOAT) 
	       for (i = 0; i < idim*jdim; i++) *(((float*)c)+i)=0;
	     else if(atype ==  C_DBL)
	       for (i = 0; i < idim*jdim; i++) *(((double*)c)+i)=0;
	     else if(atype ==  C_SCPL)
	       for (i = 0; i < idim*jdim; i++){
                 ((SingleComplex*)c)[i].real=0;
                 ((SingleComplex*)c)[i].imag=0;
               }
	     else
	       for (i = 0; i < idim*jdim; i++){ c[i].real=0;c[i].imag=0;}
	     
	     if (*transa == 'n' || *transa == 'N'){ 
	       adim = idim;
	       i0= ailo+ilo; i1= ailo+ihi;   
	       j0= ajlo+klo; j1= ajlo+khi;
          clo[0] = i0;
          clo[1] = j0;
          chi[0] = i1;
          chi[1] = j1;
	       pnga_get(g_a, clo, chi, a, &idim);
	     }else{
	       adim = kdim;
	       i0= ajlo+klo; i1= ajlo+khi;   
	       j0= ailo+ilo; j1= ailo+ihi;
          clo[0] = i0;
          clo[1] = j0;
          chi[0] = i1;
          chi[1] = j1;
	       pnga_get(g_a, clo, chi, a, &kdim);
	     }


	     /* Avoid rereading B if it is the same patch as last time. */
        if(get_new_B) { 
          if (*transb == 'n' || *transb == 'N'){ 
            bdim = kdim;
            i0= bilo+klo; i1= bilo+khi;   
            j0= bjlo+jlo; j1= bjlo+jhi;
            clo[0] = i0;
            clo[1] = j0;
            chi[0] = i1;
            chi[1] = j1;
            pnga_get(g_b, clo, chi, b, &kdim);
          }else{
            bdim = jdim;
            i0= bjlo+jlo; i1= bjlo+jhi;   
            j0= bilo+klo; j1= bilo+khi;
            clo[0] = i0;
            clo[1] = j0;
            chi[0] = i1;
            chi[1] = j1;
            pnga_get(g_b, clo, chi, b, &jdim);
          }
          get_new_B = FALSE; /* Until J or K change again */
        }

	     
	     idim_t=idim; jdim_t=jdim; kdim_t=kdim;
	     adim_t=adim; bdim_t=bdim; cdim_t=cdim;

	     switch(atype) {
	     case C_FLOAT:
	       BLAS_SGEMM(transa, transb, &idim_t, &jdim_t, &kdim_t,
		      (Real *)alpha, (Real *)a, &adim_t,
              (Real *)b, &bdim_t, (Real *)&ONE_CF,
              (Real *)c, &cdim_t);
	       break;
	     case C_DBL:
	       BLAS_DGEMM(transa, transb, &idim_t, &jdim_t, &kdim_t,
		      (DoublePrecision *)alpha, (DoublePrecision *)a, &adim_t,
              (DoublePrecision *)b, &bdim_t, (DoublePrecision *)&ONE,
              (DoublePrecision *)c, &cdim_t);
	       break;
	     case C_DCPL:
	       BLAS_ZGEMM(transa, transb, &idim_t, &jdim_t, &kdim_t,
		      (DoubleComplex *)alpha, (DoubleComplex *)a, &adim_t,
              (DoubleComplex *)b, &bdim_t, (DoubleComplex *)&ONE,
              (DoubleComplex *)c, &cdim_t);
	       break;
	     case C_SCPL:
	       BLAS_CGEMM(transa, transb, &idim_t, &jdim_t, &kdim_t,
		      (SingleComplex *)alpha, (SingleComplex *)a, &adim_t,
              (SingleComplex *)b, &bdim_t, (SingleComplex *)&ONE_CF,
              (SingleComplex *)c, &cdim_t);
	       break;
	     default:
	       pnga_error("ga_matmul_patch: wrong data type", atype);
	     }
	     
	     i0= cilo+ilo; i1= cilo+ihi;   j0= cjlo+jlo; j1= cjlo+jhi;
	     if(atype == C_FLOAT || atype == C_SCPL)  {
          clo[0] = i0;
          clo[1] = j0;
          chi[0] = i1;
          chi[1] = j1;
	       pnga_acc(g_c, clo, chi, (float *)c, &cdim, &ONE_CF);
        } else {
          clo[0] = i0;
          clo[1] = j0;
          chi[0] = i1;
          chi[1] = j1;
	       pnga_acc(g_c, clo, chi, (DoublePrecision*)c, &cdim, (DoublePrecision*)&ONE);
        }
	   }
	   ++ijk;
	 }
       }
   }
   
#ifndef STATBUF
   ga_free(a);
#endif

   GA_POP_NAME;
   if(local_sync_end)pnga_sync();

}


#if 0
void gai_matmul_patch(char *transa, char *transb, void *alpha, void *beta,
        Integer g_a,Integer ailo,Integer aihi,Integer ajlo,Integer ajhi,
        Integer g_b,Integer bilo,Integer bihi,Integer bjlo,Integer bjhi,
        Integer g_c,Integer cilo,Integer cihi,Integer cjlo,Integer cjhi)
{
    if(pnga_is_mirrored(g_a)) 
       pnga_matmul_mirrored(transa, transb, alpha, beta,
			  g_a, ailo, aihi, ajlo, ajhi,
			  g_b, bilo, bihi, bjlo, bjhi,
			  g_c, cilo, cihi, cjlo, cjhi);
    else {
       _gai_matmul_patch_flag = SET;
       pnga_matmul(transa, transb, alpha, beta,
		 g_a, ailo, aihi, ajlo, ajhi,
		 g_b, bilo, bihi, bjlo, bjhi,
		 g_c, cilo, cihi, cjlo, cjhi);
       _gai_matmul_patch_flag = UNSET;
    }
}
#endif


/*\ select the 2d plane to be used in matrix multiplication                     
  \*/
  static void  gai_setup_2d_patch(Integer rank, char *trans, Integer dims[],
                                Integer lo[], Integer hi[],
                                Integer* ilo, Integer* ihi,
                                Integer* jlo, Integer* jhi,
                                Integer* dim1, Integer* dim2,
                                int* ipos, int* jpos, int *vpos)
{
    int d,e=0;
    char t='n';
    
    for(d=0; d<rank; d++)
       if( (hi[d]-lo[d])>0 && ++e>2 ) pnga_error("3-D Patch Detected", 0L);
    *ipos = *jpos = *vpos = -1;
    for(d=0; d<rank; d++){
       if( (*ipos <0) && (hi[d]>lo[d]) ) { *ipos =d; continue; }
       if( (*ipos >=0) && (hi[d]>lo[d])) { *jpos =d; break; }
    }
    /* we have an ambiguous vector; mark its location */
    if (1 == e && *ipos != 0 && *ipos != rank-1) {
        *vpos = *ipos;
    }

    /*    if(*ipos >*jpos){Integer t=*ipos; *ipos=*jpos; *jpos=t;} 
     */
    
    /* single element case (trivial) */
    if((*ipos <0) && (*jpos <0)){ *ipos =0; *jpos=1; }
    else{
       
       /* handle almost trivial case of only one dimension with >1 elements */
       if(trans == NULL) trans = &t;      
       if(*trans == 'n' || *trans == 'N') {
          if(*ipos == rank-1) (*ipos)--; /* i cannot be the last dimension */
          if(*ipos <0) *ipos = *jpos-1; /* select i dimension based on j */
          if(*jpos <0) *jpos = *ipos+1; /* select j dimenison based on i */
       }
       else {
          if(*ipos <0) *ipos = *jpos-1; 
          if(*jpos <0) {
             if(*ipos==0) *jpos = *ipos + 1; 
             else         *jpos = (*ipos)--;
          }
       }
    }
    
    *ilo = lo[*ipos]; *ihi = hi[*ipos];
    *jlo = lo[*jpos]; *jhi = hi[*jpos];
    *dim1 = dims[*ipos];
    *dim2 = dims[*jpos];
#if 0
    printf("lo/hi=[%ld:%ld", lo[0], hi[0]);
    for (d=1; d<rank; ++d) {
        printf(",%ld:%ld", lo[d], hi[d]);
    }
    printf("]\n");
    printf("size=[%ld", hi[0]-lo[0]+1);
    for (d=1; d<rank; ++d) {
        printf(",%ld", hi[d]-lo[d]+1);
    }
    printf("]\n");
    printf("gai_setup_2d_patch(%ld, T, X, X, X, %ld, %ld, %ld, %ld, %ld, %ld, %d, %d, %d)\n",
            rank, *ilo, *ihi, *jlo, *jhi, *dim1, *dim2, *ipos, *jpos, *vpos);
#endif
}

#define  SETINT(tmp,val,n) {int _i; for(_i=0;_i<n; _i++)tmp[_i]=val;}

/*\ MATRIX MULTIPLICATION for 2d patches of multi-dimensional arrays 
 *  
 *  C[lo:hi,lo:hi] = alpha*op(A)[lo:hi,lo:hi] * op(B)[lo:hi,lo:hi]        
 *                 + beta *C[lo:hi,lo:hi]
 *
 *  where:
 *          op(A) = A or A' depending on the transpose flag
 *  [lo:hi,lo:hi] - patch indices _after_ op() operator was applied
 *
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_matmul_patch = pnga_matmul_patch
#endif
void pnga_matmul_patch(char *transa, char *transb, void *alpha, void *beta, 
		      Integer g_a, Integer alo[], Integer ahi[], 
                      Integer g_b, Integer blo[], Integer bhi[], 
		      Integer g_c, Integer clo[], Integer chi[])
{
#ifdef STATBUF
   DoubleComplex a[ICHUNK*KCHUNK], b[KCHUNK*JCHUNK], c[ICHUNK*JCHUNK];
#else
   DoubleComplex *a, *b, *c;
#endif
Integer atype, btype, ctype, adim1, adim2, bdim1, bdim2, cdim1, cdim2;
Integer me= pnga_nodeid(), nproc, inode, iproc;
Integer i, ijk = 0, i0, i1, j0, j1;
Integer ilo, ihi, idim, jlo, jhi, jdim, klo, khi, kdim;
Integer n, m, k, k2, cm, cn, adim, bdim, cdim, arank, brank, crank;
int aipos, ajpos, bipos, bjpos,cipos, cjpos, need_scaling=1;
int avpos, bvpos, cvpos;
Integer Ichunk, Kchunk, Jchunk;
Integer ailo, aihi, ajlo, ajhi;    /* 2d plane of g_a */
Integer bilo, bihi, bjlo, bjhi;    /* 2d plane of g_b */
Integer cilo, cihi, cjlo, cjhi;    /* 2d plane of g_c */
Integer adims[GA_MAX_DIM],bdims[GA_MAX_DIM],cdims[GA_MAX_DIM],tmpld[GA_MAX_DIM];
Integer *tmplo = adims, *tmphi =bdims; 
DoubleComplex ONE;
SingleComplex ONE_CF;
Integer ZERO_I = 0;
Integer get_new_B;
DoublePrecision chunk_cube;
Integer min_tasks = 10, max_chunk;
int local_sync_begin,local_sync_end;
BlasInt idim_t, jdim_t, kdim_t, adim_t, bdim_t, cdim_t;

   ONE.real =1.; ONE.imag =0.;
   ONE_CF.real =1.; ONE_CF.imag =0.;
   
   local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
   _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
   if(local_sync_begin)pnga_sync();

   GA_PUSH_NAME("nga_matmul_patch");

   /* Check to make sure all global arrays are of the same type */
   if (!(pnga_is_mirrored(g_a) == pnga_is_mirrored(g_b) &&
        pnga_is_mirrored(g_a) == pnga_is_mirrored(g_c))) {
     pnga_error("Processors do not match for all arrays",pnga_nnodes());
   }
   if (pnga_is_mirrored(g_a)) {
     inode = pnga_cluster_nodeid();
     nproc = pnga_cluster_nprocs(inode);
     iproc = me - pnga_cluster_procid(inode, ZERO_I);
   } else {
     nproc = pnga_nnodes();
     iproc = me;
   }

   pnga_inquire(g_a, &atype, &arank, adims);
   pnga_inquire(g_b, &btype, &brank, bdims);
   pnga_inquire(g_c, &ctype, &crank, cdims);

   if(arank<2)  pnga_error("rank of A must be at least 2",arank);
   if(brank<2)  pnga_error("rank of B must be at least 2",brank);
   if(crank<2)  pnga_error("rank of C must be at least 2",crank);

   if(atype != btype || atype != ctype ) pnga_error(" types mismatch ", 0L);
   if(atype != C_DCPL && atype != C_DBL && atype != C_FLOAT && atype != C_SCPL)
     pnga_error(" type error",atype);
   
   gai_setup_2d_patch(arank, transa, adims, alo, ahi, &ailo, &aihi,
                      &ajlo, &ajhi, &adim1, &adim2, &aipos, &ajpos, &avpos);
   gai_setup_2d_patch(brank, transb, bdims, blo, bhi, &bilo, &bihi,
                      &bjlo, &bjhi, &bdim1, &bdim2, &bipos, &bjpos, &bvpos);
   gai_setup_2d_patch(crank, NULL, cdims, clo, chi, &cilo, &cihi,
                      &cjlo, &cjhi, &cdim1, &cdim2, &cipos, &cjpos, &cvpos);

   /* check if patch indices and dims match */
   if (*transa == 'n' || *transa == 'N'){
      if (ailo <= 0 || aihi > adim1 || ajlo <= 0 || ajhi > adim2)
         pnga_error("  g_a indices out of range ", g_a);
   }else
      if (ailo <= 0 || aihi > adim2 || ajlo <= 0 || ajhi > adim1)
         pnga_error("  g_a indices out of range ", g_a);

   if (*transb == 'n' || *transb == 'N'){
      if (bilo <= 0 || bihi > bdim1 || bjlo <= 0 || bjhi > bdim2)
          pnga_error("  g_b indices out of range ", g_b);
   }else
      if (bilo <= 0 || bihi > bdim2 || bjlo <= 0 || bjhi > bdim1)
          pnga_error("  g_b indices out of range ", g_b);

   if (cilo <= 0 || cihi > cdim1 || cjlo <= 0 || cjhi > cdim2)
       pnga_error("  g_c indices out of range ", g_c);

   /* verify if patch dimensions are consistent */
#define RESET() do {   \
   m = aihi - ailo +1; \
   k = ajhi - ajlo +1; \
   k2= bihi - bilo +1; \
   n = bjhi - bjlo +1; \
   cm= cihi - cilo +1; \
   cn= cjhi - cjlo +1; \
} while (0)
   RESET();
#define SHIFT(L,INC) do {      \
   L##ipos+=INC;               \
   L##jpos+=INC;               \
   L##ilo = L##lo[L##ipos];    \
   L##ihi = L##hi[L##ipos];    \
   L##jlo = L##lo[L##jpos];    \
   L##jhi = L##hi[L##jpos];    \
   L##dim1 = L##dims[L##ipos]; \
   L##dim2 = L##dims[L##jpos]; \
   RESET();                    \
} while (0)
   /* gai_setup_2d_patch may produce ambiguous vectors */
   if (!(m==cm && k==k2 && n==cn)) {
       /* patches don't agree */
       if (avpos>=0 && bvpos<0 && cvpos<0) {
           /* only A is an ambiguous vector */
           SHIFT(a,-1);
       }
       else if (avpos<0 && bvpos>=0 && cvpos<0) {
           /* only B is an ambiguous vector */
           SHIFT(b,-1);
       }
       else if (avpos<0 && bvpos<0 && cvpos>=0) {
           /* only C is an ambiguous vector */
           SHIFT(c,-1);
       }
       else if (avpos>=0 && bvpos>=0 && cvpos<0) {
           /* A and B are ambiguous vectors */
           if (m != cm) {
               SHIFT(a,-1);
           }
           if (n != cn) {
               SHIFT(b,-1);
           }
       }
       else if (avpos>=0 && bvpos<0 && cvpos>=0) {
           /* A and C are ambiguous vectors */
           if (k != k2) {
               SHIFT(a,-1);
           }
           if (n != cn) {
               SHIFT(c,-1);
           }
       }
       else if (avpos<0 && bvpos>=0 && cvpos>=0) {
           /* B and C are ambiguous vectors */
           if (k != k2) {
               SHIFT(b,-1);
           }
           if (m != cm) {
               SHIFT(c,-1);
           }
       }
       else if (avpos>=0 && bvpos>=0 && cvpos>=0) {
           /* A and B and C are ambiguous vectors */
           pnga_error("a and b and c ambiguous", 1);
       }
   }
   if( (cihi - cilo +1) != m) pnga_error(" a & c dims error",m);
   if( (cjhi - cjlo +1) != n) pnga_error(" b & c dims error",n);
   if( (bihi - bilo +1) != k) pnga_error(" a & b dims error",k);
   
   chunk_cube = (k*(double)(m*n)) / (min_tasks * nproc);
   max_chunk = (Integer)pow(chunk_cube, (DoublePrecision)(1.0/3.0) );
   if (max_chunk < 32) max_chunk = 32;
   
#ifdef STATBUF
   if(atype ==  C_DBL || atype == C_FLOAT){
      Ichunk=D_CHUNK, Kchunk=D_CHUNK, Jchunk=D_CHUNK;
   }else{
      Ichunk=ICHUNK; Kchunk=KCHUNK; Jchunk=JCHUNK;
   }
#else
   {
     Integer elems, factor = sizeof(DoubleComplex)/GAsizeofM(atype);
     Ichunk = Jchunk = Kchunk = CHUNK_SIZE;
     
     if ( max_chunk > Ichunk) {       
       /*if memory if very limited, performance degrades for large matrices
	 as chunk size is very small, which leads to communication overhead)*/
       Integer avail = pnga_memory_avail_type(atype);
       pnga_gop(pnga_type_f2c(MT_F_INT), &avail, (Integer)1, "min");
       if(avail<MINMEM && pnga_nodeid()==0) pnga_error("Not enough memory",avail);
       elems = (Integer)(avail*0.9);/* Donot use every last drop */
       
       max_chunk=GA_MIN(max_chunk, (Integer)(sqrt( (double)((elems-EXTRA)/3))));
       Ichunk = GA_MIN(m,max_chunk);
       Jchunk = GA_MIN(n,max_chunk);
       Kchunk = GA_MIN(k,max_chunk);
     }
     else /* "EXTRA" elems for safety - just in case */
       elems = 3*Ichunk*Jchunk + EXTRA*factor;

     a = (DoubleComplex*) ga_malloc(elems, atype, "GA mulmat bufs");     
     b = a + (Ichunk*Kchunk)/factor + 1; 
     c = b + (Kchunk*Jchunk)/factor + 1;
   }
#endif

   if(atype==C_DCPL){if((((DoubleComplex*)beta)->real == 0) &&
	       (((DoubleComplex*)beta)->imag ==0)) need_scaling =0;} 
   else if(atype==C_SCPL){if((((SingleComplex*)beta)->real == 0) &&
	       (((SingleComplex*)beta)->imag ==0)) need_scaling =0;} 
   else if((atype==C_DBL)){if(*(DoublePrecision *)beta == 0)need_scaling =0;}
   else if( *(float*)beta ==0) need_scaling =0;

   if(need_scaling) pnga_scale_patch(g_c, clo, chi, beta);
   else      pnga_fill_patch(g_c, clo, chi, beta);
  
   for(jlo = 0; jlo < n; jlo += Jchunk){ /* loop through columns of g_c patch */
       jhi = GA_MIN(n-1, jlo+Jchunk-1);
       jdim= jhi - jlo +1;
       
       for(klo = 0; klo < k; klo += Kchunk){    /* loop cols of g_a patch */
	 khi = GA_MIN(k-1, klo+Kchunk-1);        /* loop rows of g_b patch */
	 kdim= khi - klo +1;               

	 get_new_B = TRUE;
	 
	 for(ilo = 0; ilo < m; ilo += Ichunk){ /*loop through rows of g_c patch */
	   
	   if(ijk%nproc == iproc){
	     ihi = GA_MIN(m-1, ilo+Ichunk-1);
	     idim= cdim = ihi - ilo +1;
	     
	     if(atype == C_FLOAT) 
	       for (i = 0; i < idim*jdim; i++) *(((float*)c)+i)=0;
	     else if(atype ==  C_DBL)
	       for (i = 0; i < idim*jdim; i++) *(((double*)c)+i)=0;
	     else if(atype == C_SCPL)
	       for (i = 0; i < idim*jdim; i++){
                 ((SingleComplex*)c)[i].real=0;
                 ((SingleComplex*)c)[i].imag=0;
               }
             else
	       for (i = 0; i < idim*jdim; i++){ c[i].real=0;c[i].imag=0;}
	     
	     if (*transa == 'n' || *transa == 'N'){ 
	       adim = idim;
	       i0= ailo+ilo; i1= ailo+ihi;   
	       j0= ajlo+klo; j1= ajlo+khi;
	     }else{
	       adim = kdim;
	       i0= ajlo+klo; i1= ajlo+khi;   
	       j0= ailo+ilo; j1= ailo+ihi;
	     }

	     /* ga_get_(g_a, &i0, &i1, &j0, &j1, a, &adim); */
	     memcpy(tmplo,alo,arank*sizeof(Integer));
	     memcpy(tmphi,ahi,arank*sizeof(Integer));
	     SETINT(tmpld,1,arank-1);
	     tmplo[aipos]=i0; tmphi[aipos]=i1;
	     tmplo[ajpos]=j0; tmphi[ajpos]=j1;
	     tmpld[aipos]=i1-i0+1;
	     pnga_get(g_a,tmplo,tmphi,a,tmpld);
	     
	     if(get_new_B) {
	       if (*transb == 'n' || *transb == 'N'){ 
		 bdim = kdim;
		 i0= bilo+klo; i1= bilo+khi;   
		 j0= bjlo+jlo; j1= bjlo+jhi;
	       }else{
		 bdim = jdim;
		 i0= bjlo+jlo; i1= bjlo+jhi;   
		 j0= bilo+klo; j1= bilo+khi;
	       }
	       /* ga_get_(g_b, &i0, &i1, &j0, &j1, b, &bdim); */
	       memcpy(tmplo,blo,brank*sizeof(Integer));
	       memcpy(tmphi,bhi,brank*sizeof(Integer));
	       SETINT(tmpld,1,brank-1);
	       tmplo[bipos]=i0; tmphi[bipos]=i1;
	       tmplo[bjpos]=j0; tmphi[bjpos]=j1;
	       tmpld[bipos]=i1-i0+1;
	       pnga_get(g_b,tmplo,tmphi,b,tmpld);
	       get_new_B = FALSE;
	     }

	     idim_t=idim; jdim_t=jdim; kdim_t=kdim;
	     adim_t=adim; bdim_t=bdim; cdim_t=cdim;

		  switch(atype) {
		  case C_FLOAT:
		    BLAS_SGEMM(transa, transb, &idim_t, &jdim_t, &kdim_t,
			   (Real *)alpha, (Real *)a, &adim_t,
               (Real *)b, &bdim_t, (Real *)&ONE_CF,
               (Real *)c, &cdim_t);
		    break;
		  case C_DBL:
		    BLAS_DGEMM(transa, transb, &idim_t, &jdim_t, &kdim_t,
			   (DoublePrecision *)alpha, (DoublePrecision *)a, &adim_t,
               (DoublePrecision *)b, &bdim_t, (DoublePrecision *)&ONE,
               (DoublePrecision *)c, &cdim_t);
		    break;
		  case C_DCPL:
		    BLAS_ZGEMM(transa, transb, &idim_t, &jdim_t, &kdim_t,
			   (DoubleComplex *)alpha, (DoubleComplex *)a, &adim_t,
               (DoubleComplex *)b, &bdim_t, (DoubleComplex *)&ONE,
               (DoubleComplex *)c, &cdim_t);
		    break;
		  case C_SCPL:
		    BLAS_CGEMM(transa, transb, &idim_t, &jdim_t, &kdim_t,
			   (SingleComplex *)alpha, (SingleComplex *)a, &adim_t,
               (SingleComplex *)b, &bdim_t, (SingleComplex *)&ONE_CF,
               (SingleComplex *)c, &cdim_t);
		    break;
		  default:
		    pnga_error("ga_matmul_patch: wrong data type", atype);
		  }

                  i0= cilo+ilo; i1= cilo+ihi;   j0= cjlo+jlo; j1= cjlo+jhi;
                  /* ga_acc_(g_c, &i0, &i1, &j0, &j1, (DoublePrecision*)c, 
                                            &cdim, (DoublePrecision*)&ONE); */
		  memcpy(tmplo,clo,crank*sizeof(Integer));
		  memcpy(tmphi,chi,crank*sizeof(Integer));
		  SETINT(tmpld,1,crank-1);
		  tmplo[cipos]=i0; tmphi[cipos]=i1;
		  tmplo[cjpos]=j0; tmphi[cjpos]=j1;
		  tmpld[cipos]=i1-i0+1;
		  if(atype == C_FLOAT || atype == C_SCPL) 
		    pnga_acc(g_c,tmplo,tmphi,(float *)c,tmpld, &ONE_CF);
		  else
		    pnga_acc(g_c,tmplo,tmphi,c,tmpld,(DoublePrecision*)&ONE);
               }
	   ++ijk;
	 }
       }
   }

#ifndef STATBUF
   ga_free(a);
#endif
   
   GA_POP_NAME;
   if(local_sync_end)pnga_sync(); 
}

/**
 * 1. remove STATBUF
 * 2. 
 */
