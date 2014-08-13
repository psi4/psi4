#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*************************************************************************\
 Purpose:   File global.nalg.c contains a set of linear algebra routines 
            that operate on n-dim global arrays in the SPMD mode. 
 
 Date: 10.22.98
 Author: Jarek Nieplocha
\************************************************************************/

 
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#include "message.h"
#include "globalp.h"
#include "armci.h"
#include "ga-papi.h"
#include "ga-wapi.h"

#ifdef MPI
extern ARMCI_Group* ga_get_armci_group_(int);
#endif

/* work arrays used in all routines */
static Integer dims[MAXDIM], ld[MAXDIM-1];
static Integer lo[MAXDIM],hi[MAXDIM];
static Integer one_arr[MAXDIM]={1,1,1,1,1,1,1};

#define GET_ELEMS(ndim,lo,hi,ld,pelems){\
int _i;\
      for(_i=0, *pelems = hi[ndim-1]-lo[ndim-1]+1; _i< ndim-1;_i++) {\
         if(ld[_i] != (hi[_i]-lo[_i]+1)) pnga_error("layout problem",_i);\
         *pelems *= hi[_i]-lo[_i]+1;\
      }\
}

#define GET_ELEMS_W_GHOSTS(ndim,lo,hi,ld,pelems){\
int _i;\
      for(_i=0, *pelems = hi[ndim-1]-lo[ndim-1]+1; _i< ndim-1;_i++) {\
         if(ld[_i] < (hi[_i]-lo[_i]+1))\
           pnga_error("layout problem with ghosts",_i);\
         *pelems *= hi[_i]-lo[_i]+1;\
      }\
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_zero = pnga_zero
#endif
void pnga_zero(Integer g_a)
{
  Integer ndim, type, me, elems, p_handle;
  Integer num_blocks;
  void *ptr;
  /*register Integer i;*/
  int local_sync_begin,local_sync_end;

  local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
  _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
  p_handle = pnga_get_pgroup(g_a);

  if(local_sync_begin) pnga_pgroup_sync(p_handle);

  me = pnga_pgroup_nodeid(p_handle);

  pnga_check_handle(g_a, "ga_zero");
  GA_PUSH_NAME("ga_zero");

  num_blocks = pnga_total_blocks(g_a);

  pnga_inquire(g_a, &type, &ndim, dims);
  if (num_blocks < 0) {
    pnga_distribution(g_a, me, lo, hi);

    if ( lo[0]> 0 ){ /* base index is 1: we get 0 if no elements stored on p */

      if (pnga_has_ghosts(g_a)) {
        pnga_zero_patch(g_a,lo,hi);
        return;
      }
      pnga_access_ptr(g_a, lo, hi, &ptr, ld);
      GET_ELEMS(ndim,lo,hi,ld,&elems);

      /* switch (type){ */
/*         int *ia; */
/*         double *da; */
/*         float *fa; */
/*         long *la; */
/*         long long *lla; */
/*         case C_INT: */
/*         ia = (int*)ptr; */
/*         for(i=0;i<elems;i++) ia[i]  = 0; */
/*         break; */
/*         case C_DCPL: */
/*         elems *=2; */
/*         case C_DBL: */
/*         da = (double*)ptr; */
/*         for(i=0;i<elems;i++) da[i] = 0; */
/*         break; */
/*         case C_SCPL: */
/*         elems *=2; */
/*         case C_FLOAT: */
/*         fa = (float*)ptr; */
/*         for(i=0;i<elems;i++) fa[i]  = 0; */
/*         break; */
/*         case C_LONG: */
/*         la = (long*)ptr; */
/*         for(i=0;i<elems;i++) la[i]  = 0; */
/*         break; */
/*         case C_LONGLONG: */
/*         lla = (long long*)ptr; */
/*         for(i=0;i<elems;i++) lla[i]  = 0; */
/*         break; */
/*         default: pnga_error(" wrong data type ",type); */
/*       } */
      memset(ptr, 0, GAsizeofM(type)*elems);

      /* release access to the data */
      pnga_release_update(g_a, lo, hi);
    } 
  } else {
    pnga_access_block_segment_ptr(g_a, me, &ptr, &elems);
/*     switch (type){ */
/*       int *ia; */
/*       double *da; */
/*       float *fa; */
/*       long *la; */
/*       long long *lla; */
/*       case C_INT: */
/*         ia = (int*)ptr; */
/*         for(i=0;i<elems;i++) ia[i]  = 0; */
/*         break; */
/*       case C_DCPL: */
/*         elems *=2; */
/*       case C_DBL: */
/*         da = (double*)ptr; */
/*         for(i=0;i<elems;i++) da[i] = 0; */
/*         break; */
/*       case C_SCPL: */
/*         elems *=2; */
/*       case C_FLOAT: */
/*         fa = (float*)ptr; */
/*         for(i=0;i<elems;i++) fa[i]  = 0; */
/*         break; */
/*       case C_LONG: */
/*       la = (long*)ptr; */
/*       for(i=0;i<elems;i++) la[i]  = 0; */
/*       break; */
/*       case C_LONGLONG: */
/*       lla = (long long*)ptr; */
/*       for(i=0;i<elems;i++) lla[i]  = 0; */
/*       break; */
/*       default: pnga_error(" wrong data type ",type); */
/*     } */
      memset(ptr, 0, GAsizeofM(type)*elems);

    /* release access to the data */
    pnga_release_update_block_segment(g_a, me);
  }
  if(local_sync_end)pnga_pgroup_sync(p_handle);
  GA_POP_NAME;
}



#if 0
/*\ COPY ONE GLOBAL ARRAY INTO ANOTHER
\*/
static void snga_copy_old(Integer g_a, Integer g_b)
{
Integer  ndim, ndimb, type, typeb, me, elems=0, elemsb=0;
Integer dimsb[MAXDIM];
void *ptr_a, *ptr_b;

   me = pnga_nodeid();

   GA_PUSH_NAME("ga_copy");

   if(g_a == g_b) pnga_error("arrays have to be different ", 0L);

   pnga_inquire(g_a,  &type, &ndim, dims);
   pnga_inquire(g_b,  &typeb, &ndimb, dimsb);

   if(type != typeb) pnga_error("types not the same", g_b);

   if(!pnga_compare_distr(g_a,g_b))

      pnga_copy_patch("n",g_a, one_arr, dims, g_b, one_arr, dimsb);

   else {

     pnga_sync();

     pnga_distribution(g_a, me, lo, hi);
     if(lo[0]>0){
        pnga_access_ptr(g_a, lo, hi, &ptr_a, ld);
        if (pnga_has_ghosts(g_a)) {
          GET_ELEMS_W_GHOSTS(ndim,lo,hi,ld,&elems);
        } else {
          GET_ELEMS(ndim,lo,hi,ld,&elems);
        }
     }

     pnga_distribution(g_b, me, lo, hi);
     if(lo[0]>0){
        pnga_access_ptr(g_b, lo, hi, &ptr_b, ld);
        if (pnga_has_ghosts(g_b)) {
          GET_ELEMS_W_GHOSTS(ndim,lo,hi,ld,&elems);
        } else {
          GET_ELEMS(ndim,lo,hi,ld,&elems);
        }
     }
  
     if(elems!= elemsb)pnga_error("inconsistent number of elements",elems-elemsb);

     if(elems>0){
        ARMCI_Copy(ptr_a, ptr_b, (int)elems*GAsizeofM(type));
        pnga_release(g_a,lo,hi);
        pnga_release(g_b,lo,hi);
     }

     pnga_sync();
   }

   GA_POP_NAME;
}
#endif



/*\ COPY ONE GLOBAL ARRAY INTO ANOTHER
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_copy = pnga_copy
#endif
void pnga_copy(Integer g_a, Integer g_b)
{
Integer  ndim, ndimb, type, typeb, me_a, me_b;
Integer dimsb[MAXDIM],i;
Integer a_grp, b_grp, anproc, bnproc;
Integer num_blocks_a, num_blocks_b;
Integer blocks[MAXDIM], block_dims[MAXDIM];
void *ptr_a, *ptr_b;
int local_sync_begin,local_sync_end,use_put;

   GA_PUSH_NAME("ga_copy");

   local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
   _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
   a_grp = pnga_get_pgroup(g_a);
   b_grp = pnga_get_pgroup(g_b);
   me_a = pnga_pgroup_nodeid(a_grp);
   me_b = pnga_pgroup_nodeid(b_grp);
   anproc = pnga_get_pgroup_size(a_grp);
   bnproc = pnga_get_pgroup_size(b_grp);
   num_blocks_a = pnga_total_blocks(g_a);
   num_blocks_b = pnga_total_blocks(g_b);
   if (anproc <= bnproc) {
     use_put = 1;
   } else {
     use_put = 0;
   }
   /*if (a_grp != b_grp)
     pnga_error("Both arrays must be defined on same group",0L); */
   if(local_sync_begin) {
     if (anproc <= bnproc) {
       pnga_pgroup_sync(a_grp);
     } else if (a_grp == pnga_pgroup_get_world() &&
                b_grp == pnga_pgroup_get_world()) {
       pnga_sync();
     } else {
       pnga_pgroup_sync(b_grp);
     }
   }

   if(g_a == g_b) pnga_error("arrays have to be different ", 0L);

   pnga_inquire(g_a,  &type, &ndim, dims);
   pnga_inquire(g_b,  &typeb, &ndimb, dimsb);

   if(type != typeb) pnga_error("types not the same", g_b);
   if(ndim != ndimb) pnga_error("dimensions not the same", ndimb);

   for(i=0; i< ndim; i++)if(dims[i]!=dimsb[i]) 
                          pnga_error("dimensions not the same",i);

   if ((pnga_is_mirrored(g_a) && pnga_is_mirrored(g_b)) ||
       (!pnga_is_mirrored(g_a) && !pnga_is_mirrored(g_b))) {
     /* Both global arrays are mirrored or both global arrays are not mirrored.
        Copy operation is straightforward */

     if (use_put) {
       if (num_blocks_a < 0) {
         pnga_distribution(g_a, me_a, lo, hi);
         if(lo[0]>0){
           pnga_access_ptr(g_a, lo, hi, &ptr_a, ld);
           pnga_put(g_b, lo, hi, ptr_a, ld);
         }
       } else {
         if (!pnga_uses_proc_grid(g_a)) {
           for (i=me_a; i<num_blocks_a; i += anproc) {
             pnga_distribution(g_a, i, lo, hi);
             if (lo[0]>0) {
               pnga_access_block_ptr(g_a, i, &ptr_a, ld);
               pnga_put(g_b, lo, hi, ptr_a, ld);
             }
           }
         } else {
           Integer proc_index[MAXDIM], index[MAXDIM];
           Integer topology[MAXDIM], chk;
           pnga_get_proc_index(g_a, me_a, proc_index);
           pnga_get_proc_index(g_a, me_a, index);
           pnga_get_block_info(g_a, blocks, block_dims);
           pnga_get_proc_grid(g_a, topology);
           while (index[ndim-1] < blocks[ndim-1]) {
             /* find bounding coordinates of block */
             chk = 1;
             for (i = 0; i < ndim; i++) {
               lo[i] = index[i]*block_dims[i]+1;
               hi[i] = (index[i] + 1)*block_dims[i];
               if (hi[i] > dims[i]) hi[i] = dims[i];
               if (hi[i] < lo[i]) chk = 0;
             }
             if (chk) {
               pnga_access_block_grid_ptr(g_a, index, &ptr_a, ld);
               pnga_put(g_b, lo, hi, ptr_a, ld);
             }
             /* increment index to get next block on processor */
             index[0] += topology[0];
             for (i = 0; i < ndim; i++) {
               if (index[i] >= blocks[i] && i<ndim-1) {
                 index[i] = proc_index[i];
                 index[i+1] += topology[i+1];
               }
             }
           }
         }
       }
     } else {
       if (num_blocks_b < 0) {
         pnga_distribution(g_b, me_b, lo, hi);
         if(lo[0]>0){
           pnga_access_ptr(g_b, lo, hi, &ptr_b, ld);
           pnga_get(g_a, lo, hi, ptr_b, ld);
         }
       } else {
         if (!pnga_uses_proc_grid(g_a)) {
           for (i=me_b; i<num_blocks_b; i += bnproc) {
             pnga_distribution(g_b, i, lo, hi);
             if (lo[0]>0) {
               pnga_access_block_ptr(g_b, i, &ptr_b, ld);
               pnga_get(g_a, lo, hi, ptr_b, ld);
             }
           }
         } else {
           Integer proc_index[MAXDIM], index[MAXDIM];
           Integer topology[MAXDIM], chk;
           pnga_get_proc_index(g_b, me_b, proc_index);
           pnga_get_proc_index(g_b, me_b, index);
           pnga_get_block_info(g_b, blocks, block_dims);
           pnga_get_proc_grid(g_b, topology);
           while (index[ndim-1] < blocks[ndim-1]) {
             /* find bounding coordinates of block */
             chk = 1;
             for (i = 0; i < ndim; i++) {
               lo[i] = index[i]*block_dims[i]+1;
               hi[i] = (index[i] + 1)*block_dims[i];
               if (hi[i] > dims[i]) hi[i] = dims[i];
               if (hi[i] < lo[i]) chk = 0;
             }
             if (chk) {
               pnga_access_block_grid_ptr(g_b, index, &ptr_b, ld);
               pnga_get(g_a, lo, hi, ptr_b, ld);
             }
             /* increment index to get next block on processor */
             index[0] += topology[0];
             for (i = 0; i < ndim; i++) {
               if (index[i] >= blocks[i] && i<ndim-1) {
                 index[i] = proc_index[i];
                 index[i+1] += topology[i+1];
               }
             }
           }
         }
       }
     }
   } else {
     /* One global array is mirrored and the other is not */
     if (pnga_is_mirrored(g_a)) {
       /* Source array is mirrored and destination
          array is distributed. Assume source array is consistent */
       pnga_distribution(g_b, me_b, lo, hi);
       if (lo[0]>0) {
         pnga_access_ptr(g_b, lo, hi, &ptr_b, ld);
         pnga_get(g_a, lo, hi, ptr_b, ld);
       } 
     } else {
       /* source array is distributed and destination
          array is mirrored */
       pnga_zero(g_b);
       pnga_distribution(g_a, me_a, lo, hi);
       if (lo[0] > 0) {
         pnga_access_ptr(g_a, lo, hi, &ptr_a, ld);
         pnga_put(g_b, lo, hi, ptr_a, ld);
       }
       pnga_merge_mirrored(g_b);
     }
   }

   if(local_sync_end) {
     if (anproc <= bnproc) {
       pnga_pgroup_sync(a_grp);
     } else if (a_grp == pnga_pgroup_get_world() &&
                b_grp == pnga_pgroup_get_world()) {
       pnga_sync();
     } else {
       pnga_pgroup_sync(b_grp);
     }
   }
   GA_POP_NAME;
}



/*\ internal version of dot product
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_dot = pnga_dot
#endif
void pnga_dot(int Type, Integer g_a, Integer g_b, void *value)
{
Integer  ndim=0, type=0, atype=0, me=0, elems=0, elemsb=0;
register Integer i=0;
int isum=0;
long lsum=0;
long long llsum=0;
DoubleComplex zsum ={0.,0.};
SingleComplex csum ={0.,0.};
float fsum=0.0;
void *ptr_a=NULL, *ptr_b=NULL;
int alen=0;
Integer a_grp=0, b_grp=0;
Integer num_blocks_a=0, num_blocks_b=0;

Integer andim, adims[MAXDIM];
Integer bndim, bdims[MAXDIM];

   _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/

   GA_PUSH_NAME("ga_dot");
   a_grp = pnga_get_pgroup(g_a);
   b_grp = pnga_get_pgroup(g_b);
   if (a_grp != b_grp)
     pnga_error("Both arrays must be defined on same group",0L);
   me = pnga_pgroup_nodeid(a_grp);

   /* Check to see if either GA is block cyclic distributed */
   num_blocks_a = pnga_total_blocks(g_a);
   num_blocks_b = pnga_total_blocks(g_b);
   if (num_blocks_a >= 0 || num_blocks_b >= 0) {
     pnga_inquire(g_a, &type, &andim, adims);
     pnga_inquire(g_b, &type, &bndim, bdims);
     pnga_dot_patch(g_a, "n", one_arr, adims, g_b, "n", one_arr, bdims,
         value);
     GA_POP_NAME;
     return;
   }

   if(pnga_compare_distr(g_a,g_b) == FALSE ||
      pnga_has_ghosts(g_a) || pnga_has_ghosts(g_b)) {
       /* distributions not identical */
       pnga_inquire(g_a, &type, &andim, adims);
       pnga_inquire(g_b, &type, &bndim, bdims);

       pnga_dot_patch(g_a, "n", one_arr, adims, g_b, "n", one_arr, bdims,
                      value);
       
       GA_POP_NAME;
       return;
   }
   
   pnga_pgroup_sync(a_grp);
   pnga_inquire(g_a,  &type, &ndim, dims);
   if(type != Type) pnga_error("type not correct", g_a);
   pnga_distribution(g_a, me, lo, hi);
   if(lo[0]>0){
      pnga_access_ptr(g_a, lo, hi, &ptr_a, ld);
      if (pnga_has_ghosts(g_a)) {
        GET_ELEMS_W_GHOSTS(ndim,lo,hi,ld,&elems);
      } else {
        GET_ELEMS(ndim,lo,hi,ld,&elems);
      }
   }

   if(g_a == g_b){
     elemsb = elems;
     ptr_b = ptr_a;
   }else {  
     pnga_inquire(g_b,  &type, &ndim, dims);
     if(type != Type) pnga_error("type not correct", g_b);
     pnga_distribution(g_b, me, lo, hi);
     if(lo[0]>0){
        pnga_access_ptr(g_b, lo, hi, &ptr_b, ld);
        if (pnga_has_ghosts(g_b)) {
          GET_ELEMS_W_GHOSTS(ndim,lo,hi,ld,&elemsb);
        } else {
          GET_ELEMS(ndim,lo,hi,ld,&elemsb);
        }
     }
   }

   if(elems!= elemsb)pnga_error("inconsistent number of elements",elems-elemsb); 


      /* compute "local" contribution to the dot product */
      switch (type){
	int *ia, *ib;
	double *da,*db;
        float *fa, *fb;
        long *la,*lb;
        long long *lla,*llb;
        case C_INT:
           ia = (int*)ptr_a;
           ib = (int*)ptr_b;
           for(i=0;i<elems;i++) 
                 isum += ia[i]  * ib[i];
           *(int*)value = isum; 
           type = C_INT;
           alen = 1;
           break;

        case C_DCPL:
           for(i=0;i<elems;i++){
               DoubleComplex a = ((DoubleComplex*)ptr_a)[i];
               DoubleComplex b = ((DoubleComplex*)ptr_b)[i];
               zsum.real += a.real*b.real  - b.imag * a.imag;
               zsum.imag += a.imag*b.real  + b.imag * a.real;
           }
           *(DoubleComplex*)value = zsum; 
           type = C_DCPL;
           alen = 2;
           break;

        case C_SCPL:
           for(i=0;i<elems;i++){
               SingleComplex a = ((SingleComplex*)ptr_a)[i];
               SingleComplex b = ((SingleComplex*)ptr_b)[i];
               csum.real += a.real*b.real  - b.imag * a.imag;
               csum.imag += a.imag*b.real  + b.imag * a.real;
           }
           *(SingleComplex*)value = csum; 
           type = C_SCPL;
           alen = 2;
           break;

        case C_DBL:
           da = (double*)ptr_a;
           db = (double*)ptr_b;
           for(i=0;i<elems;i++) 
                 zsum.real += da[i]  * db[i];
           *(double*)value = zsum.real; 
           type = C_DBL;
           alen = 1;
           break;
        case C_FLOAT:
           fa = (float*)ptr_a;
           fb = (float*)ptr_b;
           for(i=0;i<elems;i++)
                 fsum += fa[i]  * fb[i];
           *(float*)value = fsum;
           type = C_FLOAT;
           alen = 1;
           break;         
        case C_LONG:
           la = (long*)ptr_a;
           lb = (long*)ptr_b;
           for(i=0;i<elems;i++)
		lsum += la[i]  * lb[i];
           *(long*)value = lsum;
           type = C_LONG;
           alen = 1;
           break;               
        case C_LONGLONG:
           lla = (long long*)ptr_a;
           llb = (long long*)ptr_b;
           for(i=0;i<elems;i++)
		llsum += lla[i]  * llb[i];
           *(long long*)value = llsum;
           type = C_LONGLONG;
           alen = 1;
           break;               
        default: pnga_error(" wrong data type ",type);
      }
   
      /* release access to the data */
      if(elems>0){
         pnga_release(g_a, lo, hi);
         if(g_a != g_b)pnga_release(g_b, lo, hi);
      }

    /*convert from C data type to ARMCI type */
    switch(type) {
      case C_FLOAT: atype=ARMCI_FLOAT; break;
      case C_DBL: atype=ARMCI_DOUBLE; break;
      case C_INT: atype=ARMCI_INT; break;
      case C_LONG: atype=ARMCI_LONG; break;
      case C_LONGLONG: atype=ARMCI_LONG_LONG; break;
      case C_DCPL: atype=ARMCI_DOUBLE; break;
      case C_SCPL: atype=ARMCI_FLOAT; break;
      default: pnga_error("pnga_dot: type not supported",type);
    }

   if (pnga_is_mirrored(g_a) && pnga_is_mirrored(g_b)) {
     armci_msg_gop_scope(SCOPE_NODE,value,alen,"+",atype);
   } else {
     if (a_grp == -1) {
       armci_msg_gop_scope(SCOPE_ALL,value,alen,"+",atype);
#ifdef MPI
     } else {
       armci_msg_group_gop_scope(SCOPE_ALL,value,alen,"+",atype,
           ga_get_armci_group_((int)a_grp));
#endif
     }
   }
    
   GA_POP_NAME;

}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_scale = pnga_scale
#endif
void pnga_scale(Integer g_a, void* alpha)
{
  Integer ndim, type, me, elems, grp_id;
  register Integer i;
  Integer num_blocks;
  void *ptr;
  int local_sync_begin,local_sync_end;

  local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
  _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
  grp_id = pnga_get_pgroup(g_a);
  if(local_sync_begin)pnga_pgroup_sync(grp_id);

  me = pnga_pgroup_nodeid(grp_id);

  pnga_check_handle(g_a, "ga_scale");
  GA_PUSH_NAME("ga_scale");
  num_blocks = pnga_total_blocks(g_a);

  pnga_inquire(g_a, &type, &ndim, dims);
  if (num_blocks < 0) {
    pnga_distribution(g_a, me, lo, hi);
    if (pnga_has_ghosts(g_a)) {
      pnga_scale_patch(g_a, lo, hi, alpha);
      return;
    }

    if ( lo[0]> 0 ){ /* base index is 1: we get 0 if no elements stored on p */

      pnga_access_ptr(g_a, lo, hi, &ptr, ld);
      GET_ELEMS(ndim,lo,hi,ld,&elems);

      switch (type){
        int *ia;
        double *da;
        DoubleComplex *ca, scale;
        SingleComplex *cfa, cfscale;
        long *la;
        long long *lla;
        float *fa;
        case C_INT:
        ia = (int*)ptr;
        for(i=0;i<elems;i++) ia[i]  *= *(int*)alpha;
        break;
        case C_LONG:
        la = (long*)ptr;
        for(i=0;i<elems;i++) la[i]  *= *(long*)alpha;
        break;
        case C_LONGLONG:
        lla = (long long*)ptr;
        for(i=0;i<elems;i++) lla[i]  *= *(long long*)alpha;
        break;
        case C_DCPL:
        ca = (DoubleComplex*)ptr;
        scale= *(DoubleComplex*)alpha;
        for(i=0;i<elems;i++){
          DoubleComplex val = ca[i]; 
          ca[i].real = scale.real*val.real  - val.imag * scale.imag;
          ca[i].imag = scale.imag*val.real  + val.imag * scale.real;
        }
        break;
        case C_SCPL:
        cfa = (SingleComplex*)ptr;
        cfscale= *(SingleComplex*)alpha;
        for(i=0;i<elems;i++){
          SingleComplex val = cfa[i]; 
          cfa[i].real = cfscale.real*val.real  - val.imag * cfscale.imag;
          cfa[i].imag = cfscale.imag*val.real  + val.imag * cfscale.real;
        }
        break;
        case C_DBL:
        da = (double*)ptr;
        for(i=0;i<elems;i++) da[i] *= *(double*)alpha;
        break;
        case C_FLOAT:
        fa = (float*)ptr;
        for(i=0;i<elems;i++) fa[i]  *= *(float*)alpha;
        break;       
        default: pnga_error(" wrong data type ",type);
      }

      /* release access to the data */
      pnga_release_update(g_a, lo, hi);
    }
  } else {
    pnga_access_block_segment_ptr(g_a, me, &ptr, &elems);
    switch (type){
      int *ia;
      double *da;
      DoubleComplex *ca, scale;
      SingleComplex *cfa, cfscale;
      long *la;
      long long *lla;
      float *fa;
      case C_INT:
      ia = (int*)ptr;
      for(i=0;i<elems;i++) ia[i]  *= *(int*)alpha;
      break;
      case C_LONG:
      la = (long*)ptr;
      for(i=0;i<elems;i++) la[i]  *= *(long*)alpha;
      break;
      case C_LONGLONG:
      lla = (long long*)ptr;
      for(i=0;i<elems;i++) lla[i]  *= *(long long*)alpha;
      break;
      case C_DCPL:
      ca = (DoubleComplex*)ptr;
      scale= *(DoubleComplex*)alpha;
      for(i=0;i<elems;i++){
        DoubleComplex val = ca[i]; 
        ca[i].real = scale.real*val.real  - val.imag * scale.imag;
        ca[i].imag = scale.imag*val.real  + val.imag * scale.real;
      }
      break;
      case C_SCPL:
      cfa = (SingleComplex*)ptr;
      cfscale= *(SingleComplex*)alpha;
      for(i=0;i<elems;i++){
        SingleComplex val = cfa[i]; 
        cfa[i].real = cfscale.real*val.real  - val.imag * cfscale.imag;
        cfa[i].imag = cfscale.imag*val.real  + val.imag * cfscale.real;
      }
      break;
      case C_DBL:
      da = (double*)ptr;
      for(i=0;i<elems;i++) da[i] *= *(double*)alpha;
      break;
      case C_FLOAT:
      fa = (float*)ptr;
      for(i=0;i<elems;i++) fa[i]  *= *(float*)alpha;
      break;       
      default: pnga_error(" wrong data type ",type);
    }
    /* release access to the data */
    pnga_release_update_block_segment(g_a, me);
  }
  GA_POP_NAME;
  if(local_sync_end)pnga_pgroup_sync(grp_id); 
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_add = pnga_add
#endif
void pnga_add(void *alpha, Integer g_a, void* beta, Integer g_b, Integer g_c)
{
Integer  ndim, type, typeC, me, elems=0, elemsb=0, elemsa=0;
register Integer i;
void *ptr_a, *ptr_b, *ptr_c;
Integer a_grp, b_grp, c_grp;
int local_sync_begin,local_sync_end;

 Integer andim, adims[MAXDIM];
 Integer bndim, bdims[MAXDIM];
 Integer cndim, cdims[MAXDIM];
 
   local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
   _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/

   GA_PUSH_NAME("ga_add");
   a_grp = pnga_get_pgroup(g_a);
   b_grp = pnga_get_pgroup(g_b);
   c_grp = pnga_get_pgroup(g_c);
   if (a_grp != b_grp || b_grp != c_grp)
     pnga_error("All three arrays must be on same group for ga_add",0L);

   if(local_sync_begin)pnga_pgroup_sync(a_grp);

   me = pnga_pgroup_nodeid(a_grp);
   if((pnga_compare_distr(g_a,g_b) == FALSE) ||
      (pnga_compare_distr(g_a,g_c) == FALSE) ||
       pnga_has_ghosts(g_a) || pnga_has_ghosts(g_b) || pnga_has_ghosts(g_c) ||
       pnga_total_blocks(g_a) > 0 || pnga_total_blocks(g_b) > 0 ||
       pnga_total_blocks(g_c) > 0) {
       /* distributions not identical */
       pnga_inquire(g_a, &type, &andim, adims);
       pnga_inquire(g_b, &type, &bndim, bdims);
       pnga_inquire(g_b, &type, &cndim, cdims);

       pnga_add_patch(alpha, g_a, one_arr, adims, beta, g_b, one_arr, bdims,
                      g_c, one_arr, cdims);
       
       GA_POP_NAME;
       return;
   }

   pnga_pgroup_sync(a_grp);
   pnga_inquire(g_c,  &typeC, &ndim, dims);
   pnga_distribution(g_c, me, lo, hi);
   if (  lo[0]>0 ){
     pnga_access_ptr(g_c, lo, hi, &ptr_c, ld);
     GET_ELEMS(ndim,lo,hi,ld,&elems);
   }

   if(g_a == g_c){
     ptr_a  = ptr_c;
     elemsa = elems;
   }else { 
     pnga_inquire(g_a,  &type, &ndim, dims);
     if(type != typeC) pnga_error("types not consistent", g_a);
     pnga_distribution(g_a, me, lo, hi);
     if (  lo[0]>0 ){
       pnga_access_ptr(g_a, lo, hi, &ptr_a, ld);
       GET_ELEMS(ndim,lo,hi,ld,&elemsa);
     }
   }

   if(g_b == g_c){
     ptr_b  = ptr_c;
     elemsb = elems;
   }else {
     pnga_inquire(g_b,  &type, &ndim, dims);
     if(type != typeC) pnga_error("types not consistent", g_b);
     pnga_distribution(g_b, me, lo, hi);
     if (  lo[0]>0 ){
       pnga_access_ptr(g_b, lo, hi, &ptr_b, ld);
       GET_ELEMS(ndim,lo,hi,ld,&elemsb);
     }
   }

   if(elems!= elemsb)pnga_error("inconsistent number of elements a",elems-elemsb);
   if(elems!= elemsa)pnga_error("inconsistent number of elements b",elems-elemsa);

   if (  lo[0]>0 ){

       /* operation on the "local" piece of data */
       switch(type){
         int *ia, *ib, *ic;
         double *da,*db,*dc;
         float *fa, *fb, *fc;
         long *la,*lb,*lc;
         long long *lla,*llb,*llc;
         case C_DBL:
                  da = (double*)ptr_a;
                  db = (double*)ptr_b;
                  dc = (double*)ptr_c;
                  for(i=0; i<elems; i++)
                      dc[i] = *(double*)alpha *da[i] +
                              *(double*)beta * db[i];
              break;
         case C_DCPL:
                  for(i=0; i<elems; i++){
                     DoubleComplex a = ((DoubleComplex*)ptr_a)[i];
                     DoubleComplex b = ((DoubleComplex*)ptr_b)[i];
                     DoubleComplex *ac = (DoubleComplex*)ptr_c;
                     DoubleComplex x= *(DoubleComplex*)alpha;
                     DoubleComplex y= *(DoubleComplex*)beta;
                     /* c = x*a + y*b */
                     ac[i].real = x.real*a.real - 
                              x.imag*a.imag + y.real*b.real - y.imag*b.imag;
                     ac[i].imag = x.real*a.imag + 
                              x.imag*a.real + y.real*b.imag + y.imag*b.real;
                  }
              break;
         case C_SCPL:
                  for(i=0; i<elems; i++){
                     SingleComplex a = ((SingleComplex*)ptr_a)[i];
                     SingleComplex b = ((SingleComplex*)ptr_b)[i];
                     SingleComplex *ac = (SingleComplex*)ptr_c;
                     SingleComplex x= *(SingleComplex*)alpha;
                     SingleComplex y= *(SingleComplex*)beta;
                     /* c = x*a + y*b */
                     ac[i].real = x.real*a.real - 
                              x.imag*a.imag + y.real*b.real - y.imag*b.imag;
                     ac[i].imag = x.real*a.imag + 
                              x.imag*a.real + y.real*b.imag + y.imag*b.real;
                  }
              break;
         case C_FLOAT:
                  fa = (float*)ptr_a;
                  fb = (float*)ptr_b;
                  fc = (float*)ptr_c;
                  for(i=0; i<elems; i++)
                      fc[i] = *(float*)alpha *fa[i] + *(float*)beta *fb[i]; 
              break;
         case C_INT:
                  ia = (int*)ptr_a;
                  ib = (int*)ptr_b;
                  ic = (int*)ptr_c;
                  for(i=0; i<elems; i++) 
                      ic[i] = *(int*)alpha *ia[i] + *(int*)beta *ib[i];
              break;    
         case C_LONG:
                  la = (long*)ptr_a;
		  lb = (long*)ptr_b;
		  lc = (long*)ptr_c;
                  for(i=0; i<elems; i++)
                      lc[i] = *(long*)alpha *la[i] + *(long*)beta *lb[i];
              break;
         case C_LONGLONG:
                  lla = (long long*)ptr_a;
		  llb = (long long*)ptr_b;
		  llc = (long long*)ptr_c;
                  for(i=0; i<elems; i++)
                      llc[i] = ( *(long long*)alpha *lla[i] +
                                 *(long long*)beta  * llb[i] );
                
       }

       /* release access to the data */
       pnga_release_update(g_c, lo, hi);
       if(g_c != g_a)pnga_release(g_a, lo, hi);
       if(g_c != g_b)pnga_release(g_b, lo, hi);
   }


   GA_POP_NAME;
   if(local_sync_end)pnga_pgroup_sync(a_grp);
}


static 
void snga_local_transpose(Integer type, char *ptra, Integer n, Integer stride, char *ptrb)
{
int i;
    switch(type){

       case C_INT:
            for(i = 0; i< n; i++, ptrb+= stride) 
               *(int*)ptrb= ((int*)ptra)[i];
            break;
       case C_DCPL:
            for(i = 0; i< n; i++, ptrb+= stride) 
               *(DoubleComplex*)ptrb= ((DoubleComplex*)ptra)[i];
            break;
       case C_SCPL:
            for(i = 0; i< n; i++, ptrb+= stride) 
               *(SingleComplex*)ptrb= ((SingleComplex*)ptra)[i];
            break;
       case C_DBL:
            for(i = 0; i< n; i++, ptrb+= stride) 
               *(double*)ptrb= ((double*)ptra)[i];
            break;
       case C_FLOAT:
            for(i = 0; i< n; i++, ptrb+= stride)
               *(float*)ptrb= ((float*)ptra)[i];
            break;      
       case C_LONG:
            for(i = 0; i< n; i++, ptrb+= stride)
               *(long*)ptrb= ((long*)ptra)[i];
            break;                                 
       case C_LONGLONG:
            for(i = 0; i< n; i++, ptrb+= stride)
               *(long long*)ptrb= ((long long*)ptra)[i];
            break;                                 
       default: pnga_error("bad type:",type);
    }
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_transpose = pnga_transpose
#endif
void pnga_transpose(Integer g_a, Integer g_b)
{
Integer me = pnga_nodeid();
Integer nproc = pnga_nnodes(); 
Integer atype, btype, andim, adims[MAXDIM], bndim, bdims[MAXDIM];
Integer lo[2],hi[2];
int local_sync_begin,local_sync_end;
Integer num_blocks_a;
char *ptr_tmp, *ptr_a;

    GA_PUSH_NAME("ga_transpose");
    
    local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
    _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
    if(local_sync_begin)pnga_sync();

    if(g_a == g_b) pnga_error("arrays have to be different ", 0L);

    pnga_inquire(g_a, &atype, &andim, adims);
    pnga_inquire(g_b, &btype, &bndim, bdims);

    if(bndim != 2 || andim != 2) pnga_error("dimension must be 2",0);
    if(atype != btype ) pnga_error("array type mismatch ", 0L);

    num_blocks_a = pnga_total_blocks(g_a);

    if (num_blocks_a < 0) {
      pnga_distribution(g_a, me, lo, hi);

      if(lo[0]>0){
        Integer nelem, lob[2], hib[2], nrow, ncol;
        int i, size=GAsizeofM(atype);

        nrow   = hi[0] -lo[0]+1;
        ncol   = hi[1] -lo[1]+1; 
        nelem  = nrow*ncol;
        lob[0] = lo[1]; lob[1] = lo[0];
        hib[0] = hi[1]; hib[1] = hi[0];

        /* allocate memory for transposing elements locally */
        ptr_tmp = (char *) ga_malloc(nelem, atype, "transpose_tmp");

        /* get access to local data */
        pnga_access_ptr(g_a, lo, hi, &ptr_a, ld);

        for(i = 0; i < ncol; i++){
          char *ptr = ptr_tmp + i*size;

          snga_local_transpose(atype, ptr_a, nrow, ncol*size, ptr);
          ptr_a += ld[0]*size;
        }

        pnga_release(g_a, lo, hi); 

        pnga_put(g_b, lob, hib, ptr_tmp ,&ncol);

        ga_free(ptr_tmp);
      }
    } else {
      Integer idx;
      Integer blocks[MAXDIM], block_dims[MAXDIM];
      Integer nelem, lob[2], hib[2], nrow, ncol;
      int i, size=GAsizeofM(atype);

      /* allocate memory for transposing elements locally */
      pnga_get_block_info(g_a, blocks, block_dims);

      /* Simple block-cyclic data distribution */
      nelem = block_dims[0]*block_dims[1];
      ptr_tmp = (char *) ga_malloc(nelem, atype, "transpose_tmp");
      if (!pnga_uses_proc_grid(g_a)) {
        for (idx = me; idx < num_blocks_a; idx += nproc) {
          pnga_distribution(g_a, idx, lo, hi);
          pnga_access_block_ptr(g_a, idx, &ptr_a, ld);

          nrow   = hi[0] -lo[0]+1;
          ncol   = hi[1] -lo[1]+1; 
          nelem  = nrow*ncol;
          lob[0] = lo[1]; lob[1] = lo[0];
          hib[0] = hi[1]; hib[1] = hi[0];
          for(i = 0; i < ncol; i++){
            char *ptr = ptr_tmp + i*size;

            snga_local_transpose(atype, ptr_a, nrow, ncol*size, ptr);
            ptr_a += ld[0]*size;
          }
          pnga_put(g_b, lob, hib, ptr_tmp ,&ncol);

          pnga_release_update_block(g_a, idx);
        }
      } else {
        /* Uses scalapack block-cyclic data distribution */
        Integer chk;
        Integer proc_index[MAXDIM], index[MAXDIM];
        Integer topology[MAXDIM], ichk;

        pnga_get_proc_index(g_a, me, proc_index);
        pnga_get_proc_index(g_a, me, index);
        pnga_get_block_info(g_a, blocks, block_dims);
        pnga_get_proc_grid(g_a, topology);
        /* Verify that processor has data */
        ichk = 1;
        for (i=0; i<andim; i++) {
          if (index[i]<0 || index[i] >= blocks[i]) {
            ichk = 0;
          }
        }

        if (ichk) {
          pnga_access_block_grid_ptr(g_a, index, &ptr_a, ld);
          while (index[andim-1] < blocks[andim-1]) {
            /* find bounding coordinates of block */
            chk = 1;
            for (i = 0; i < andim; i++) {
              lo[i] = index[i]*block_dims[i]+1;
              hi[i] = (index[i] + 1)*block_dims[i];
              if (hi[i] > adims[i]) hi[i] = adims[i];
              if (hi[i] < lo[i]) chk = 0;
            }
            if (chk) {
              pnga_access_block_grid_ptr(g_a, index, &ptr_a, ld);
              nrow   = hi[0] -lo[0]+1;
              ncol   = hi[1] -lo[1]+1; 
              nelem  = nrow*ncol;
              lob[0] = lo[1]; lob[1] = lo[0];
              hib[0] = hi[1]; hib[1] = hi[0];
              for(i = 0; i < ncol; i++){
                char *ptr = ptr_tmp + i*size;
                snga_local_transpose(atype, ptr_a, nrow, block_dims[0]*size, ptr);
                ptr_a += ld[0]*size;
              }
              pnga_put(g_b, lob, hib, ptr_tmp ,&block_dims[0]);
              pnga_release_update_block_grid(g_a, index);
            }
            /* increment index to get next block on processor */
            index[0] += topology[0];
            for (i = 0; i < andim; i++) {
              if (index[i] >= blocks[i] && i<andim-1) {
                index[i] = proc_index[i];
                index[i+1] += topology[i+1];
              }
            }
          }
        }
      }
      ga_free(ptr_tmp);
    }

    if(local_sync_end)pnga_sync();
    GA_POP_NAME;
}
