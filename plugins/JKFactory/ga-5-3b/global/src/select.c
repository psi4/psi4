#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#include "globalp.h"
#include "message.h"
#include "ga-papi.h"
#include "ga-wapi.h"
  
#define GET_ELEMS(ndim,lo,hi,ld,pelems){\
int _i;\
      for(_i=0, *pelems = hi[ndim-1]-lo[ndim-1]+1; _i< ndim-1;_i++) {\
         if(ld[_i] != (hi[_i]-lo[_i]+1)) pnga_error("layout problem",_i);\
         *pelems *= hi[_i]-lo[_i]+1;\
      }\
}

typedef struct { 
        union val_t {double dval; int ival; long lval; long long llval; float fval;}v; 
        Integer subscr[MAXDIM];
        DoubleComplex extra;
        SingleComplex extra2;      
} elem_info_t;

static void snga_select_elem(Integer type, char* op, void *ptr, Integer elems, elem_info_t *info,
                      Integer *ind)
{
  Integer i;
  switch (type){
    int *ia,ival;
    double *da,dval;
    DoubleComplex *ca;
    SingleComplex *cfa;
    float *fa,fval;
    long *la,lval;
    long long *lla,llval;

    case C_INT:
    ia = (int*)ptr;
    ival = *ia;
    if (strncmp(op,"min",3) == 0)
      for(i=0;i<elems;i++){ if(ival > ia[i]) {ival=ia[i];*ind=i; } } 
    else
      for(i=0;i<elems;i++){ if(ival < ia[i]) {ival=ia[i];*ind=i; } }

    info->v.ival = (int) ival;
    break;

    case C_DCPL:
    ca = (DoubleComplex*)ptr;
    dval=ca->real*ca->real + ca->imag*ca->imag;
    if (strncmp(op,"min",3) == 0)
      for(i=0;i<elems;i++, ca+=1 ){
        DoublePrecision tmp = ca->real*ca->real + ca->imag*ca->imag; 
        if(dval > tmp){dval = tmp; *ind = i;}
      }
    else
      for(i=0;i<elems;i++, ca+=1 ){
        DoublePrecision tmp = ca->real*ca->real + ca->imag*ca->imag; 
        if(dval < tmp){dval = tmp; *ind = i;}
      }

    info->v.dval = dval; /* use abs value  for comparison*/
    info->extra = ((DoubleComplex*)ptr)[*ind]; /* append the actual val */
    break;

    case C_SCPL:
       cfa = (SingleComplex*)ptr;
       fval=cfa->real*cfa->real + cfa->imag*cfa->imag;
       if (strncmp(op,"min",3) == 0)
          for(i=0;i<elems;i++, cfa+=1 ){
             float tmp = cfa->real*cfa->real + cfa->imag*cfa->imag;
             if(fval > tmp){fval = tmp; *ind = i;}
          }
       else
          for(i=0;i<elems;i++, cfa+=1 ){
             float tmp = cfa->real*cfa->real + cfa->imag*cfa->imag;
             if(fval < tmp){fval = tmp; *ind = i;}
          }

       info->v.fval = fval; /* use abs value  for comparison*/
       info->extra2 = ((SingleComplex*)ptr)[*ind]; /* append the actual val */
       break;
                                                               
    case C_DBL:
    da = (double*)ptr;
    dval = *da;
    if (strncmp(op,"min",3) == 0)
      for(i=0;i<elems;i++){ if(dval > da[i]) {dval=da[i];*ind=i; } }
    else
      for(i=0;i<elems;i++){ if(dval < da[i]) {dval=da[i];*ind=i; } }

    info->v.dval = dval; 
    break;

    case C_FLOAT:
    fa = (float*)ptr;
    fval = *fa;

    if (strncmp(op,"min",3) == 0)
      for(i=0;i<elems;i++){ if(fval > fa[i]) {fval=fa[i];*ind=i; } }
    else
      for(i=0;i<elems;i++){ if(fval < fa[i]) {fval=fa[i];*ind=i; } }

    info->v.fval = fval;
    break;
    case C_LONG:
    la = (long*)ptr;
    lval = *la;

    if (strncmp(op,"min",3) == 0)
      for(i=0;i<elems;i++){ if(lval > la[i]) {lval=la[i];*ind=i; } }
    else
      for(i=0;i<elems;i++){ if(lval < la[i]) {lval=la[i];*ind=i; } }

    info->v.lval = lval;
    break;
    case C_LONGLONG:
    lla = (long long*)ptr;
    llval = *lla;

    if (strncmp(op,"min",3) == 0)
      for(i=0;i<elems;i++){ if(llval > lla[i]) {llval=lla[i];*ind=i; } }
    else
      for(i=0;i<elems;i++){ if(llval < lla[i]) {llval=lla[i];*ind=i; } }

    info->v.llval = llval;
    break;

    default: pnga_error(" wrong data type ",type);
  }
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_select_elem = pnga_select_elem
#endif
void pnga_select_elem(Integer g_a, char* op, void* val, Integer *subscript)
{
  Integer ndim, type, me, elems, ind=0, i;
  Integer lo[MAXDIM],hi[MAXDIM],dims[MAXDIM],ld[MAXDIM-1];
  elem_info_t info;
  Integer num_blocks;
  int     participate=0;
  int local_sync_begin;

  local_sync_begin = _ga_sync_begin; 
  _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
  if(local_sync_begin)pnga_sync();

  me = pnga_nodeid();

  pnga_check_handle(g_a, "ga_select_elem");
  GA_PUSH_NAME("ga_elem_op");

  if (strncmp(op,"min",3) == 0);
  else if (strncmp(op,"max",3) == 0);
  else pnga_error("operator not recognized",0);

  pnga_inquire(g_a, &type, &ndim, dims);
  num_blocks = pnga_total_blocks(g_a);

  if (num_blocks < 0) {
    pnga_distribution(g_a, me, lo, hi);

    if ( lo[0]> 0 ){ /* base index is 1: we get 0 if no elements stored on p */

      /******************* calculate local result ************************/
      void    *ptr;
      pnga_access_ptr(g_a, lo, hi, &ptr, ld);
      GET_ELEMS(ndim,lo,hi,ld,&elems);
      participate =1;

      /* select local element */
      snga_select_elem(type, op, ptr, elems, &info, &ind);

      /* release access to the data */
      pnga_release(g_a, lo, hi);

      /* determine element subscript in the ndim-array */
      for(i = 0; i < ndim; i++){
        int elems = (int)( hi[i]-lo[i]+1);
        info.subscr[i] = ind%elems + lo[i] ;
        ind /= elems;
      }
    } 
  } else {
    void *ptr;
    Integer j, offset, jtot, upper;
    Integer nproc = pnga_nnodes();
    pnga_access_block_segment_ptr(g_a, me, &ptr, &elems);
    if (elems > 0) {
      participate =1;

      /* select local element */
      snga_select_elem(type, op, ptr, elems, &info, &ind);

      /* release access to the data */
      pnga_release_block_segment(g_a, me);

      /* convert local index back into a global array index */
      if (!pnga_uses_proc_grid(g_a)) {
        offset = 0;
        for (i=me; i<num_blocks; i += nproc) {
          pnga_distribution(g_a, i, lo, hi);
          jtot = 1;
          for (j=0; j<ndim; j++) {
            jtot *= (hi[j]-lo[j]+1);
          }
          upper = offset + jtot;
          if (ind >= offset && ind < upper) {
            break;
          }  else {
            offset += jtot;
          }
        }
        /* determine element subscript in the ndim-array */
        ind -= offset;
        for(i = 0; i < ndim; i++){
          int elems = (int)( hi[i]-lo[i]+1);
          info.subscr[i] = ind%elems + lo[i] ;
          ind /= elems;
        }
      } else {
        Integer stride[MAXDIM], index[MAXDIM];
        Integer blocks[MAXDIM], block_dims[MAXDIM];
        Integer proc_index[MAXDIM], topology[MAXDIM];
        Integer l_index[MAXDIM];
        Integer min, max;
        pnga_get_proc_index(g_a, me, proc_index);
        pnga_get_block_info(g_a, blocks, block_dims);
        pnga_get_proc_grid(g_a, topology);
        /* figure out strides for locally held block of data */
        for (i=0; i<ndim; i++) {
          stride[i] = 0;
          for (j=proc_index[i]; j<blocks[i]; j += topology[i]) {
            min = j*block_dims[i] + 1;
            max = (j+1)*block_dims[i];
            if (max > dims[i])
              max = dims[i];
            stride[i] += (max - min + 1);
          }
        }
        /* use strides to figure out local index */
        l_index[0] = ind%stride[0];
        for (i=1; i<ndim; i++) {
          ind = (ind-l_index[i-1])/stride[i-1];
          l_index[i] = ind%stride[i];
        }
        /* figure out block index for block holding data element */
        for (i=0; i<ndim; i++) {
          index[i] = l_index[i]/block_dims[i];
        }
        for (i=0; i<ndim; i++) {
          lo[i] = (topology[i]*index[i] + proc_index[i])*block_dims[i];
          info.subscr[i] = l_index[i]%block_dims[i] + lo[i];
        }
      }
    }
  }
  /* calculate global result */
  if(type==C_INT){
    int size = sizeof(double) + sizeof(Integer)*(int)ndim;
    armci_msg_sel_scope(SCOPE_ALL,&info,size,op,ARMCI_INT,participate);
    *(int*)val = (int)info.v.ival;
  }else if(type==C_LONG){
    int size = sizeof(double) + sizeof(Integer)*(int)ndim;
    armci_msg_sel_scope(SCOPE_ALL,&info,size,op,ARMCI_LONG,participate);
    *(long*)val = info.v.lval;
  }else if(type==C_LONGLONG){
    int size = sizeof(double) + sizeof(Integer)*(int)ndim;
    armci_msg_sel_scope(SCOPE_ALL,&info,size,op,ARMCI_LONG_LONG,participate);
    *(long long*)val = info.v.llval;
  }else if(type==C_DBL){
    int size = sizeof(double) + sizeof(Integer)*(int)ndim;
    armci_msg_sel_scope(SCOPE_ALL,&info,size,op,ARMCI_DOUBLE,participate);
    *(DoublePrecision*)val = info.v.dval;
  }else if(type==C_FLOAT){
    int size = sizeof(double) + sizeof(Integer)*ndim;
    armci_msg_sel_scope(SCOPE_ALL,&info,size,op,ARMCI_FLOAT,participate);
    *(float*)val = info.v.fval;       
  }else if(type==C_SCPL){
    int size = sizeof(info); /* for simplicity we send entire info */
    armci_msg_sel_scope(SCOPE_ALL,&info,size,op,ARMCI_FLOAT,participate);
    *(SingleComplex*)val = info.extra2;
  }else{
    int size = sizeof(info); /* for simplicity we send entire info */
    armci_msg_sel_scope(SCOPE_ALL,&info,size,op,ARMCI_DOUBLE,participate);
    *(DoubleComplex*)val = info.extra;
  }

  for(i = 0; i < ndim; i++) subscript[i]= info.subscr[i];
  GA_POP_NAME;
}
