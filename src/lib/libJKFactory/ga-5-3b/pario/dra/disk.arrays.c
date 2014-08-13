#if HAVE_CONFIG_H
#   include "config.h"
#endif
/** @file
 * DRA operations with a buffer manager layer, modified by Bilash
 * The buffer manager provides functionalities related to buffers
 *
 * $Id: disk.arrays.c,v 1.79.2.2 2007-03-24 01:19:28 manoj Exp $
 *
 ************************** DISK ARRAYS **************************************
 *         Jarek Nieplocha, Fri May 12 11:26:38 PDT 1995                     *
 *****************************************************************************
 *
 * DISCLAIMER
 *
 * This material was prepared as an account of work sponsored by an
 * agency of the United States Government.  Neither the United States
 * Government nor the United States Department of Energy, nor Battelle,
 * nor any of their employees, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
 * ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
 * COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
 * SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT
 * INFRINGE PRIVATELY OWNED RIGHTS.
 *
 *
 * ACKNOWLEDGMENT
 *
 * This software and its documentation were produced with United States
 * Government support under Contract Number DE-AC06-76RLO-1830 awarded by
 * the United States Department of Energy.  The United States Government
 * retains a paid-up non-exclusive, irrevocable worldwide license to
 * reproduce, prepare derivative works, perform publicly and display
 * publicly by or for the US Government, including the right to
 * distribute to other US Government contractors.
 */
#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STDIO_H
#   include <stdio.h>
#endif

#include "buffers.h"
#include "dra.h"
#include "draf2c.h"
#include "drap.h"
#include "global.h"
#include "ga-papi.h"
#include "macdecls.h"

#define WALLTIME 0
#if WALLTIME
#   include "walltime.c"
#endif

/************************** constants ****************************************/

#define DRA_FAIL  (Integer)1
#define COLUMN    1
#define ROW       0
#define ON        1
#define OFF       0

#define ILO       0
#define IHI       1
#define JLO       2
#define JHI       3

#define DRA_OP_WRITE 777
#define DRA_OP_READ  888
#define PROBE 111
#define WAIT 222

#define MAX_REQ   5

/** message type/tag used by DRA */
#define  DRA_TYPE_GSM 32760 - 6
#define  DRA_TYPE_GOP 32760 - 7

#define INFINITE_NUM_PROCS  8094

#define CLIENT_TO_SERVER 2

/* #define DRA_DBLE_BUFFER */

#if defined(SP1)|| defined(SP) || defined(LAPI)
#   define DRA_NUM_IOPROCS 8 
#else
#   define DRA_NUM_IOPROCS 1
#endif

#ifndef DRA_NUM_FILE_MGR
#  define DRA_NUM_FILE_MGR DRA_NUM_IOPROCS
#endif

#define DRA_BUF_SIZE BUF_SIZE

#define DEF_MAX_ARRAYS 16
#define DRA_MAX_ARRAYS 1024
#define EPS_SEARCH 100
#define LOAD 1
#define STORE 2
#define TRANS 1
#define NOTRANS 0
/*#define DEBUG 1*/
/*#define CLEAR_BUF 1*/


/***************************** Global Data ***********************************/

/** Default number of IO procs */
Integer _dra_io_procs;
/** Default number of files to use on open file systems */
Integer _dra_number_of_files;
/** array of struct for basic info about DRA arrays*/
disk_array_t *DRA;

buf_context_t buf_ctxt; /**< buffer context handle */
int nbuf = 4; /**< number of buffers to be used */

Integer _max_disk_array; /**< max number of disk arrays open at a time */
logical dra_debug_flag;  /**< globally defined debug parameter */

request_t Requests[MAX_REQ];
int num_pending_requests=0;
Integer _dra_turn=0;
int Dra_num_serv=DRA_NUM_IOPROCS;


/****************************** Macros ***************************************/

#define dai_sizeofM(_type) MA_sizeof(_type, 1, MT_C_CHAR)

#define dai_check_typeM(_type)  if (_type != C_DBL && _type != C_INT \
     && _type != C_LONG && _type != C_DCPL && _type != C_FLOAT && _type != C_SCPL) \
                                  dai_error("invalid type ",_type)  
#define dai_check_handleM(_handle, msg)                                    \
{\
        if((_handle+DRA_OFFSET)>=_max_disk_array || (_handle+DRA_OFFSET)<0) \
        {fprintf(stderr,"%s, %ld --",msg, (long)_max_disk_array);\
        dai_error("invalid DRA handle",_handle);}                           \
        if( DRA[(_handle+DRA_OFFSET)].actv == 0)                            \
        {fprintf(stderr,"%s:",msg);\
        dai_error("disk array not active",_handle);}                       \
}
        
#define dai_check_rangeM(_lo, _hi, _dim, _err_msg)                         \
        if(_lo < (Integer)1   || _lo > _dim ||_hi < _lo || _hi > _dim)     \
        dai_error(_err_msg, _dim)
 
#define ga_get_sectM(sect, _buf, _ld)\
   pnga_get(sect.handle, sect.lo, sect.hi, _buf, &_ld)

#define ga_put_sectM(sect, _buf, _ld)\
   pnga_put(sect.handle, sect.lo, sect.hi, _buf, &_ld)

#define fill_sectionM(sect, _hndl, _ilo, _ihi, _jlo, _jhi) \
{ \
        sect.handle = _hndl;\
        sect.ndim   = 2; \
        sect.lo[0]    = _ilo;\
        sect.hi[0]    = _ihi;\
        sect.lo[1]    = _jlo;\
        sect.hi[1]    = _jhi;\
}

#define sect_to_blockM(ds_a, CR)\
{\
      Integer   hndl = (ds_a).handle+DRA_OFFSET;\
      Integer   br   = ((ds_a).lo[0]-1)/DRA[hndl].chunk[0];\
      Integer   bc   = ((ds_a).lo[1]-1)/DRA[hndl].chunk[1];\
      Integer   R    = (DRA[hndl].dims[0] + DRA[hndl].chunk[0] -1)/DRA[hndl].chunk[0];\
               *(CR) = bc * R + br;\
}

#define block_to_sectM(ds_a, CR)\
{\
      Integer   hndl = (ds_a)->handle+DRA_OFFSET;\
      Integer   R    = (DRA[hndl].dims[0] + DRA[hndl].chunk[0]-1)/DRA[hndl].chunk[0];\
      Integer   br = (CR)%R;\
      Integer   bc = ((CR) - br)/R;\
      (ds_a)->  lo[0]= br * DRA[hndl].chunk[0] +1;\
      (ds_a)->  lo[1]= bc * DRA[hndl].chunk[1] +1;\
      (ds_a)->  hi[0]= (ds_a)->lo[0] + DRA[hndl].chunk[0];\
      (ds_a)->  hi[1]= (ds_a)->lo[1] + DRA[hndl].chunk[1];\
      if( (ds_a)->hi[0] > DRA[hndl].dims[0]) (ds_a)->hi[0] = DRA[hndl].dims[0];\
      if( (ds_a)->hi[1] > DRA[hndl].dims[1]) (ds_a)->hi[1] = DRA[hndl].dims[1];\
}
      
#define INDEPFILES(x) (DRA[(x)+DRA_OFFSET].indep)

#define ITERATOR_2D(i,j, base, ds_chunk)\
    for(j = ds_chunk.lo[1], base=0, jj=0; j<= ds_chunk.hi[1]; j++,jj++)\
        for(i = ds_chunk.lo[0], ii=0; i<= ds_chunk.hi[0]; i++,ii++,base++)

#define COPY_SCATTER(ADDR_BASE, TYPE, ds_chunk)\
    ITERATOR_2D(i,j, base, ds_chunk) \
    ADDR_BASE[base+vindex] = ((TYPE*)buffer)[ldb*jj + ii]

#define COPY_GATHER(ADDR_BASE, TYPE, ds_chunk)\
    for(i=0; i< nelem; i++){\
        Integer ldc = ds_chunk.hi[0] - ds_chunk.lo[0]+1;\
        base = INT_MB[pindex+i]; jj = base/ldc; ii = base%ldc;\
        ((TYPE*)buffer)[ldb*jj + ii] = ADDR_BASE[i+vindex];\
    }

#define COPY_TYPE(OPERATION, MATYPE, ds_chunk)\
    switch(MATYPE){\
        case C_DBL:   COPY_ ## OPERATION(DBL_MB,double,ds_chunk);break;\
        case C_INT:   COPY_ ## OPERATION(INT_MB,int,ds_chunk);break;\
        case C_DCPL:  COPY_ ## OPERATION(DCPL_MB,DoubleComplex,ds_chunk);break;\
        case C_SCPL:  COPY_ ## OPERATION(SCPL_MB,SingleComplex,ds_chunk);break;\
        case C_FLOAT: COPY_ ## OPERATION(FLT_MB,float,ds_chunk);\
    }

#define nsect_to_blockM(ds_a, CR) \
{ \
  Integer hndl = (ds_a).handle+DRA_OFFSET;\
  Integer _i, _ndim = DRA[hndl].ndim; \
  Integer _R, _b; \
  *(CR) = 0; \
  _R = 0; \
  for (_i=_ndim-1; _i >= 0; _i--) { \
    _b = ((ds_a).lo[_i]-1)/DRA[hndl].chunk[_i]; \
    _R = (DRA[hndl].dims[_i]+DRA[hndl].chunk[_i]-1)/DRA[hndl].chunk[_i];\
    *(CR) = *(CR) * _R + _b; \
  } \
}


#define dai_dest_indices_1d_M(index, id, jd, ilod, jlod, ldd) \
{ \
    Integer _index_;\
    _index_ = (is)-(ilos);\
    *(id) = (_index_)%(ldd) + (ilod);\
    *(jd) = (_index_)/(ldd) + (jlod);\
}


#define dai_dest_indicesM(is, js, ilos, jlos, lds, id, jd, ilod, jlod, ldd)\
{ \
    Integer _index_;\
    _index_ = (lds)*((js)-(jlos)) + (is)-(ilos);\
    *(id) = (_index_)%(ldd) + (ilod);\
    *(jd) = (_index_)/(ldd) + (jlod);\
}


#define nga_get_sectM(sect, _buf, _ld, hdl)\
 if (hdl != NULL)\
 pnga_nbget(sect.handle, sect.lo, sect.hi, _buf, _ld, hdl);\
 else\
 pnga_get(sect.handle, sect.lo, sect.hi, _buf, _ld);


#define nga_put_sectM(sect, _buf, _ld, hdl)\
 if (hdl != NULL)\
 pnga_nbput(sect.handle, sect.lo, sect.hi, _buf, _ld, hdl);\
 else\
 pnga_put(sect.handle, sect.lo, sect.hi, _buf, _ld);


#define ndai_dest_indicesM(ds_chunk, ds_a, gs_chunk, gs_a)   \
{\
  Integer _i; \
  Integer _ndim = ds_a.ndim; \
  for (_i=0; _i<_ndim; _i++) { \
    gs_chunk.lo[_i] = gs_a.lo[_i] + ds_chunk.lo[_i]- ds_a.lo[_i]; \
    gs_chunk.hi[_i] = gs_a.lo[_i] + ds_chunk.hi[_i]- ds_a.lo[_i]; \
  } \
}


#define ndai_trnsp_dest_indicesM(ds_chunk, ds_a, gs_chunk, gs_a)   \
{\
  Integer _i; \
  Integer _ndim = ds_a.ndim; \
  for (_i=0; _i<_ndim; _i++) { \
    gs_chunk.lo[_ndim-1-_i] = gs_a.lo[_ndim-1-_i] \
                            + ds_chunk.lo[_i]- ds_a.lo[_i]; \
    gs_chunk.hi[_ndim-1-_i] = gs_a.lo[_ndim-1-_i] \
                            + ds_chunk.hi[_i]- ds_a.lo[_i]; \
  } \
}


/** Simple sort using straight insertion */
#define block_sortM(_ndim, _block_orig, _block_map) \
{\
  Integer _i,_j,_it,_bt; \
  Integer _block_tmp[MAXDIM]; \
  for (_i=0; _i < (_ndim); _i++) { \
    _block_map[_i] = _i; \
    _block_tmp[_i] = _block_orig[_i]; \
  } \
  for (_j=(_ndim)-2; _j >= 0; _j--) { \
    _i = _j + 1; \
    _bt = _block_tmp[_j]; \
    _it = _block_map[_j]; \
    while (_i < (_ndim) && _bt < _block_tmp[_i]) { \
      _block_tmp[_i-1] = _block_tmp[_i]; \
      _block_map[_i-1] = _block_map[_i]; \
      _i++; \
    }\
    _block_tmp[_i-1] = _bt; \
    _block_map[_i-1] = _it; \
  }\
}


#define nfill_sectionM(sect, _hndl, _ndim, _lo, _hi) \
{ \
  Integer _i; \
  sect.handle = _hndl; \
  sect.ndim = _ndim; \
  for (_i=0; _i<_ndim; _i++) { \
    sect.lo[_i]    = _lo[_i]; \
    sect.hi[_i]    = _hi[_i]; \
  } \
}


#define nblock_to_sectM(ds_a, _CR) \
{\
  Integer _i, _b[MAXDIM], _C = (_CR); \
  Integer _hndl = (ds_a)->handle+DRA_OFFSET; \
  Integer _R = (DRA[_hndl].dims[0]+DRA[_hndl].chunk[0]-1)/DRA[_hndl].chunk[0]; \
  (ds_a)->ndim = DRA[_hndl].ndim; \
  _b[0] = _C%_R; \
  for (_i=1; _i<DRA[_hndl].ndim; _i++) { \
    _C = (_C - _b[_i-1])/_R; \
    _R = (DRA[_hndl].dims[_i]+DRA[_hndl].chunk[_i]-1)/DRA[_hndl].chunk[_i]; \
    _b[_i] = (_C)%_R; \
  } \
  for (_i=0; _i<DRA[_hndl].ndim; _i++) { \
    (ds_a)->lo[_i] = _b[_i]*DRA[_hndl].chunk[_i] + 1; \
    (ds_a)->hi[_i] = (ds_a)->lo[_i] + DRA[_hndl].chunk[_i] - 1; \
    if ((ds_a)->hi[_i] > DRA[_hndl].dims[_i]) \
      (ds_a)->hi[_i] = DRA[_hndl].dims[_i]; \
  } \
}


#define nblock_to_indicesM(_index,_ndim,_block_dims,_CC) \
{ \
  Integer _i, _C=(_CC); \
  _index[0] = _C%_block_dims[0]; \
  for (_i=1; _i<(_ndim); _i++) { \
    _C = (_C - _index[_i-1])/_block_dims[_i-1]; \
    _index[_i] = _C%_block_dims[_i]; \
  } \
}


#define ndai_check_rangeM(_lo, _hi, _ndim, _dims, _err_msg) \
{ \
  int _range_ok = 1, _i; \
  for (_i=0; _i < (_ndim); _i++) { \
    if (_lo[_i] < 1 || _lo[_i] > _dims[_i] || _hi[_i] < _lo[_i] \
                || _hi[_i] > _dims[_i] ) _range_ok = 0; \
  } \
  if(!_range_ok) dai_error(_err_msg, _dim); \
}


char dummy_fname[DRA_MAX_FNAME];

/*****************************************************************************/

/**
 * determines if write operation to a disk array is allowed
 */
int dai_write_allowed(Integer d_a)
{
    Integer handle = d_a+DRA_OFFSET;
    if(DRA[handle].mode == DRA_W || DRA[handle].mode == DRA_RW) return 1;
    else return 0;
}


/**
 * determines if read operation from a disk array is allowed
 */
int dai_read_allowed(Integer d_a)
{
    Integer handle = d_a+DRA_OFFSET;
    if(DRA[handle].mode == DRA_R || DRA[handle].mode == DRA_RW) return 1;
    else return 0;
}


/**
 * number of processes that could perform I/O
 */
Integer dai_io_procs(Integer d_a)
{
    Integer handle = d_a+DRA_OFFSET;
    Integer num;

    /* this one of many possibilities -- depends on the system */
    /*
#ifdef _CRAYMPP
num = DRA_NUM_IOPROCS;
#else
num = (INDEPFILES(d_a)) ? INFINITE_NUM_PROCS: DRA_NUM_IOPROCS; 
#endif
*/
    if (INDEPFILES(d_a)) {
        num = pnga_cluster_nnodes();
    } else {
        num = DRA[handle].ioprocs;
    }

    return( PARIO_MIN( pnga_nnodes(), num));
}

/**
 * Translation of DRA create/opening modes to ELIO create/open
 * mode. DRA modes map directly to ELIO modes, except that write-only
 * DRAs are backed by read-write ELIO files.
 */
int dai_elio_mode(int dra_mode) {
  int emode = dra_mode; /* dra modes map to elio mode*/
  if(dra_mode == DRA_W) {
    /*except W, which translate to read-write files*/
    emode = ELIO_RW;
  }
  return emode;
}


/**
 * rank of calling process in group of processes that could perform I/O
 * a negative value means that this process doesn't do I/O
 */
Integer dai_io_nodeid(Integer d_a)
{
    Integer handle = d_a+DRA_OFFSET;
    Integer me = pnga_nodeid();
    Integer pid, id, nid, nnodes,nprocs;
    Integer nodeid = pnga_cluster_nodeid();
    Integer zero = 0;

    /* again, one of many possibilities: 
     * if proc id beyond I/O procs number, negate it
     */

    if (INDEPFILES(d_a)) {
        if(me == pnga_cluster_procid(nodeid, zero)) me = nodeid;
        else me = -1;
    } else {
        if (DRA[handle].ioprocs == 1) {
            if (me == 0) return me;
            else return -1;
        } else {
            nnodes = pnga_cluster_nnodes();
            nprocs = pnga_cluster_nprocs(nodeid);
            pid = me % nprocs;
            nid = (me - pid)/nprocs;
            id = pid * nnodes + nid;
            if (id < DRA[handle].ioprocs) return id;
            else return -1;
        }
    }

    /*        if (me >= dai_io_procs(d_a)) me = -me;*/
    return (me);
}


/**
 * determines if I/O process participates in file management (create/delete)
 */
Integer dai_io_manage(Integer d_a)
{
    Integer me = dai_io_nodeid(d_a);

    if(me >= 0 )
        return (1);
    else
        return (0);
}


/**
 * select one master process for each file associated with d_a
 */
Integer dai_file_master(Integer d_a)
{
    Integer handle = d_a+DRA_OFFSET;

    if(dai_io_nodeid(d_a)<0)return 0;

    /* for indep files each I/O process has its own file
     * for shared file 0 is the master
     */

    if(INDEPFILES(d_a) || DRA[handle].numfiles > 1 ||
            dai_io_nodeid(d_a) == 0 ) return 1;
    else return 0;

}


void dai_callback(int op, int transp, section_t gs_a, section_t ds_a, section_t ds_chunk, Integer ld[], char *buf, Integer req)
{
    int i;
    buf_info *bi;

    bi = (buf_info*) buf;
    if (bi->callback==ON)
        dai_error("DRA: callback not cleared for a buffer",0);

    bi->callback = ON;
    bi->args.op = op;
    bi->args.transp = transp;
    bi->args.gs_a = gs_a;
    bi->args.ds_a = ds_a;
    bi->args.ds_chunk = ds_chunk;
    for (i=0; i<gs_a.ndim-1; i++)
        bi->args.ld[i] = ld[i];


}


/**
 * function to release buffers by completing the transfers
 * this function will be passed on as a parameter to a buffer management layer
 */
void wait_buf(char *buf);


/**
 * INITIALIZE DISK ARRAY DATA STRUCTURES
 *
 * @param max_arrays[in]
 * @param max_array_size[in]
 * @param tot_disk_space[in]
 * @param max_memory[in]
 */
Integer FATR dra_init_(
        Integer *max_arrays,
        DoublePrecision *max_array_size,
        DoublePrecision *tot_disk_space,
        DoublePrecision *max_memory)
{
    int i, buf_size;
    pnga_sync();

    if(*max_arrays<-1 || *max_arrays> DRA_MAX_ARRAYS)
        dai_error("dra_init: incorrect max number of arrays",*max_arrays);
    _max_disk_array = (*max_arrays==-1) ? DEF_MAX_ARRAYS: *max_arrays;

    Dra_num_serv = drai_get_num_serv();

    DRA = (disk_array_t*)malloc(sizeof(disk_array_t)* (int)*max_arrays);
    if(!DRA) dai_error("dra_init: memory alocation failed\n",0);
    for(i=0; i<_max_disk_array ; i++)DRA[i].actv=0;

    for(i=0; i<MAX_REQ; i++) Requests[i].num_pending=0;

    dra_debug_flag = FALSE;
    _dra_io_procs = pnga_cluster_nnodes();
    _dra_number_of_files = pnga_cluster_nnodes();

    /* initialize Buffer Manager */
    buf_size = sizeof (buf_info) + (int) DBL_BUF_SIZE;
    buffer_init(&buf_ctxt, nbuf, buf_size, &wait_buf);

    pnga_sync();

    return(ELIO_OK);
}


/**
 * correct chunk size to fit the buffer and introduce allignment
 *
 * @param a[in,out] Initial guess for larger chunk
 * @param b[in,out] Initial guess for smaller chunk
 * @param prod[in]  Total number of elements that will fit in buffer [input]
 * @param ratio[in] Ratio of size of current patch formed by a and b
 *                  to size of I/O buffer
 */
void dai_correct_chunking(Integer* a, Integer* b, Integer prod, double ratio)
{
    Integer b0, bt, eps0, eps1, trial;
    double  da = (double)*a, db = (double)*b;
    double  ch_ratio = da/db;

    db = sqrt(prod /ch_ratio); 
    trial = b0 = (Integer)db; eps0 = prod%b0;
    trial = PARIO_MAX(trial,EPS_SEARCH); /******/
    /* search in the neighborhood for a better solution */
    for (bt = trial+ EPS_SEARCH; bt> PARIO_MIN(1,trial-EPS_SEARCH); bt--){
        eps1 =  prod%bt;
        if(eps1 < eps0){
            /* better solution found */
            b0 = bt; eps0 = eps1;
        }
    } 
    *a = prod/b0;
    *b = b0;
}    


/**
 * compute chunk parameters for layout of arrays on the disk
 *   ---- a very simple algorithm to be refined later ----
 * @param elem_size[in] Size of stored data
 * @param block1[in]    Estimated size of request in dimension 1
 * @param block2[in]    Estimated size of request in dimension 2
 * @param dim1[in]      Size of DRA in dimension 1
 * @param dim2[in]      Size of DRA in dimension 2
 * @param chunk1[out]   Data block size in dimension 1?
 * @param chunk2[out]   Data block size in dimension 2?
 */
void dai_chunking(Integer elem_size, Integer block1, Integer block2, 
        Integer dim1, Integer dim2, Integer *chunk1, Integer *chunk2)
{
    Integer patch_size;

    *chunk1 = *chunk2 =0; 
    if(block1 <= 0 && block2 <= 0){

        *chunk1 = dim1;
        *chunk2 = dim2;

    }else if(block1 <= 0){
        *chunk2 = block2;
        *chunk1 = PARIO_MAX(1, DRA_BUF_SIZE/(elem_size* (*chunk2)));
    }else if(block2 <= 0){
        *chunk1 = block1;
        *chunk2 = PARIO_MAX(1, DRA_BUF_SIZE/(elem_size* (*chunk1)));
    }else{
        *chunk1 = block1;
        *chunk2 = block2;
    }

    /* need to correct chunk size to fit chunk1 x chunk2 request in buffer*/
    patch_size = (*chunk1)* (*chunk2)*elem_size;

    if (patch_size > ((Integer)DRA_BUF_SIZE)){

        if( *chunk1 == 1) *chunk2  = DRA_BUF_SIZE/elem_size;
        else if( *chunk2 == 1) *chunk1  = DRA_BUF_SIZE/elem_size;
        else {
            double  ratio = ((double)patch_size)/((double)DRA_BUF_SIZE); 
            /* smaller chunk to be scaled first */
            if(*chunk1 < *chunk2){
                dai_correct_chunking(chunk2,chunk1,DRA_BUF_SIZE/elem_size,ratio);
            }else{
                dai_correct_chunking(chunk1,chunk2,DRA_BUF_SIZE/elem_size,ratio);
            }
        }
    }
#ifdef DEBUG
    printf("\n%d:CREATE chunk=(%d,%d) elem_size=%d req=(%d,%d) buf=%d\n",
            pnga_nodeid(),*chunk1, *chunk2, elem_size, block1, block2,
            DRA_DBL_BUF_SIZE); 
    fflush(stdout);
#endif
}


/**
 * get a new handle for disk array 
 */
Integer dai_get_handle(void)
{
    Integer dra_handle =-1, candidate = 0;

    do{
        if(!DRA[candidate].actv){ 
            dra_handle=candidate;
            DRA[candidate].actv =1;
        }
        candidate++;
    }while(candidate < _max_disk_array && dra_handle == -1);
    return(dra_handle);
}      


/**
 * release handle -- makes array inactive
 */
void dai_release_handle(Integer *handle)
{
    DRA[*handle+DRA_OFFSET].actv =0;
    *handle = 0;
}


/**
 * find offset in file for (ilo,ihi) element
 */
void dai_file_location(section_t ds_a, Off_t* offset)
{
    Integer row_blocks, handle=ds_a.handle+DRA_OFFSET, offelem, cur_ld, part_chunk1;

    if((ds_a.lo[0]-1)%DRA[handle].chunk[0])
        dai_error("dai_file_location: not alligned ??",ds_a.lo[0]);

    row_blocks  = (ds_a.lo[0]-1)/DRA[handle].chunk[0];/* # row blocks from top*/
    part_chunk1 = DRA[handle].dims[0]%DRA[handle].chunk[0];/*dim1 in part block*/
    cur_ld      = (row_blocks == DRA[handle].dims[0] / DRA[handle].chunk[0]) ? 
        part_chunk1: DRA[handle].chunk[0];

    /* compute offset (in elements) */

    if(INDEPFILES(ds_a.handle) || DRA[handle].numfiles > 1) {

        Integer   CR, R; 
        Integer   i, num_part_block = 0;
        Integer   ioprocs=dai_io_procs(ds_a.handle); 
        Integer   iome = dai_io_nodeid(ds_a.handle);

        sect_to_blockM(ds_a, &CR); 

        R    = (DRA[handle].dims[0] + DRA[handle].chunk[0]-1)/DRA[handle].chunk[0];
        for(i = R -1; i< CR; i+=R) if(i%ioprocs == iome)num_part_block++;

        if(!part_chunk1) part_chunk1=DRA[handle].chunk[0];
        offelem = ((CR/ioprocs - num_part_block)*DRA[handle].chunk[0] +
                num_part_block * part_chunk1 ) * DRA[handle].chunk[1];

        /* add offset within block */
        offelem += ((ds_a.lo[1]-1) %DRA[handle].chunk[1])*cur_ld; 
    } else {

        offelem = row_blocks  * DRA[handle].dims[1] * DRA[handle].chunk[0];
        offelem += (ds_a.lo[1]-1)*cur_ld;

    }

    *offset = offelem* dai_sizeofM(DRA[handle].type); 
}


/**
 * write aligned block of data from memory buffer to d_a
 */
void dai_put(
        section_t    ds_a,
        void         *buf,
        Integer      ld,
        io_request_t *id)
{
    Integer handle = ds_a.handle + DRA_OFFSET, elem;
    Off_t   offset;
    Size_t  bytes;

    /* find location in a file where data should be written */
    dai_file_location(ds_a, &offset);

    if((ds_a.hi[0] - ds_a.lo[0] + 1) != ld) dai_error("dai_put: bad ld",ld); 

    /* since everything is aligned, write data to disk */
    elem = (ds_a.hi[0] - ds_a.lo[0] + 1) * (ds_a.hi[1] - ds_a.lo[1] + 1);
    bytes= (Size_t) elem * dai_sizeofM(DRA[handle].type);
    if( ELIO_OK != elio_awrite(DRA[handle].fd, offset, buf, bytes, id ))
        dai_error("dai_put failed", ds_a.handle);
}


/**
 * write zero at EOF
 */
void dai_zero_eof(Integer d_a)
{
    Integer handle = d_a+DRA_OFFSET, nelem;
    char byte;
    Off_t offset;

    byte = (char)0;

    if(INDEPFILES(d_a) || DRA[handle].numfiles > 1) {

        Integer   CR=0, i=0, nblocks=0; 
        section_t ds_a;
        /* number of processors that do io */
        Integer   ioprocs=dai_io_procs(d_a); 
        /* node id of current process (if it does io) */
        Integer   iome = dai_io_nodeid(d_a);

        /* total number of blocks in the disk resident array */
        nblocks = ((DRA[handle].dims[0]
                    + DRA[handle].chunk[0]-1)/DRA[handle].chunk[0])
            * ((DRA[handle].dims[1]
                        + DRA[handle].chunk[1]-1)/DRA[handle].chunk[1]);
        fill_sectionM(ds_a, d_a, 0, 0, 0, 0); 

        /* search for the last block for each I/O processor */
        for(i = 0; i <ioprocs; i++){
            CR = nblocks - 1 -i;
            if(CR % ioprocs == iome) break;
        }
        if(CR<0) return; /* no blocks owned */

        block_to_sectM(&ds_a, CR); /* convert block number to section */
        dai_file_location(ds_a, &offset);
        nelem = (ds_a.hi[0] - ds_a.lo[0] +1)*(ds_a.hi[1] - ds_a.lo[1] +1); 
        offset += ((Off_t)nelem) * dai_sizeofM(DRA[handle].type);

#ifdef DEBUG
        printf("me=%d zeroing EOF (%d) at %ld bytes \n",iome,CR,offset);
#endif
    } else {

        nelem = DRA[handle].dims[0]*DRA[handle].dims[1];
        offset = ((Off_t)nelem) * dai_sizeofM(DRA[handle].type);
    }

    if(elio_write(DRA[handle].fd, offset-1, &byte, 1) != (Size_t)1)
        dai_error("dai_zero_eof: write error ",0);
}


#ifdef CLEAR_BUF
static void dai_clear_buffer()
{
    /* int i, j; */
    /*
       for (j = 0; j < DRA_NBUF; j++) 
       for (i=0;i<DRA_DBL_BUF_SIZE;i++)
       ((double*)_dra_buffer_state[j].buffer)[i]=0.;
       */
}
#endif


/**
 * read aligned block of data from d_a to memory buffer
 */
void dai_get(section_t ds_a, void *buf, Integer ld, io_request_t *id)
{
    Integer handle = ds_a.handle + DRA_OFFSET, elem, rc;
    Off_t   offset;
    Size_t  bytes;

    /* find location in a file where data should be read from */
    dai_file_location(ds_a, &offset);

#       ifdef CLEAR_BUF
    dai_clear_buffer();
#       endif

    if((ds_a.hi[0] - ds_a.lo[0] + 1) != ld) dai_error("dai_get: bad ld",ld); 
    /* since everything is aligned, read data from disk */
    elem = (ds_a.hi[0] - ds_a.lo[0] + 1) * (ds_a.hi[1] - ds_a.lo[1] + 1);
    bytes= (Size_t) elem * dai_sizeofM(DRA[handle].type);
    rc= elio_aread(DRA[handle].fd, offset, buf, bytes, id );
    if(rc !=  ELIO_OK) dai_error("dai_get failed", rc);
}


void dai_assign_request_handle(Integer* request)
{
    int      i;

    *request = -1;
    for(i=0;i<MAX_REQ;i++)if(Requests[i].num_pending==0){
        *request = i;
        break;
    }

    if(*request ==-1) 
        dai_error("DRA: number of pending I/O requests exceeded",MAX_REQ);

    Requests[*request].na=0;
    Requests[*request].nu=0;
    Requests[*request].call_id = *request;
    Requests[*request].num_pending = ON;
    Requests[*request].call_id = *request;
}


/**
 * OPEN AN ARRAY THAT EXISTS ON THE DISK
 *
 * @param filename[in]
 * @param mode[in]
 * @param d_a[out]
 */
Integer drai_open(char *filename, Integer *mode, Integer *d_a)
{
    Integer handle;
    int emode;

    pnga_sync();

    /*** Get next free DRA handle ***/
    if( (handle = dai_get_handle()) == -1)
        dai_error("dra_open: too many disk arrays ", _max_disk_array);
    *d_a = handle - DRA_OFFSET;

    DRA[handle].mode = (int)*mode;
    strncpy (DRA[handle].fname, filename,  DRA_MAX_FNAME);

    /*translate DRA mode into ELIO mode*/
    emode = dai_elio_mode((int)*mode);

    if(dai_read_param(DRA[handle].fname, *d_a))return((Integer)-1);

    DRA[handle].indep = dai_file_config(filename); /*check file configuration*/

    if(dai_io_manage(*d_a)){ 

        if (INDEPFILES(*d_a) || DRA[handle].numfiles > 1) {

            sprintf(dummy_fname,"%s.%ld",DRA[handle].fname,(long)dai_io_nodeid(*d_a));
            DRA[handle].fd = elio_open(dummy_fname,emode, ELIO_PRIVATE);

        }else{

            DRA[handle].fd = elio_open(DRA[handle].fname,emode, ELIO_SHARED);
        }

        if(DRA[handle].fd ==NULL)dai_error("dra_open failed (null)",
                pnga_nodeid());
        if(DRA[handle].fd->fd ==-1)dai_error("dra_open failed (-1)",
                pnga_nodeid());  
    }


#ifdef DEBUG
    printf("\n%d:OPEN chunking=(%d,%d) type=%d buf=%d\n",
            pnga_nodeid(),DRA[handle].chunk[0], DRA[handle].chunk[1], 
            DRA[handle].type, DRA_DBL_BUF_SIZE);
    fflush(stdout);
#endif

    pnga_sync();

    /* printf("FILE OPENED!!\n"); */
    return(ELIO_OK);
}


/**
 * CLOSE AN ARRAY AND SAVE IT ON THE DISK
 */
Integer FATR dra_close_(Integer* d_a) /* input:DRA handle*/ 
{
    Integer handle = *d_a+DRA_OFFSET;
    int rc;

    pnga_sync();

    dai_check_handleM(*d_a, "dra_close");
    if(dai_io_manage(*d_a)) if(ELIO_OK != (rc=elio_close(DRA[handle].fd)))
        dai_error("dra_close: close failed",rc);
    dai_release_handle(d_a); 

    pnga_sync();

    return(ELIO_OK);
}


/**
 * decompose [ilo:ihi, jlo:jhi] into aligned and unaligned DRA subsections
 *
 * 
 * section [ilo:ihi, jlo:jhi] is decomposed into a number of
 * 'aligned' and 'unaligned' (on chunk1/chunk2 bounday) subsections
 * depending on the layout of the 2D array on the disk;
 *
 * 'cover' subsections correspond to 'unaligned' subsections and
 * extend them to aligned on chunk1/chunk2 boundaries;
 *
 * disk I/O will be actually performed on 'aligned' and
 * 'cover' instead of 'unaligned' subsections
 */
void dai_decomp_section(
        section_t ds_a,
        Integer aligned[][2*MAXDIM], 
        int *na,
        Integer cover[][2*MAXDIM],
        Integer unaligned[][2*MAXDIM], 
        int *nu) 
{
    Integer a=0, u=0, handle = ds_a.handle+DRA_OFFSET, off, chunk_units, algn_flag;


    aligned[a][ ILO ] = ds_a.lo[0]; aligned[a][ IHI ] = ds_a.hi[0];
    aligned[a][ JLO ] = ds_a.lo[1]; aligned[a][ JHI ] = ds_a.hi[1];

    switch   (DRA[handle].layout){
        case COLUMN : /* need to check row alignment only */

            algn_flag = ON; /* has at least one aligned subsection */

            /* top of section */
            off = (ds_a.lo[0] -1) % DRA[handle].chunk[0]; 
            if(off){ 

                if(MAX_UNLG<= u) 
                    dai_error("dai_decomp_sect:insufficient nu",u);

                chunk_units = (ds_a.lo[0] -1) / DRA[handle].chunk[0];

                cover[u][ ILO ] = chunk_units*DRA[handle].chunk[0] + 1;
                cover[u][ IHI ] = PARIO_MIN(cover[u][ ILO ] + DRA[handle].chunk[0]-1,
                        DRA[handle].dims[0]);

                unaligned[u][ ILO ] = ds_a.lo[0];
                unaligned[u][ IHI ] = PARIO_MIN(ds_a.hi[0],cover[u][ IHI ]);
                unaligned[u][ JLO ] = cover[u][ JLO ] = ds_a.lo[1];
                unaligned[u][ JHI ] = cover[u][ JHI ] = ds_a.hi[1];

                if(cover[u][ IHI ] < ds_a.hi[0]){
                    /* cover subsection ends above ihi */
                    if(MAX_ALGN<=a)
                        dai_error("dai_decomp_sect: na too small",a);
                    aligned[a][ ILO ] = cover[u][ IHI ]+1; 
                }else{
                    /* cover subsection includes ihi */
                    algn_flag = OFF;
                }
                u++;
            }

            /* bottom of section */
            off = ds_a.hi[0] % DRA[handle].chunk[0]; 
            if(off && (ds_a.hi[0] != DRA[handle].dims[0]) && (algn_flag == ON)){

                if(MAX_UNLG<=u) 
                    dai_error("dai_decomp_sect:insufficient nu",u); 
                chunk_units = ds_a.hi[0] / DRA[handle].chunk[0];

                cover[u][ ILO ] = chunk_units*DRA[handle].chunk[0] + 1;
                cover[u][ IHI ] = PARIO_MIN(cover[u][ ILO ] + DRA[handle].chunk[0]-1,
                        DRA[handle].dims[0]);

                unaligned[u][ ILO ] = cover[u][ ILO ];
                unaligned[u][ IHI ] = ds_a.hi[0];
                unaligned[u][ JLO ] = cover[u][ JLO ] = ds_a.lo[1];
                unaligned[u][ JHI ] = cover[u][ JHI ] = ds_a.hi[1];

                aligned[a][ IHI ] = PARIO_MAX(1,unaligned[u][ ILO ]-1);
                algn_flag=(DRA[handle].chunk[0] == DRA[handle].dims[0])?OFF:ON;
                u++;
            }
            *nu = (int)u;
            if(aligned[0][ IHI ]-aligned[0][ ILO ] < 0) algn_flag= OFF;
            *na = (algn_flag== OFF)? 0: 1;
            break;

        case ROW : /* we need to check column alignment only */

        default: dai_error("dai_decomp_sect: ROW layout not yet implemented",
                         DRA[handle].layout);
    }
}


/**
 * given current (i,j) compute (ni, nj) - next loop index
 * o - outermost loop, i- innermost loop
 * iinc increment for i
 * oinc increment for o
 */
int dai_next2d(Integer* i, Integer imin, Integer imax, Integer iinc, 
        Integer* o, Integer omin, Integer omax, Integer oinc)
{
    int retval;
    if (*o == 0  || *i == 0) {
        /* to handle initial out-of range indices */
        *o = omin;
        *i = imin;
    } else {
        *i = *i + iinc;
    }
    if (*i > imax) {
        *i = imin;
        *o += oinc;
    }
    retval = (*o <= omax);
    return retval;
}


/**
 * compute next chunk of array to process
 */
int dai_next_chunk(Integer req, Integer* list, section_t* ds_chunk)
{
    Integer   handle = ds_chunk->handle+DRA_OFFSET;
    int       retval;

    if(INDEPFILES(ds_chunk->handle) || DRA[handle].numfiles > 1)
        if(ds_chunk->lo[1] && DRA[handle].chunk[1]>1) 
            ds_chunk->lo[1] -= (ds_chunk->lo[1] -1) % DRA[handle].chunk[1];

    retval = dai_next2d(&ds_chunk->lo[0], list[ ILO ], list[ IHI ],
            DRA[handle].chunk[0],
            &ds_chunk->lo[1], list[ JLO ], list[ JHI ],
            DRA[handle].chunk[1]);
    if(!retval) return(retval);

    ds_chunk->hi[0] = PARIO_MIN(list[ IHI ], ds_chunk->lo[0] + DRA[handle].chunk[0] -1);
    ds_chunk->hi[1] = PARIO_MIN(list[ JHI ], ds_chunk->lo[1] + DRA[handle].chunk[1] -1);

    if(INDEPFILES(ds_chunk->handle) || DRA[handle].numfiles > 1) { 
        Integer jhi_temp =  ds_chunk->lo[1] + DRA[handle].chunk[1] -1;
        jhi_temp -= jhi_temp % DRA[handle].chunk[1];
        ds_chunk->hi[1] = PARIO_MIN(ds_chunk->hi[1], jhi_temp); 

        /*this line was absent from older version on bonnie that worked */
        if(ds_chunk->lo[1] < list[ JLO ]) ds_chunk->lo[1] = list[ JLO ]; 
    }

    return 1;
}


int dai_myturn(section_t ds_chunk)
{
    /* Integer   handle = ds_chunk.handle+DRA_OFFSET; */
    Integer   ioprocs = dai_io_procs(ds_chunk.handle);
    Integer   iome    = dai_io_nodeid(ds_chunk.handle);

    /*    if(INDEPFILES(ds_chunk.handle) || DRA[handle].numfiles > 1){ */

    /* compute cardinal number for the current chunk */
    nsect_to_blockM(ds_chunk, &_dra_turn);

    /*    }else{
          _dra_turn++;
          } */

    return ((_dra_turn%ioprocs) == iome);
}


/**
 * print routine for debugging purposes only (double)
 */
void dai_print_buf(double *buf, Integer ld, Integer rows, Integer cols)
{
    int i,j;
    printf("\n ld=%ld rows=%ld cols=%ld\n",(long)ld,(long)rows,(long)cols);

    for (i=0; i<rows; i++){
        for (j=0; j<cols; j++)
            printf("%f ", buf[j*ld+i]);
        printf("\n");
    }
}


void dra_set_mode_(Integer* val)
{
}



void ga_move_1d(int op, section_t gs_a, section_t ds_a,
        section_t ds_chunk, void* buffer, Integer ldb)
{
    Integer index, ldd = gs_a.hi[0] - gs_a.lo[0] + 1, one=1;
    Integer atype, cols, rows, elemsize, ilo, ihi;
    Integer istart, iend, jstart, jend;
    void  (*f)(Integer,Integer*,Integer*,void*,Integer*); 
    char *buf = (char*)buffer;

    if(op==LOAD) f = pnga_get;
    else f = pnga_put;

    pnga_inquire(gs_a.handle, &atype, &rows, &cols);     
    elemsize = MA_sizeof(atype, 1, MT_C_CHAR);

    /* find where in global array the first dra chunk element in buffer goes*/
    index = ds_chunk.lo[0] - ds_a.lo[0];
    istart = index%ldd + gs_a.lo[0]; 
    jstart = index/ldd + gs_a.lo[1];

    /* find where in global array the last dra chunk element in buffer goes*/
    index = ds_chunk.hi[0] - ds_a.lo[0];
    iend = index%ldd + gs_a.lo[0]; 
    jend = index/ldd + gs_a.lo[1];

    /* we have up to 3 rectangle chunks corresponding to gs_chunk 
       .|' incomplete first column, full complete middle column, and
       incomplete last column */
    if(istart != gs_a.lo[0] || jstart==jend ){
        Integer lo[2], hi[2];
        ilo = istart; 
        ihi = gs_a.hi[0]; 
        if(jstart==jend) ihi=iend;
        lo[0] = ilo; lo[1] = jstart;
        hi[0] = ihi; hi[1] = jstart;
        f(gs_a.handle, lo, hi, buf, &one); 
        buf += elemsize*(ihi -ilo+1);
        if(jstart==jend) return;
        jstart++;
    }

    if(iend != gs_a.hi[0]) jend--;

    if(jstart <= jend) { 
        Integer lo[2], hi[2];
        lo[0] = gs_a.lo[0]; lo[1] = jstart;
        hi[0] = gs_a.hi[0]; hi[1] = jend;
        f(gs_a.handle, lo, hi, buf, &ldd);
        buf += elemsize*ldd*(jend-jstart+1); 
    } 

    if(iend != gs_a.hi[0]){
        Integer lo[2], hi[2];
        jend++; /* Since decremented above */  
        lo[0] = gs_a.lo[0]; lo[1] = jend;
        hi[0] = iend;       hi[1] = jend;
        f(gs_a.handle, lo, hi, buf, &one);
    }
}


void ga_move(int op, int trans, section_t gs_a, section_t ds_a, 
        section_t ds_chunk, void* buffer, Integer ldb)
{
    if(!trans && (gs_a.lo[0]- gs_a.hi[0] ==  ds_a.lo[0]- ds_a.hi[0]) ){
        /*** straight copy possible if there's no reshaping or transpose ***/

        /* determine gs_chunk corresponding to ds_chunk */
        section_t gs_chunk = gs_a;
        dai_dest_indicesM(ds_chunk.lo[0], ds_chunk.lo[1], ds_a.lo[0], ds_a.lo[1], 
                ds_a.hi[0]-ds_a.lo[0]+1, &gs_chunk.lo[0], &gs_chunk.lo[1], 
                gs_a.lo[0], gs_a.lo[1],   gs_a.hi[0] - gs_a.lo[0] + 1);
        dai_dest_indicesM(ds_chunk.hi[0], ds_chunk.hi[1], ds_a.lo[0], ds_a.lo[1], 
                ds_a.hi[0]-ds_a.lo[0]+1, &gs_chunk.hi[0], &gs_chunk.hi[1],
                gs_a.lo[0], gs_a.lo[1],  gs_a.hi[0] - gs_a.lo[0] + 1);

        /* move data */
        if(op==LOAD) ga_get_sectM(gs_chunk, buffer, ldb);
        else         ga_put_sectM(gs_chunk, buffer, ldb);

#ifdef MOVE1D_ENABLED
    }else if(!trans && (ds_a.lo[1]==ds_a.hi[1]) ){

        /* for a 1-dim section (column) some optimization possible */
        ga_move_1d(op, gs_a, ds_a, ds_chunk, buffer, ldb);        
#endif
    }else{
        /** due to generality of this transformation scatter/gather is required **/

        MA_AccessIndex iindex, jindex, vindex, pindex;
        Integer ihandle, jhandle, vhandle, phandle;
        int type = DRA[ds_a.handle+DRA_OFFSET].type;
        Integer i, j, ii, jj, base,nelem;  
        char    *base_addr;

        if(pnga_nodeid()==0) printf("DRA warning: using scatter/gather\n");

        nelem = (ds_chunk.hi[0]-ds_chunk.lo[0]+1)
            * (ds_chunk.hi[1]-ds_chunk.lo[1]+1);
        if(!MA_push_get(C_INT, nelem, "i_", &ihandle, &iindex))
            dai_error("DRA move: MA failed-i ", 0L);
        if(!MA_push_get(C_INT, nelem, "j_", &jhandle, &jindex))
            dai_error("DRA move: MA failed-j ", 0L);
        if(!MA_push_get(type, nelem, "v_", &vhandle, &vindex))
            dai_error("DRA move: MA failed-v ", 0L);
        if(!MA_get_pointer(vhandle, &base_addr))
            dai_error("DRA move: MA get_pointer failed ", 0L);

        if(trans==TRANS) 
            ITERATOR_2D(i,j, base, ds_chunk) {
                dai_dest_indicesM(j, i, ds_a.lo[0], ds_a.lo[1], ds_a.hi[0]-ds_a.lo[0]+1, 
                        INT_MB+base+iindex, INT_MB+base+jindex,
                        gs_a.lo[0], gs_a.lo[1],  gs_a.hi[0] - gs_a.lo[0] + 1);
            }
        else
            ITERATOR_2D(i,j, base, ds_chunk) {
                dai_dest_indicesM(i, j, ds_a.lo[0], ds_a.lo[1], ds_a.hi[0]-ds_a.lo[0]+1, 
                        INT_MB+base+iindex, INT_MB+base+jindex,
                        gs_a.lo[0], gs_a.lo[1],  gs_a.hi[0] - gs_a.lo[0] + 1);
            }

        /* move data */
        if(op==LOAD){

            if(!MA_push_get(C_INT, nelem, "pindex", &phandle, &pindex))
                dai_error("DRA move: MA failed-p ", 0L);
            for(i=0; i< nelem; i++) INT_MB[pindex+i] = i; 
            pnga_gather2d(gs_a.handle, base_addr, INT_MB+iindex, INT_MB+jindex, nelem);
            COPY_TYPE(GATHER, type, ds_chunk);
            MA_pop_stack(phandle);

        }else{ 

            COPY_TYPE(SCATTER, type, ds_chunk);
            pnga_scatter2d(gs_a.handle, base_addr, INT_MB+iindex, INT_MB+jindex, nelem);
        }

        MA_pop_stack(vhandle);
        MA_pop_stack(jhandle);
        MA_pop_stack(ihandle);
    }
}


/**
 * @param op[in]       flag for read or write
 * @param trans[in]    flag for transpose
 * @param gs_a[in]     complete section of global array
 * @param ds_a[in]     complete section of DRA
 * @param ds_chunk[in] actual DRA chunk
 * @param buffer[in]   pointer to io buffer containing DRA cover section
 * @param ldb
 * @param ga_movhdl
 */
void nga_move(int op, int trans, section_t gs_a, section_t ds_a,
        section_t ds_chunk, void* buffer, Integer ldb[], Integer *ga_movhdl)
{
    Integer ndim = gs_a.ndim, i;
    logical consistent = TRUE;

#if WALLTIME
    double ss0,tt0;
    walltime_(&ss0,&tt0);
    printf("p[%d] Beginning nga_move: %16.6f\n",pnga_nodeid(),tt0);
#endif
    if (!trans) {
        for (i=0; i<ndim-1; i++) 
            if (gs_a.lo[i]-gs_a.hi[i] != ds_a.lo[i]-ds_a.hi[i]) consistent = FALSE;
    } else {
        for (i=0; i<ndim-1; i++) 
            if (gs_a.lo[ndim-1-i]-gs_a.hi[ndim-1-i]
                    != ds_a.lo[i]-ds_a.hi[i]) consistent = FALSE;
    }
    if (!trans && consistent){

        /*** straight copy possible if there's no reshaping or transpose ***/

        /* determine gs_chunk corresponding to ds_chunk */
        section_t gs_chunk = gs_a;
        ndai_dest_indicesM(ds_chunk, ds_a, gs_chunk, gs_a);
        consistent = TRUE;
        for (i=0; i<ndim; i++) {
            if (gs_chunk.hi[i]<gs_chunk.lo[i] || gs_chunk.lo[i]<0) {
                consistent = FALSE;
            }
        }
        if (!consistent) {
            for(i=0; i<ndim; i++) {
                printf("gs_chunk[%d] %5d:%5d  ds_chunk[%d] %5d:%5d",
                        (int)i,(int)gs_chunk.lo[i],(int)gs_chunk.hi[i],
                        (int)i,(int)ds_chunk.lo[i],(int)ds_chunk.hi[i]);
                printf(" gs_a[%d] %5d:%5d  ds_a[%d] %5d:%5d\n",
                        (int)i,(int)gs_a.lo[i],(int)gs_a.hi[i],
                        (int)i,(int)ds_a.lo[i],(int)ds_a.hi[i]);
            }
        }
        /* move data */
        if (op==LOAD) {
            nga_get_sectM(gs_chunk, buffer, ldb, ga_movhdl);
        } else {
            nga_put_sectM(gs_chunk, buffer, ldb, ga_movhdl);
        }

    } else if (trans && consistent) {
        /* Only transpose is supported, not reshaping, so scatter/gather is not
           required */
        MA_AccessIndex vindex;
        Integer vhandle, index[MAXDIM];
        Integer i, j, itmp, jtmp, nelem, ldt[MAXDIM], ldg[MAXDIM];
        Integer nelem1, nelem2, nelem3;
        int type = DRA[ds_a.handle+DRA_OFFSET].type;
        char    *base_addr;
        section_t gs_chunk = gs_a;

        /* create space to copy transpose of DRA section into */
        nelem = 1;
        for (i=0; i<ndim; i++) {
            nelem *= (ds_chunk.hi[i] - ds_chunk.lo[i] + 1);
        }
        nelem1 = 1;
        ndai_trnsp_dest_indicesM(ds_chunk, ds_a, gs_chunk, gs_a);
        for (i=1; i<ndim; i++) nelem1 *= (gs_chunk.hi[i] - gs_chunk.lo[i] + 1);
        nelem2 = 1;
        for (i=0; i<ndim-1; i++) nelem2 *= ldb[i];
        if(!MA_push_get(type, nelem, "v_", &vhandle, &vindex))
            dai_error("DRA move: MA failed-v ", 0L);
        if(!MA_get_pointer(vhandle, &base_addr))
            dai_error("DRA move: MA get_pointer failed ", 0L);

        /* copy and transpose relevant numbers from IO buffer to temporary array */
        for (i=1; i<ndim; i++) ldt[ndim-1-i] = ds_chunk.hi[i] - ds_chunk.lo[i] + 1;
        if (op == LOAD) {
            /* transpose buffer with data from global array */
            for (i=0; i<ndim; i++) ldg[i] = gs_chunk.hi[i] - gs_chunk.lo[i] + 1;
            /* copy data from global array to temporary buffer */
            nga_get_sectM(gs_chunk, base_addr, ldg, ga_movhdl); 
            for (i=0; i<nelem1; i++ ) {
                /* find indices of elements in MA buffer */
                if (ndim > 1) {
                    itmp = i;
                    index[1] = itmp%ldg[1];
                    for (j=2; j<ndim; j++) {
                        itmp = (itmp-index[j-1])/ldg[j-1];
                        if (j != ndim-1) {
                            index[j] = itmp%ldg[j];
                        } else {
                            index[j] = itmp;
                        }
                    }
                    nelem3 = index[1];
                    for (j=2; j<ndim; j++) {
                        nelem3 *= ldb[ndim-1-j];
                        nelem3 += index[j];
                    }
                } else {
                    nelem2 = 1;
                    nelem3 = 0;
                }
                /* find corresponding indices of element from IO buffer */
                itmp = ldg[0]*i;
                jtmp = nelem3;
                for (j=0; j<ldg[0]; j++) {
                    switch(type){
                        case C_DBL:
                            ((double*)buffer)[jtmp] = ((double*)base_addr)[itmp];
                            break;
                        case C_INT:
                            ((int*)buffer)[jtmp] = ((int*)base_addr)[itmp];
                            break;
                        case C_LONG:
                            ((long*)buffer)[jtmp] = ((long*)base_addr)[itmp];
                            break;
                        case C_DCPL:
                            ((double*)buffer)[2*jtmp] = ((double*)base_addr)[2*itmp];
                            ((double*)buffer)[2*jtmp+1] = ((double*)base_addr)[2*itmp+1];
                            break;
                        case C_SCPL:
                            ((float*)buffer)[2*jtmp] = ((float*)base_addr)[2*itmp];
                            ((float*)buffer)[2*jtmp+1] = ((float*)base_addr)[2*itmp+1];
                            break;
                        case C_FLOAT:
                            ((float*)buffer)[jtmp] = ((float*)base_addr)[itmp];
                            break;
                    }
                    itmp++;
                    jtmp += nelem2;
                }
            }
        } else {
            /* get transposed indices */
            for (i=0; i<ndim; i++) ldg[i] = gs_chunk.hi[i] - gs_chunk.lo[i] + 1;
            for (i=0; i<nelem1; i++ ) {
                /* find indices of elements in MA buffer */
                if (ndim > 1) {
                    itmp = i;
                    index[1] = itmp%ldg[1];
                    for (j=2; j<ndim; j++) {
                        itmp = (itmp-index[j-1])/ldg[j-1];
                        if (j != ndim-1) {
                            index[j] = itmp%ldg[j];
                        } else {
                            index[j] = itmp;
                        }
                    }
                    nelem3 = index[1];
                    for (j=2; j<ndim; j++) {
                        nelem3 *= ldb[ndim-1-j];
                        nelem3 += index[j];
                    }
                } else {
                    nelem2 = 1;
                    nelem3 = 0;
                }
                /* find corresponding indices of element from IO buffer */
                itmp = ldg[0]*i;
                jtmp = nelem3;
                for (j=0; j<ldg[0]; j++) {
                    switch(type){
                        case C_DBL:
                            ((double*)base_addr)[itmp] = ((double*)buffer)[jtmp];
                            break;
                        case C_INT:
                            ((int*)base_addr)[itmp] = ((int*)buffer)[jtmp];
                            break;
                        case C_LONG:
                            ((long*)base_addr)[itmp] = ((long*)buffer)[jtmp];
                            break;
                        case C_DCPL:
                            ((double*)base_addr)[2*itmp] = ((double*)buffer)[2*jtmp];
                            ((double*)base_addr)[2*itmp+1] = ((double*)buffer)[2*jtmp+1];
                            break;
                        case C_SCPL:
                            ((float*)base_addr)[2*itmp] = ((float*)buffer)[2*jtmp];
                            ((float*)base_addr)[2*itmp+1] = ((float*)buffer)[2*jtmp+1];
                            break;
                        case C_FLOAT:
                            ((float*)base_addr)[itmp] = ((float*)buffer)[jtmp];
                            break;
                    }
                    itmp++;
                    jtmp += nelem2;
                }
            }
            nga_put_sectM(gs_chunk, base_addr, ldt, ga_movhdl); 
        }
        MA_pop_stack(vhandle);
    } else {
        dai_error("DRA move: Inconsistent dimensions found ", 0L);
    }
#if WALLTIME
    walltime_(&ss0,&tt0);
    printf("p[%d] Ending nga_move: %16.6f\n",pnga_nodeid(),tt0);
#endif
}


/**
 * executes callback function associated with completion of asynch. I/\O
 */
void dai_exec_callback(char *buf, int caller)
{
    args_t   *arg;
    char *buffer;
    buf_info *bi;

#if WALLTIME
    double ss0,tt0;
#endif

    bi = (buf_info*) buf;
    if(bi->callback==OFF)
        return;

    bi->callback = OFF;

    arg = &(bi->args);
    /* bail if there is no valid global array handle */
    if (arg->gs_a.handle == 0)
        return;

    buffer = (char*) (buf + sizeof(buf_info));
    if (caller == WAIT) {/* call blocking nga_move() */
        nga_move(arg->op, arg->transp, arg->gs_a, arg->ds_a, arg->ds_chunk, buffer, arg->ld, NULL);
        free_buf(&buf_ctxt, buf);
    }
    else if (caller == PROBE) /* call non-blocking nga_move() */
        nga_move(arg->op, arg->transp, arg->gs_a, arg->ds_a, arg->ds_chunk, buffer, arg->ld, &(bi->ga_movhdl));

}


/**
 * wait until buffer space associated with request is avilable
 */
void dai_wait(Integer req0)
{
    /*
    Integer req;
    int ibuf;
    */

    /* wait for all requests to complete on buffer Requests[req].ibuf */

    /*        ibuf = Requests[req0].ibuf;
              for(req=0; req<MAX_REQ; req++)
              if (Requests[req].num_pending && Requests[req].ibuf == ibuf)
              if (elio_wait(&_dra_buffer_state[_dra_cur_buf].id)==ELIO_OK)
              dai_exec_callback(Requests + req);
              else
              dai_error("dai_wait: DRA internal error",0);
              */
}


/**
 * WAIT FOR COMPLETION OF DRA OPERATION ASSOCIATED WITH request
 */ 
Integer FATR dra_wait_(Integer* request)
{
#if WALLTIME
    double ss0, tt0;
#endif
    if(*request == DRA_REQ_INVALID) return(ELIO_OK);
#if WALLTIME
    walltime(&ss0,&tt0);
    printf("p[%d] executing dra_wait: %16.6f\n",pnga_nodeid(),tt0);
#endif

    /* complete all outstanding operations and release the corresponding buffers invloved with this request */

    buf_complete_call(&buf_ctxt, Requests[*request].call_id);

    /* mark this request to be no longer pending */
    Requests[*request].num_pending=0;

    pnga_sync();

    return(ELIO_OK);

}


/**
 * TEST FOR COMPLETION OF DRA OPERATION ASSOCIATED WITH request
 */
Integer FATR dra_probe_(
        Integer *request, /* [in] */
        Integer *status)  /* [out] */
{
    Integer done;
    int  stat, i, k, call_id, op_code, n_buf, ret;
    io_request_t *io_req;
    Integer *ga_movhdl;
    char *op = "*", *bufs[MAXBUF];
    int cb[MAXBUF]; /* will mark id of bufs for which we'll call a callback */
    buf_info *bi;

    k = 1; 
    ret = ELIO_OK;
    for (i = 0; i < MAXBUF; i++)
    {
        cb[i] = 0;
        bufs[i] = NULL;
    }

    if(*request == DRA_REQ_INVALID || Requests[*request].num_pending ==0) {
        *status = ELIO_DONE;
        ret = ELIO_OK;
    }

    call_id = Requests[*request].call_id;  
    /* get the buffers associated with this call_id */
    if (get_bufs_of_call_id(&buf_ctxt, call_id, &n_buf, bufs) != 0)
        ret = ELIO_OK;

    for (i = 0; i < n_buf; i++) {
        bi = (buf_info*) bufs[i]; 

        op_code = bi->op;
        io_req = &(bi->io_req);
        ga_movhdl = &(bi->ga_movhdl);

        if (op_code == DRA_OP_WRITE) {
            /* last op is a disk write */
            if(elio_probe(io_req, &stat) != ELIO_OK) {
                ret = DRA_FAIL;
                k = 0;
                break;
            }
            if (stat != ELIO_DONE) {
                k = 0;
            }
            else {
                free_buf(&buf_ctxt, bufs[i]);
            }
        }
        else if (op_code == DRA_OP_READ) {
            /* last op depends on aligned or unaligned transfer */
            if (bi->align == 0) { /* unaligned read */
                /* last op is a ga move */
                if (NGA_NbTest(ga_movhdl) == 0) { /* ga op not complete */
                    k = 0;
                }
                else { /* ga op complete, free this buf */
                    free_buf(&buf_ctxt, bufs[i]);
                }
            }
            else { /* if aligned read, last op is a disk read */
                if(elio_probe(io_req, &stat) != ELIO_OK) {
                    ret = DRA_FAIL;
                    k = 0;
                    break;
                }

                if (stat != ELIO_DONE) 
                    k = 0;
                else { /* disk read done, initiate/test ga move */
                    /* callback=OFF means ga move done/being done */
                    if (bi->callback == OFF && NGA_NbTest(ga_movhdl) == 0)
                        k = 0;
                    else if (bi->callback == OFF && NGA_NbTest(ga_movhdl) ==1) {
                        free_buf(&buf_ctxt, bufs[i]);
                    }
                    else if (bi->callback == ON) {/* need to call callback */
                        k = 0;
                        cb[i] = 1; /* mark for a ga move */
                    }
                }
            }
        }
    }

    done = (Integer) k;

    /* determine global status */
    pnga_gop(pnga_type_f2c(MT_F_INT), &done, (Integer)1, op);

    if(done){
        *status = ELIO_DONE;
        Requests[*request].num_pending = 0;
    }
    else {
        *status = 0;
        for (i = 0; i < n_buf; i++)
            if (cb[i]) 
                dai_exec_callback(bufs[i], PROBE);
    }
    if (ret == DRA_FAIL)
        *status = 0; /* basically value of status is irrelevant/undetermined in this case */
    return ((Integer) ret);
}


/**
 * Returns control to DRA for a VERY short time to improve progress
 */
void dra_flick_()
{
    Integer req, stat;

    for (req = 0; req < MAX_REQ; req++) {
        if (Requests[req].num_pending) {
            dra_probe_(&req, &stat);
        }
    }
}


/**
 * INQUIRE PARAMETERS OF EXISTING DISK ARRAY
 * @param d_a[in] DRA handle
 * @param type[out]
 * @param dim1[out]
 * @param dim2[out]
 * @param name[out]
 * @param filename[out]
 */
Integer drai_inquire(Integer *d_a, Integer *type, Integer *dim1, Integer *dim2,
        char *name, char *filename)
{
    Integer handle=*d_a+DRA_OFFSET;

    dai_check_handleM(*d_a,"dra_inquire");

    *type = (Integer)DRA[handle].type;
    *dim1 = DRA[handle].dims[0];
    *dim2 = DRA[handle].dims[1];
    strcpy(name, DRA[handle].name);
    strcpy(filename, DRA[handle].fname);

    return(ELIO_OK);
}


/**
 * DELETE DISK ARRAY -- relevant file(s) gone
 *
 * @param d_a[in] DRA handle
 */
Integer FATR dra_delete_(Integer* d_a)
{
    Integer handle = *d_a+DRA_OFFSET;
    int rc;

    pnga_sync();

    dai_check_handleM(*d_a,"dra_delete");
    dai_delete_param(DRA[handle].fname,*d_a);

    if(dai_io_manage(*d_a)) if(ELIO_OK != (rc=elio_close(DRA[handle].fd)))
        dai_error("dra_close: close failed",rc);

    if(dai_file_master(*d_a)) {
        if(INDEPFILES(*d_a) || DRA[handle].numfiles > 1){ 
            sprintf(dummy_fname,"%s.%ld",DRA[handle].fname,(long)dai_io_nodeid(*d_a));
            elio_delete(dummy_fname);
        } else {
            elio_delete(DRA[handle].fname);
        }
    }

    dai_release_handle(d_a); 

    pnga_sync();
    return(ELIO_OK);
}


/**
 * TERMINATE DRA DRATA STRUCTURES
 */
Integer FATR dra_terminate_()
{
    free(DRA);
    buf_terminate(&buf_ctxt);

    pnga_sync();
    return(ELIO_OK);
}


/**
 * compute chunk parameters for layout of arrays on the disk
 *   ---- a very simple algorithm to be refined later ----
 *
 * @param elem_size[in]  Size of individual data element in bytes
 * @param ndim[in]       Dimension of DRA
 * @param block_orig[in] Estimated size of request in each coordinate
 *                       direction. If size is unknown then use -1.
 * @param dims[in]       Size of DRA in each coordinate direction
 * @param chunk[out]     Size of data block size (in elements) in each
 *                       coordinate direction
 */
void ndai_chunking(Integer elem_size, Integer ndim, Integer block_orig[], 
        Integer dims[], Integer chunk[])
{
    long patch_size, tmp_patch;
    Integer i, j, block[MAXDIM], block_map[MAXDIM];
    double ratio;
    logical full_buf, some_neg, overfull_buf;
    /* copy block_orig so that original guesses are not destroyed */
    for (i=0; i<ndim; i++) block[i] = block_orig[i];

    /* do some preliminary checks on block to make sure initial guesses
       are less than corresponding DRA dimensions */
    for (i=0; i<ndim; i++) {
        if (block[i] > dims[i]) block[i] = dims[i];
    }
    /* do additional adjustments to see if initial guesses are near some
       perfect factors of DRA dimensions */
    for (i=0; i<ndim; i++) {
        if (block[i] > 0 && block[i]<dims[i]) {
            if (dims[i]%block[i] != 0) {
                ratio = (double)dims[i]/(double)block[i];
                j = (int)(ratio+0.5);
                if (dims[i]%j ==0) block[i] = dims[i]/j;
            }
        }
    }

    /* initialize chunk array to zero and find out how big patch is based
       on specified block dimensions */
    patch_size = 1;
    some_neg = FALSE;
    full_buf = FALSE;
    overfull_buf = FALSE;
    for (i=0; i<ndim; i++) {
        if (block[i] > 0) patch_size *= (long)block[i];
        else some_neg = TRUE;
    }
    if (patch_size*((long)elem_size) > ((long)DRA_BUF_SIZE))
        overfull_buf = TRUE;

    /* map dimension sizes from highest to lowest */
    block_sortM(ndim, dims, block_map);

    /* IO buffer is not full and there are some unspecied chunk dimensions.
       Set unspecified dimensions equal to block dimensions until buffer
       is filled. */
    if (!full_buf && !overfull_buf && some_neg) {
        for (i=ndim-1; i>=0; i--) {
            if (block[block_map[i]] < 0) {
                tmp_patch = patch_size * ((long)dims[block_map[i]]);
                if (tmp_patch*elem_size < ((long)DRA_BUF_SIZE)) {
                    patch_size *= (long)dims[block_map[i]];
                    block[block_map[i]] = dims[block_map[i]];
                } else {
                    block[block_map[i]] = (Integer)(DRA_BUF_SIZE/(patch_size*((long)elem_size))); 
                    patch_size *= ((long)block[block_map[i]]);
                    full_buf = TRUE;
                }
            }
        }
    }

    /* copy block array to chunk array */
    for (i=0; i<ndim; i++) {
        if (block[i] > 0) chunk[i] = block[i];
        else chunk[i] = 1;
    }

    /* If patch overfills buffer, scale patch down until it fits */
    if (overfull_buf) {
        ratio = ((double)DRA_BUF_SIZE)
            / ((double)(patch_size*((long)elem_size)));
        ratio = pow(ratio,1.0/((double)ndim));
        patch_size = 1;
        for (i=0; i<ndim; i++) {
            chunk[i] = (int)(((double)chunk[i])*ratio);
            if (chunk[i] < 1) chunk[i] = 1;
            patch_size *= ((long)chunk[i]);
        }
    }

#ifdef DEBUG
    printf("Current patch at 2 is %d\n",(int)patch_size*elem_size);
#endif
    /* set remaining block sizes equal to 1 */
    for (i=0; i<ndim; i++) {
        if (chunk[i] == 0) chunk[i] = 1;
    }
    /* Patch size may be slightly larger than buffer. If so, nudge
       size down until patch is smaller than buffer. */
    if (((long)elem_size)*patch_size > ((long)DRA_BUF_SIZE)) {
        /* map chunks from highest to lowest */
        block_sortM(ndim, chunk, block_map);
        for (i=0; i < ndim; i++) {
            while (chunk[block_map[i]] > 1 &&
                    ((long)elem_size)*patch_size > ((long)DRA_BUF_SIZE)) {
                patch_size /= ((long)chunk[block_map[i]]);
                chunk[block_map[i]]--;
                patch_size *= ((long)chunk[block_map[i]]);
            }
        }
    }
}


/**
 * find offset in file for (lo,hi) element
 */
void ndai_file_location(section_t ds_a, Off_t* offset)
{
    Integer handle=ds_a.handle+DRA_OFFSET, ndim, i, j;
    Integer blocks[MAXDIM], part_chunk[MAXDIM], cur_ld[MAXDIM];
    long par_block[MAXDIM];
    long offelem=0;


    ndim = DRA[handle].ndim;
    for (i=0; i<ndim-1; i++) {
        if((ds_a.lo[i]-1)%DRA[handle].chunk[i])
            dai_error("ndai_file_location: not alligned ??",ds_a.lo[i]);
    }

    for (i=0; i<ndim; i++) {
        /* number of blocks from edge */
        blocks[i] = (ds_a.lo[i]-1)/DRA[handle].chunk[i];
        /* size of incomplete chunk */
        part_chunk[i] = DRA[handle].dims[i]%DRA[handle].chunk[i];
        /* stride for this block of data in this direction */
        cur_ld[i] = (blocks[i] == DRA[handle].dims[i]/DRA[handle].chunk[i]) ?
            part_chunk[i]: DRA[handle].chunk[i];
    }

    /* compute offset (in elements) */

    if (INDEPFILES(ds_a.handle) || DRA[handle].numfiles > 1) {
        Integer   CR, block_dims[MAXDIM]; 
        Integer   index[MAXDIM];
        long      nelem;
        Integer   i, j;
        Integer   ioprocs = dai_io_procs(ds_a.handle); 
        Integer   iome = dai_io_nodeid(ds_a.handle);

        /* Find index of current block and find number of chunks in
           each dimension of DRA */
        nsect_to_blockM(ds_a, &CR); 
        for (i=0; i<ndim; i++) {
            block_dims[i] = (DRA[handle].dims[i]+DRA[handle].chunk[i]-1)
                / DRA[handle].chunk[i];
        }
        if (iome >= 0) {
            offelem = 0;
            for (i=iome; i<CR; i+=ioprocs) {
                /* Copy i because macro destroys i */
                nblock_to_indicesM(index,ndim,block_dims,i);
                nelem = 1;
                for (j=0; j<ndim; j++) {
                    if (index[j]<block_dims[j]-1) {
                        nelem *= (long)DRA[handle].chunk[j];
                    } else {
                        if (part_chunk[j] != 0) {
                            nelem *= (long)part_chunk[j];
                        } else {
                            nelem *= (long)DRA[handle].chunk[j];
                        }
                    }
                }
                offelem += nelem;
            }
            /* add fractional offset for current block */
            nelem = 1;
            nblock_to_indicesM(index,ndim,block_dims,CR);
            for (i=0; i<ndim-1; i++) {
                if (index[i]<block_dims[i]-1) {
                    nelem *= (long)DRA[handle].chunk[i];
                } else {
                    if (part_chunk[i] != 0) {
                        nelem *= (long)part_chunk[i];
                    } else {
                        nelem *= (long)DRA[handle].chunk[i];
                    }
                }
            }
            nelem *= (long)(ds_a.lo[ndim-1]-1)%DRA[handle].chunk[ndim-1];
            offelem += (long)nelem;
        }
    } else {
        /* Find offset by calculating the number of chunks that must be
         * traversed to get to the corner of block containing the lower
         * coordinate index ds_a.lo[]. Then move into the block along
         * the last dimension to the point ds_a.lo[ndim-1]. */
        for (i=0; i<ndim; i++) {
            par_block[i] = 1;
            for (j=0; j<ndim; j++) {
                if (j < i) {
                    par_block[i] *= (long)cur_ld[j];
                } else if (j == i) {
                    if (i == ndim-1) {
                        /* special case for last dimension, which may represent
                         * a fraction of a chunk */
                        par_block[i] *= (long)(ds_a.lo[i]-1);
                    } else {
                        par_block[i] *= (long)(blocks[j]*DRA[handle].chunk[j]);
                    }
                } else {
                    par_block[i] *= (long)(DRA[handle].dims[j]);
                }
            }
        }
        offelem = 0;
        for (i=0; i<ndim; i++) offelem += (long)par_block[i];
    }

    *offset = (Off_t)offelem * dai_sizeofM(DRA[handle].type); 
}


/**
 * write zero at EOF for NDRA
 */
void ndai_zero_eof(Integer d_a)
{
    Integer handle = d_a+DRA_OFFSET, nelem, i;
    Integer zero[MAXDIM];
    char byte;
    Off_t offset;

    byte = (char)0;

    if(INDEPFILES(d_a) || DRA[handle].numfiles > 1) {

        Integer   CR=0, i=0, nblocks=0; 
        section_t ds_a;
        /* number of processors that do io */
        Integer   ioprocs=dai_io_procs(d_a); 
        /* node id of current process (if it does io) */
        Integer   iome = dai_io_nodeid(d_a);

        /* total number of blocks in the disk resident array */
        nblocks = 1;
        for (i=0; i<DRA[handle].ndim; i++) {
            nblocks *= (DRA[handle].dims[i]+DRA[handle].chunk[i]-1)
                / DRA[handle].chunk[i];
            zero[i] = 0;
        }
        nfill_sectionM(ds_a, d_a, DRA[handle].ndim, zero, zero); 

        /* search for the last block for each I/O processor */
        for(i = 0; i <ioprocs; i++){
            CR = nblocks - 1 -i;
            if(CR % ioprocs == iome) break;
        }
        if(CR<0) return; /* no blocks owned */

        nblock_to_sectM(&ds_a, CR); /* convert block number to section */
        ndai_file_location(ds_a, &offset);
        nelem = 1;
        for (i=0; i<DRA[handle].ndim; i++) nelem *= (ds_a.hi[i] - ds_a.lo[i] + 1);
        offset += ((Off_t)nelem) * dai_sizeofM(DRA[handle].type);

#         ifdef DEBUG
        printf("me=%d zeroing EOF (%d) at %ld bytes \n",iome,CR,offset);
#         endif
    } else {

        nelem = 1;
        for (i=0; i<DRA[handle].ndim; i++) nelem *= DRA[handle].dims[i];
        offset = ((Off_t)nelem) * dai_sizeofM(DRA[handle].type);
    }

    if(elio_write(DRA[handle].fd, offset-1, &byte, 1) != (Size_t)1)
        dai_error("ndai_zero_eof: write error ",0);

    /* This is a modification added by Sriram. Not sure what it is suppose
     * to do for you so I'm commenting it out for now. This function is
     * strictly an addition to existing code.
     elio_zero_eof(DRA[handle].fd);
     */
}


/**
 * SET CONFIGURATION FOR HANDLING DRAs STORED ON OPEN FILE SYSTEMS
 */
void dai_set_config(Integer numfiles, Integer numioprocs,
        Integer *number_of_files, Integer *io_procs)
{
    if (numfiles < 1) {
        if (numioprocs > 0) {
            numfiles = numioprocs;
        } else {
            numfiles = pnga_cluster_nnodes();
        }
    }
    if (numioprocs < 1) {
        numioprocs = numfiles;
    }
    *number_of_files = numfiles;
    *io_procs = numioprocs;
    if (*number_of_files > pnga_nnodes()) {
        if (pnga_nodeid() == 0) {
            printf("WARNING: Number of files requested exceeds number of\n");
            printf("processors. Value is reset to number of processors: %ld\n",
                    (long)pnga_nnodes());
        }
        *number_of_files = pnga_nnodes();
    }
    if (*io_procs > 1 && *number_of_files > 1) {
        if (*io_procs != *number_of_files) {
            if (pnga_nodeid() == 0) {
                printf("WARNING: Number of IO processors is not equal to the\n");
                printf("number of files requested. Number of IO processors\n");
                printf("is reset to number of files: %ld\n",(long)*number_of_files);
            }
            *io_procs = *number_of_files;
        }
    }
    if (*number_of_files == 1) {
        if (*io_procs > pnga_nnodes()) {
            if (pnga_nodeid() == 0) {
                printf("WARNING: Number of requested IO processors\n");
                printf("exceeds number of available processors. Number of IO\n");
                printf("processors reset to the number of available processors %ld\n",
                        (long)pnga_nnodes());
            }
            *io_procs = pnga_nnodes();
        }
    }
    if (*number_of_files > *io_procs) {
        if (pnga_nodeid() == 0) {
            printf("WARNING: Number of files is greater than\n");
            printf("number of IO processors. Number of files reset to number of\n");
            printf("IO processors: %ld",(long)*io_procs);
        }
        *number_of_files = *io_procs;
    }
}


/**
 * CREATE AN N-DIMENSIONAL DISK ARRAY WITH USER SPECIFIED IO CONFIGURATION
 *
 * @param type[in]
 * @param ndim[in]       dimension of DRA
 * @param dims[in]       dimensions of DRA
 * @param name[in]
 * @param filename[in]
 * @param mode[in]
 * @param reqdims[in]    dimension of typical request
 * @param numfiles[in]   number of files for DRA
 * @param numioprocs[in] number of IO procs to use
 * @param d_a[out]       DRA handle
 */
Integer ndrai_create_config(Integer *type, Integer *ndim, Integer dims[],
        char *name, char *filename, Integer *mode, Integer reqdims[],
        Integer *numfiles, Integer *numioprocs, Integer *d_a)
{
    Integer handle, elem_size, ctype, i;
    int emode;

    /* convert Fortran to C data type */
    ctype = pnga_type_f2c(*type);
    pnga_sync();

    /* if we have an error here, it is fatal */       
    dai_check_typeM(ctype);    
    for (i=0; i<*ndim; i++) if (dims[i] <=0)
        dai_error("ndra_create: disk array dimension invalid ", dims[i]);
    if(strlen(filename)>DRA_MAX_FNAME)
        dai_error("ndra_create: filename too long", DRA_MAX_FNAME);

    /*** Get next free DRA handle ***/
    if( (handle = dai_get_handle()) == -1)
        dai_error("ndra_create: too many disk arrays ", _max_disk_array);
    *d_a = handle - DRA_OFFSET;

    /*translate DRA mode into ELIO mode*/
    emode = dai_elio_mode((int)*mode);

    /* Determine array configuration */
    dai_set_config(*numfiles, *numioprocs, &DRA[handle].numfiles,
            &DRA[handle].ioprocs);

    /* determine disk array decomposition */ 
    elem_size = dai_sizeofM(ctype);
    ndai_chunking( elem_size, *ndim, reqdims, dims, DRA[handle].chunk);

    /* determine layout -- by row or column */
    DRA[handle].layout = COLUMN;

    /* complete initialization */
    for (i=0; i<*ndim; i++) DRA[handle].dims[i] = dims[i];
    DRA[handle].ndim = *ndim;
    DRA[handle].type = ctype;
    DRA[handle].mode = (int)*mode;
    strncpy (DRA[handle].fname, filename,  DRA_MAX_FNAME);
    strncpy(DRA[handle].name, name, DRA_MAX_NAME );

    dai_write_param(DRA[handle].fname, *d_a);      /* create param file */
    DRA[handle].indep = dai_file_config(filename); /*check file configuration*/

    /* create file */
    if(dai_io_manage(*d_a)){ 

        if (INDEPFILES(*d_a) || DRA[handle].numfiles > 1) {

            sprintf(dummy_fname,"%s.%ld",DRA[handle].fname,(long)dai_io_nodeid(*d_a));
            DRA[handle].fd = elio_open(dummy_fname,emode, ELIO_PRIVATE);
        } else{

            DRA[handle].fd = elio_open(DRA[handle].fname,emode, ELIO_SHARED); 
        }

        if(DRA[handle].fd==NULL)dai_error("ndra_create:failed to open file",0);
        if(DRA[handle].fd->fd==-1)dai_error("ndra_create:failed to open file",-1);
    }

    /*
     *  Need to zero the last element in the array on the disk so that
     *  we never read beyond EOF.
     *
     *  For multiple component files will stamp every one of them.
     *
     */
    pnga_sync();

    if(dai_file_master(*d_a) && dai_write_allowed(*d_a)) ndai_zero_eof(*d_a);

    pnga_sync();

    return(ELIO_OK);
}


/**
 * CREATE AN N-DIMENSIONAL DISK ARRAY
 *
 * @param type[in]
 * @param ndim[in]     dimension of DRA
 * @param dims[][in]   dimensions of DRA
 * @param name[in]
 * @param filename[in]
 * @param mode[in]
 * @param reqdims[in]  dimension of typical request
 * @param d_a[out]     DRA handle
 */
Integer ndrai_create(Integer *type, Integer *ndim, Integer dims[], char *name,
        char *filename, Integer *mode, Integer reqdims[], Integer *d_a)
{
    Integer ret;
    Integer files = _dra_number_of_files;
    Integer procs = _dra_io_procs;
    ret = ndrai_create_config(type, ndim, dims, name, filename, mode, reqdims,
            &files, &procs, d_a);
    return ret;
}


/**
 * CREATE A 2-D DISK ARRAY
 *
 * @param type[in]
 * @param dim1[in]
 * @param dim2[in]
 * @param name[in]
 * @param filename[in]
 * @param mode[in]
 * @param reqdim1[in] dim1 of typical request
 * @param reqdim2[in] dim2 of typical request
 * @param d_a[out]    DRA handle
 */
Integer drai_create(Integer *type, Integer *dim1, Integer *dim2, char *name,
        char *filename, Integer *mode, Integer *reqdim1, Integer *reqdim2,
        Integer *d_a)
{
    Integer ndim = 2;
    Integer dims[2], reqdims[2];
    dims[0] = *dim1; dims[1] = *dim2;
    reqdims[0] = *reqdim1; reqdims[1] = *reqdim2;
    return ndrai_create(type, &ndim, dims, name, filename, mode, reqdims, d_a);

#if 0
    Integer handle, elem_size, ctype;
    int emode;

    /* convert Fortran to C data type */
    ctype = pnga_type_f2c(*type);
    pnga_sync();

    /* if we have an error here, it is fatal        */
    dai_check_typeM(ctype);    
    if( *dim1 <= 0 )
        dai_error("dra_create: disk array dimension1 invalid ",  *dim1);
    else if( *dim2 <= 0)
        dai_error("dra_create: disk array dimension2 invalid ",  *dim2);
    if(strlen(filename)>DRA_MAX_FNAME)
        dai_error("dra_create: filename too long", DRA_MAX_FNAME);

    /*  Get next free DRA handle */
    if( (handle = dai_get_handle()) == -1)
        dai_error("dai_create: too many disk arrays ", _max_disk_array);
    *d_a = handle - DRA_OFFSET;

    /*translate DRA mode into ELIO mode*/
    emode = dai_elio_mode((int)*mode);

    /* determine disk array decomposition  */
    elem_size = dai_sizeofM(ctype);
    dai_chunking( elem_size, *reqdim1, *reqdim2, *dim1, *dim2, 
            &DRA[handle].chunk[0], &DRA[handle].chunk[1]);

    /* determine layout -- by row or column */
    DRA[handle].layout = COLUMN;

    /* complete initialization */
    DRA[handle].dims[0] = *dim1;
    DRA[handle].dims[1] = *dim2;
    DRA[handle].ndim = 2;
    DRA[handle].type = ctype;
    DRA[handle].mode = (int)*mode;
    strncpy (DRA[handle].fname, filename,  DRA_MAX_FNAME);
    strncpy(DRA[handle].name, name, DRA_MAX_NAME );

    DRA[handle].ioprocs = _dra_io_procs;
    DRA[handle].numfiles = _dra_number_of_files;

    dai_write_param(DRA[handle].fname, *d_a);      
    DRA[handle].indep = dai_file_config(filename); 
    /* create file */
    if(dai_io_manage(*d_a)){ 

        if (INDEPFILES(*d_a) || DRA[handle].numfiles > 1) {

            sprintf(dummy_fname,"%s.%ld",DRA[handle].fname,(long)dai_io_nodeid(*d_a));
            DRA[handle].fd = elio_open(dummy_fname,emode, ELIO_PRIVATE);
        } else {

            DRA[handle].fd = elio_open(DRA[handle].fname,emode, ELIO_SHARED); 
        }

        if(DRA[handle].fd==NULL)dai_error("dra_create:failed to open file",0);
        if(DRA[handle].fd->fd==-1)dai_error("dra_create:failed to open file",0);
    }


    pnga_sync();

    if(dai_file_master(*d_a) && dai_write_allowed(*d_a)) dai_zero_eof(*d_a);

    pnga_sync();

    return(ELIO_OK);
#endif

}


/**
 * write N-dimensional aligned block of data from memory buffer to d_a
 *
 * @param ds_a[in] section of DRA written to disk
 * @param buf[in]  pointer to io buffer
 * @param ld[in]   array of strides
 * @param id
 */
void ndai_put(section_t ds_a, void *buf, Integer ld[], io_request_t *id)
{
    Integer handle = ds_a.handle + DRA_OFFSET, elem, i;
    Integer ndim = ds_a.ndim;
    Off_t   offset;
    Size_t  bytes;
#if WALLTIME
    double ss0,tt0,tt1;
#endif

    /* find location in a file where data should be written */
    ndai_file_location(ds_a, &offset);
    for (i=0; i<ndim-1; i++) if ((ds_a.hi[i]-ds_a.lo[i]+1) != ld[i])
        dai_error("ndai_put: bad ld",ld[i]); 

    /* since everything is aligned, write data to disk */
    elem = 1;
    for (i=0; i<ndim; i++) elem *= (ds_a.hi[i]-ds_a.lo[i]+1);
    bytes= (Size_t) elem * dai_sizeofM(DRA[handle].type);
#if WALLTIME
    walltime_(&ss0,&tt0);
#endif
    if( ELIO_OK != elio_awrite(DRA[handle].fd, offset, buf, bytes, id ))
        dai_error("ndai_put failed", ds_a.handle);
#if WALLTIME
    walltime_(&ss0,&tt1);
    printf("p[%d] Beginning ndai_put: %16.6f\n",pnga_nodeid(),tt0);
    printf("p[%d] Ending ndai_put: %16.6f\n",pnga_nodeid(),tt1);
#endif
}


/**
 * read N-dimensional aligned block of data from d_a to memory buffer
 *
 * @param ds_a[in] section of DRA read from disk
 * @param buf[in]  pointer to io buffer
 * @param ld[in]   array of strides
 * @param id
 */
void ndai_get(section_t ds_a, void *buf, Integer ld[], io_request_t *id)
{
    Integer handle = ds_a.handle + DRA_OFFSET, elem, rc;
    Integer ndim = DRA[handle].ndim, i;
    Off_t   offset;
    Size_t  bytes;
#if WALLTIME
    double ss0,tt0,tt1;
#endif

    /* find location in a file where data should be read from */
    ndai_file_location(ds_a, &offset);

#ifdef CLEAR_BUF
    dai_clear_buffer();
#endif

    for (i=0; i<ndim-1; i++) if ((ds_a.hi[i] - ds_a.lo[i] + 1) != ld[i])
        dai_error("ndai_get: bad ld",ld[i]); 
    /* since everything is aligned, read data from disk */
    elem = 1;
    for (i=0; i<ndim; i++) elem *= (ds_a.hi[i]-ds_a.lo[i]+1);
    bytes= (Size_t) elem * dai_sizeofM(DRA[handle].type);
#if WALLTIME
    walltime_(&ss0,&tt0);
#endif
    rc= elio_aread(DRA[handle].fd, offset, buf, bytes, id );
#if WALLTIME
    walltime_(&ss0,&tt1);
    printf("p[%d] Beginning ndai_get: %16.6f\n",pnga_nodeid(),tt0);
    printf("p[%d] Ending ndai_get:    %16.6f\n",pnga_nodeid(),tt1);
#endif
    if(rc !=  ELIO_OK) dai_error("ndai_get failed", rc);
}


/**
 * decompose section defined by lo and hi into aligned and unaligned DRA
 * subsections
 *
 * @param ds_a[in]
 * @param aligned[out]   Indices of aligned subsections.
 * @param na[out]        Number of aligned subsections.
 * @param cover[out]     Indices of cover subsections.
 * @param unaligned[out] Indices of unaligned subsections.
 * @param nu[out]        Number of unaligned subsections.
 */ 
void ndai_decomp_section(section_t ds_a, Integer aligned[][2*MAXDIM],
        int *na, Integer cover[][2*MAXDIM], Integer unaligned[][2*MAXDIM],
        int *nu)
{
    Integer a=0, u=0, handle=ds_a.handle+DRA_OFFSET;
    Integer chunk_units;
    Integer i, j, idir, ndim = DRA[handle].ndim;
    Integer off_low[MAXDIM], off_hi[MAXDIM];
    Integer cover_lo[MAXDIM], cover_hi[MAXDIM];
    Integer check, chunk_lo, chunk_hi;

    /* section[lo,hi] is decomposed into 'aligned' and 'unaligned'
     * subsections.  The aligned subsections are aligned on
     * chunk[1]..chunk[ndim-1] boundaries. The unaligned subsections are
     * not completely covered by chunk[1]..chunk[ndim]-1 boundaries. These
     * are subsets of the 'cover' subsections which are aligned on chunk
     * boundaries and contain the unaligned subsections. Disk I/O will
     * actually be performed on 'aligned' and 'cover' subsections instead
     * of 'unaligned' subsections.
     *
     * The indexing of the aligned[][idir], cover[][idir], and
     * unaligned[][idir] arrays is idir = 0,1 corresponds to low and
     * high values in the 0 direction, idir = 2,3 corresponds to low and
     * high values in the 1 direction and so on up to a value of
     * idir = 2*ndim-1.
     *
     * The strategy for performing the decomposition is to first find the
     * coordinates corresponding to an aligned patch that completely covers
     * the originally requested section.
     * 
     * Begin by initializing some arrays. */

    for (i=0, j=0; i<ndim; i++) {
        aligned[a][j] = ds_a.lo[i];
        cover_lo[i] = ds_a.lo[i];
        off_low[i] = (ds_a.lo[i] - 1) % DRA[handle].chunk[i];
        j++;
        aligned[a][j] = ds_a.hi[i];
        cover_hi[i] = ds_a.hi[i];
        off_hi[i] = ds_a.hi[i] % DRA[handle].chunk[i];
        j++;
    }
    /* Find coordinates of aligned patch that completely covers the first
       ndim-1 dimensions of ds_a */
    for (i=0; i<ndim-1; i++) {
        if (off_low[i] !=0) {
            chunk_lo = (ds_a.lo[i] - 1) / DRA[handle].chunk[i];
            cover_lo[i] = chunk_lo * DRA[handle].chunk[i] + 1;
        }
        if (off_hi[i] !=0) {
            chunk_hi = ds_a.hi[i] / DRA[handle].chunk[i] + 1;
            cover_hi[i] = chunk_hi * DRA[handle].chunk[i];
            if (cover_hi[i] > DRA[handle].dims[i])
                cover_hi[i] = DRA[handle].dims[i];
        }
    }
    /* Find coordinates of aligned chunk (if there is one) */
    j = 0;
    check = 1;
    for (i=0; i<ndim-1; i++) {
        if (off_low[i] != 0) {
            chunk_lo = (ds_a.lo[i] - 1) / DRA[handle].chunk[i] + 1;
            aligned[a][j] = chunk_lo * DRA[handle].chunk[i] + 1;
        }
        j++;
        if (off_hi[i] !=0) {
            chunk_hi = ds_a.hi[i] / DRA[handle].chunk[i];
            aligned[a][j] = chunk_hi * DRA[handle].chunk[i];
        }
        if (aligned[a][j] < aligned[a][j-1]) check = 0;
        j++;
    }
    *na = (check == 1)  ?1 :0;

    /* evaluate cover sections and unaligned chunk dimensions. We
       break the evaluation of chunks into the following two
cases:
1) There is no aligned section
2) There is an aligned section */

    if (*na == 0) {
        /* There is no aligned block. Just return with one cover
           section */
        for (i=0, j=0; i<ndim; i++) {
            cover[u][j] = cover_lo[i];
            unaligned[u][j] = ds_a.lo[i];
            j++;
            cover[u][j] = cover_hi[i];
            unaligned[u][j] = ds_a.hi[i];
            j++;
        }
        *nu = 1;
        return;
    }

    /* An aligned chunk exists so we must find cover sections that
       surround it. We scan over the coordinate directions idir
       and choose cover sections such that if the coordinate
       direction of the cover section is less than idir, then the
       cover section extends to the boundary of the aligned
       section. If the coordinate direction of the cover section
       is greater than idir then the cover extends beyond the
       dimensions of the aligned chunk (if there is a nonzero
       offset). This scheme guarantees that corner pieces of the
       sections are picked up once and only once. */

    for (idir=0; idir<ndim-1; idir++) {
        check = 1;
        /* cover over lower end of patch */
        if (off_low[idir] != 0) {
            for (i=0, j=0; i<ndim-1; i++) {
                if (i < idir) {
                    if (off_low[i] != 0) {
                        chunk_units = (ds_a.lo[i] - 1) / DRA[handle].chunk[i];
                        cover[u][j] = chunk_units * DRA[handle].chunk[i] + 1;
                    } else {
                        cover[u][j] = ds_a.lo[i];
                    }
                    unaligned[u][j] = ds_a.lo[i];
                    j++;
                    if (off_hi[i] != 0) {
                        chunk_units = ds_a.hi[i] / DRA[handle].chunk[i]+1;
                        cover[u][j] = PARIO_MIN(chunk_units * DRA[handle].chunk[i],
                                DRA[handle].dims[i]);
                    } else {
                        cover[u][j] = ds_a.hi[i];
                    }
                    unaligned[u][j] = ds_a.hi[i];
                    j++;
                } else if (i == idir) {
                    chunk_units = (ds_a.lo[i] - 1) / DRA[handle].chunk[i];
                    cover[u][j] = chunk_units * DRA[handle].chunk[i] + 1;
                    unaligned[u][j] = ds_a.lo[i];
                    j++;
                    cover[u][j] = PARIO_MIN(cover[u][j-1] + DRA[handle].chunk[i]-1,
                            DRA[handle].dims[i]);
                    unaligned[u][j] = PARIO_MIN(ds_a.hi[i],cover[u][j]);
                    j++;
                } else {
                    if (off_low[i] != 0) {
                        chunk_units = (ds_a.lo[i] - 1) / DRA[handle].chunk[i]+1;
                        cover[u][j] = chunk_units * DRA[handle].chunk[i] + 1;
                    } else {
                        cover[u][j] = ds_a.lo[i];
                    }
                    unaligned[u][j] = ds_a.lo[i];
                    j++;
                    if (off_hi[i] != 0) {
                        chunk_units = ds_a.hi[i] / DRA[handle].chunk[i];
                        cover[u][j] = chunk_units * DRA[handle].chunk[i];
                    } else {
                        cover[u][j] = ds_a.hi[i];
                    }
                    unaligned[u][j] = ds_a.hi[i];
                    j++;
                }
            }
            cover[u][j] = ds_a.lo[ndim-1];
            unaligned[u][j] = ds_a.lo[ndim-1];
            j++;
            cover[u][j] = ds_a.hi[ndim-1];
            unaligned[u][j] = ds_a.hi[ndim-1];
            u++;
            check = 1;
            aligned[a][2*idir] = cover[u-1][2*idir+1]+1;
        }
        /* check to see if there is only one unaligned section covering this
           dimension */
        if (check == 1) {
            if (cover[u-1][2*idir+1] >= ds_a.hi[idir]) check = 0;
        } else {
            check = 1;
        } 
        /* handle cover over upper end of patch */
        if (off_hi[idir] != 0 && check == 1) {
            for (i=0, j=0; i<ndim-1; i++) {
                if (i < idir) {
                    if (off_low[i] != 0) {
                        chunk_units = (ds_a.lo[i] - 1) / DRA[handle].chunk[i];
                        cover[u][j] = chunk_units * DRA[handle].chunk[i] + 1;
                    } else {
                        cover[u][j] = ds_a.lo[i];
                    }
                    unaligned[u][j] = ds_a.lo[i];
                    j++;
                    if (off_hi[i] != 0) {
                        chunk_units = ds_a.hi[i] / DRA[handle].chunk[i]+1;
                        cover[u][j] = PARIO_MIN(chunk_units * DRA[handle].chunk[i],
                                DRA[handle].dims[i]);
                    } else {
                        cover[u][j] = ds_a.hi[i];
                    }
                    unaligned[u][j] = ds_a.hi[i];
                    j++;
                } else if (i == idir) {
                    chunk_units = ds_a.hi[i] / DRA[handle].chunk[i];
                    cover[u][j] = chunk_units * DRA[handle].chunk[i] + 1;
                    unaligned[u][j] = cover[u][j];
                    aligned[a][2*i+1] = PARIO_MIN(cover[u][j]-1,ds_a.hi[idir]);
                    j++;
                    cover[u][j] = PARIO_MIN(cover[u][j-1] + DRA[handle].chunk[i]-1,
                            DRA[handle].dims[i]);
                    unaligned[u][j] = PARIO_MIN(ds_a.hi[i],cover[u][j]);
                    j++;
                } else {
                    if (off_low[i] != 0) {
                        chunk_units = (ds_a.lo[i] - 1) / DRA[handle].chunk[i]+1;
                        cover[u][j] = chunk_units * DRA[handle].chunk[i] + 1;
                    } else {
                        cover[u][j] = ds_a.lo[i];
                    }
                    unaligned[u][j] = ds_a.lo[i];
                    j++;
                    if (off_hi[i] != 0) {
                        chunk_units = ds_a.hi[i] / DRA[handle].chunk[i];
                        cover[u][j] = chunk_units * DRA[handle].chunk[i];
                    } else {
                        cover[u][j] = ds_a.hi[i];
                    }
                    unaligned[u][j] = ds_a.hi[i];
                    j++;
                }
            }
            cover[u][j] = ds_a.lo[ndim-1];
            unaligned[u][j] = ds_a.lo[ndim-1];
            j++;
            cover[u][j] = ds_a.hi[ndim-1];
            unaligned[u][j] = ds_a.hi[ndim-1];
            u++;
            aligned[a][2*idir+1] = cover[u-1][2*idir]-1;
        }
    }
    *nu = (int)u;
    return;
}


/**
 * given the set of indices lo inside the patch cover, find the next
 * set of indices assuming that the area represented by cover has been
 * divided up into blocks whose size is given in inc
 */
int ndai_next(Integer *lo, Integer *cover, Integer *inc, Integer ndim)
{
    /* first check to see if any of the low indices are out of range.
       If so then reset all low indices to minimum values. */
    int retval=1;
    Integer i;
    for (i = 0; i<ndim; i++) { 
        if (lo[i] == 0) retval = 0;
    }
    if (retval == 0) {
        for (i = 0; i<ndim; i++) { 
            lo[i] = cover[2*i];
        }
    }
    /* increment all indices in lo. If index exceeds value cover[2*i+1]
       for that index, then set index back to cover[2*i] and increment
       next index. */
    if (retval != 0) {
        for (i=0; i<ndim; i++) {
            lo[i] += inc[i];
            if (lo[i] > cover[2*i+1]) {
                if (i<ndim-1) lo[i] = cover[2*i];
            } else {
                break;
            }
        }
    }
    retval = (lo[ndim-1] <= cover[2*ndim-1]);
    return retval;
}


/**
 * compute next chunk of array to process
 */
int ndai_next_chunk(Integer req, Integer* list, section_t* ds_chunk)
{
    Integer   handle = ds_chunk->handle+DRA_OFFSET;
    int       retval, ndim = DRA[handle].ndim, i;

    /* If we are writing out to multiple files then we need to consider
       chunk boundaries along last dimension */
    /*    if(INDEPFILES(ds_chunk->handle) || DRA[handle].numfiles > 1) */
    if(ds_chunk->lo[ndim-1] && DRA[handle].chunk[ndim-1]>1) 
        ds_chunk->lo[ndim-1] -= (ds_chunk->lo[ndim-1] -1) %
            DRA[handle].chunk[ndim-1];

    /* ds_chunk->lo is getting set in this call. list contains the
       the lower and upper indices of the cover section. */
    retval = ndai_next(ds_chunk->lo, list, DRA[handle].chunk, ndim);
    /*
       printf("Request %d\n",req);
       for (i=0; i<ndim; i++) {
       printf("ds_chunk.lo[%d] = %d cover.lo[%d] = %d cover.hi[%d] = %d\n", i,
       ds_chunk->lo[i], i, list[2*i], i, list[2*i+1]);
       } */
    if(!retval) {
        return(retval);
    }

    for (i=0; i<ndim; i++) {
        ds_chunk->hi[i] = PARIO_MIN(list[2*i+1],
                ds_chunk->lo[i]+DRA[handle].chunk[i]-1);
    }

    /* Again, if we are writing out to multiple files then we need to consider
       chunk boundaries along last dimension */
    /*    if(INDEPFILES(ds_chunk->handle) || DRA[handle].numfiles > 1) {  */
    if (1) {
        Integer nlo;
        Integer hi_temp =  ds_chunk->lo[ndim-1] +
            DRA[handle].chunk[ndim-1] -1;
        hi_temp -= hi_temp % DRA[handle].chunk[ndim-1];
        ds_chunk->hi[ndim-1] = PARIO_MIN(ds_chunk->hi[ndim-1], hi_temp); 

        /*this line was absent from older version on bonnie that worked */
        nlo = 2*(ndim-1);
        if(ds_chunk->lo[ndim-1] < list[nlo]) ds_chunk->lo[ndim-1] = list[nlo]; 
    } 
    /*
       for (i=0; i<ndim; i++) {
       printf("ds_chunk.hi[%d] = %d\n", i, ds_chunk->hi[i]);
       } */

    return 1;
}


/**
 * function to complete an operation and release the buffer associated
 * with a buffer id
 */
void wait_buf(char *buf)
{
    Integer *ga_movhdl;
    io_request_t *io_req;
    int op_code; /* int buf_id = nbuf; */
    buf_info *bi;

    if (buf == NULL) return;

    bi = (buf_info*) buf;
    op_code = bi->op;
    io_req = &(bi->io_req);
    ga_movhdl = &(bi->ga_movhdl);

    /*if (buf_id >= nbuf) {
      printf("Wait_buf Error: No operation is associated with this buffer\n");
      return;
      }*/

    switch(op_code) {
        case DRA_OP_WRITE:
            elio_wait(io_req);
            break;

        case DRA_OP_READ:
            if (bi->align == 0)
                pnga_nbwait(ga_movhdl);
            else {
                elio_wait(io_req);
                dai_exec_callback(buf, WAIT);
            }
            break;

        default:
            return;
    }

#ifdef BUF_DEBUG
    printf("Released a buffer\n");
#endif

}


/**
 * Write or Read Unaligned Subsections to/from disk: 
 * always read an aligned extension of a section from disk to local buffer then 
 * for read : copy requested data from buffer to global array;
 * for write: overwrite part of buffer with data from g_a and write
 * complete buffer to disk
 *
 * @param opcode[in] signal for read or write
 * @param transp[in] should data be transposed
 * @param ds_a[in] section of DRA that is to be read from or written to
 * @param gs_a[in] section of GA that is to be read from or written to
 * @param req[in] request number
 */
void ndai_transfer_unlgn(
        int opcode,
        int transp,
        section_t ds_a,
        section_t gs_a,
        Integer req)
{
    Integer   chunk_ld[MAXDIM],  next, offset, i, j;
    int   type = DRA[ds_a.handle+DRA_OFFSET].type;
    Integer   ndim = DRA[ds_a.handle+DRA_OFFSET].ndim;
    section_t ds_chunk, ds_unlg;
    char      *buf, *buffer; 
    Integer *ga_movhdl;
    io_request_t *io_req;
    buf_info *bi;

    ds_chunk =  ds_unlg = ds_a;
    if (dra_debug_flag && 0) {
        for (i=0; i<ndim; i++) {
            printf("ndai_transfer_unlgn: ds_chunk.lo[%ld] = %ld\n",
                    (long)i,(long)ds_chunk.lo[i]);
            printf("ndai_transfer_unlgn: ds_chunk.hi[%ld] = %ld\n",
                    (long)i,(long)ds_chunk.hi[i]);
        }
        printf("ndai_transfer_unlgn: number of unaligned chunks = %d\n",
                Requests[req].nu);
        for (j=0; j<Requests[req].nu; j++) {
            for (i=0; i<ndim; i++) {
                printf("ndai_transfer_unlgn: list_cover[%ld][%ld] = %ld\n",
                        (long)j,(long)2*i,(long)
                        Requests[req].list_cover[j][2*i]);
                printf("ndai_transfer_unlgn: list_cover[%ld][%ld] = %ld\n",
                        (long)j,(long)2*i+1,
                        (long)Requests[req].list_cover[j][2*i+1]);
            }
        }
    }

    for(next = 0; next < Requests[req].nu; next++){

        for (i=0; i<ndim; i++) ds_chunk.lo[i] = 0;   /* initialize */
        while(ndai_next_chunk(req, Requests[req].list_cover[next],&ds_chunk)){

            if(dai_myturn(ds_chunk)){

                /*find corresponding to chunk of 'cover' unaligned sub-subsection*/
                for (i=0; i<ndim; i++) {
                    ds_unlg.lo[i] = Requests[req].list_unlgn[next][2*i];
                    ds_unlg.hi[i] = Requests[req].list_unlgn[next][2*i+1];
                }

                if (dra_debug_flag && 0) {
                    for (i=0; i<ndim; i++) {
                        printf("ndai_transfer_unlgn: ds_chunk.lo[%ld] = %ld\n",
                                (long)i,(long)ds_chunk.lo[i]);
                        printf("ndai_transfer_unlgn: ds_chunk.hi[%ld] = %ld\n",
                                (long)i,(long)ds_chunk.hi[i]);
                    }
                    for (i=0; i<ndim; i++) {
                        printf("ndai_transfer_unlgn: ds_unlg.lo[%ld] = %ld\n",
                                (long)i,(long)ds_unlg.lo[i]);
                        printf("ndai_transfer_unlgn: ds_unlg.hi[%ld] = %ld\n",
                                (long)i,(long)ds_unlg.hi[i]);
                    }
                }

                if(!dai_section_intersect(ds_chunk, &ds_unlg))
                    dai_error("ndai_transfer_unlgn: inconsistent cover", 0);

                /* copy data from disk to DRA buffer */
                for (i=0; i<ndim-1; i++) chunk_ld[i] = ds_chunk.hi[i] - ds_chunk.lo[i] + 1;
                /* get a free buffer */
                buf = get_buf(&buf_ctxt, Requests[req].call_id);

                bi = (buf_info*) buf;
                io_req = &(bi->io_req);
                ga_movhdl = &(bi->ga_movhdl);
                bi->align = 0;

                buf = (char*) (buf + sizeof(buf_info));

                ndai_get(ds_chunk, buf, chunk_ld, io_req);
                elio_wait(io_req); 
                /* determine location in the buffer where GA data should be */
                offset = ds_unlg.lo[ndim-1]-ds_chunk.lo[ndim-1];
                for (i=ndim-2; i>=0; i--)  {
                    offset = offset*chunk_ld[i];
                    offset += ds_unlg.lo[i] - ds_chunk.lo[i];
                }
                buffer  = (char*)buf;
                buffer += offset * dai_sizeofM(type);

                switch (opcode){
                    case DRA_OP_WRITE: 
                        bi->op = DRA_OP_WRITE;
                        /* overwrite a part of buffer with data from g_a */  
                        nga_move(LOAD, transp, gs_a, ds_a, ds_unlg, buffer, chunk_ld, ga_movhdl);
                        pnga_nbwait(ga_movhdl);

                        /* write ENTIRE updated buffer back to disk */
                        ndai_put(ds_chunk, buf, chunk_ld, io_req);

                        break;

                    case DRA_OP_READ:
                        bi->op = DRA_OP_READ;
                        /* copy requested data from buffer to g_a */
                        nga_move(STORE, transp, gs_a, ds_a, ds_unlg, buffer, chunk_ld, ga_movhdl);

                        break;

                    default:
                        dai_error("dai_transfer_unlg: invalid opcode",(Integer)opcode);
                }

#       ifdef DEBUG
                fprintf(stderr,"%d transf unlg g[%d:%d,%d:%d]-d[%d:%d,%d:%d]\n",
                        dai_io_nodeid(), gs_chunk.lo[0], gs_chunk.hi[0],
                        gs_chunk.lo[1], gs_chunk.hi[1],
                        ds_unlg.lo[0], ds_unlg.hi[0],
                        ds_unlg.lo[1], ds_unlg.hi[1]);
#       endif

            }
        }
    }
    /*
       returning from this function leaving some outstanding operations,
       so that dra_read()/write() can be non-blocking to some extent. we will
       have to call dra_wait() to make sure these operations are complete.
       */
}


/**
 * write or read aligned subsections to disk 
 */
void ndai_transfer_algn(
        int opcode, int transp, section_t ds_a, section_t gs_a, Integer req)
{
    Integer  next, chunk_ld[MAXDIM], ndim = ds_a.ndim;
    Integer i;
    section_t ds_chunk = ds_a;
    char *buf, *buffer;
    Integer *ga_movhdl;
    io_request_t *io_req;
    buf_info *bi;

    for(next = 0; next < Requests[req].na; next++){
        for (i=0; i<ndim; i++) ds_chunk.lo[i] = 0; /*initialize */

        while(ndai_next_chunk(req, Requests[req].list_algn[next], &ds_chunk)){
            if (dra_debug_flag && 0) { 
                printf("ndai_transfer_algn: Request %ld\n",(long)req);
                for (i=0; i<ndim; i++) {
                    printf("ndai_transfer_algn: ds_chunk.lo[%ld] = %ld\n",
                            (long)i,(long)ds_chunk.lo[i]);
                    printf("ndai_transfer_algn: ds_chunk.hi[%ld] = %ld\n",
                            (long)i,(long)ds_chunk.hi[i]);
                }
            }

            if(dai_myturn(ds_chunk)){

                for (i=0; i<ndim-1; i++) chunk_ld[i] = ds_chunk.hi[i] - ds_chunk.lo[i] + 1;
                /* get a free buffer */
                buf = get_buf(&buf_ctxt, Requests[req].call_id);

                bi = (buf_info*) buf;
                io_req = &(bi->io_req);
                ga_movhdl = &(bi->ga_movhdl);
                bi->align = 1;
                bi->callback = OFF;

                buffer = buf;
                buf = buf + sizeof(buf_info);

                switch (opcode){

                    case DRA_OP_WRITE:
                        bi->op = DRA_OP_WRITE;
                        /* copy data from g_a to DRA buffer */
                        nga_move(LOAD, transp, gs_a, ds_a, ds_chunk, buf, chunk_ld, ga_movhdl);
                        pnga_nbwait(ga_movhdl);

                        /* copy data from DRA buffer to disk */
                        ndai_put(ds_chunk, buf, chunk_ld, io_req);

                        break;

                    case DRA_OP_READ:
                        bi->op = DRA_OP_READ;
                        /* copy data from disk to DRA buffer */
                        ndai_get(ds_chunk, buf, chunk_ld, io_req);

                        /* copy data from DRA buffer to g_a */
                        /* nga_move(STORE, transp, gs_a, ds_a, ds_chunk, buf, chunk_ld, ga_movhdl); */
                        dai_callback(STORE, transp, gs_a, ds_a, ds_chunk, chunk_ld, buffer, req);
                        break;

                    default:
                        dai_error("dai_transfer_algn: invalid opcode",(Integer)opcode);
                }

#       ifdef DEBUG
                fprintf(stderr,"%d transf algn g[%d:%d,%d:%d]-d[%d:%d,%d:%d]\n",
                        dai_io_nodeid(), gs_chunk.lo[0], gs_chunk.hi[0],
                        gs_chunk.lo[1], gs_chunk.hi[1],
                        ds_chunk.lo[0], ds_chunk.hi[0],
                        ds_chunk.lo[1], ds_chunk.hi[1]);
#       endif

            }
        }
    }
    /*
       returning from this function leaving some outstanding operations,
       so that dra_read()/write() can be non-blocking to some extent. we will 
       have to call dra_wait() to make sure these operations are complete.
       */
}


/**
 * WRITE SECTION g_a[glo:ghi] TO d_a[dlo:dhi]
 *
 * @param transp[in] transpose operator
 * @param g_a[in] GA handle
 * @param glo[in]
 * @param ghi[in]
 * @param d_a[in] DRA handle
 * @param dlo[in]
 * @param dhi[in]
 * @param request[out] async. request id
 */
Integer FATR ndra_write_section_(logical *transp,
        Integer *g_a, Integer glo[], Integer ghi[],
        Integer *d_a, Integer dlo[], Integer dhi[],
        Integer *request)
{
    Integer gdims[MAXDIM], gtype, handle=*d_a+DRA_OFFSET;
    Integer i, gelem, delem, ndim;
    section_t d_sect, g_sect;

    pnga_sync();

    /* usual argument/type/range checking stuff */

    dai_check_handleM(*d_a,"ndra_write_sect");
    pnga_inquire(*g_a, &gtype, &ndim, gdims);
    if(!dai_write_allowed(*d_a))dai_error("ndra_write_sect: write not allowed",*d_a);
    if(DRA[handle].type != (int)gtype)dai_error("ndra_write_sect: type mismatch",gtype);
    if(DRA[handle].ndim != ndim)dai_error("ndra_write_sect: dimension mismatch", ndim);
    for (i=0; i<ndim; i++) dai_check_rangeM(glo[i], ghi[i], gdims[i],
            "ndra_write_sect: g_a dim error");
    for (i=0; i<ndim; i++) dai_check_rangeM(dlo[i], dhi[i], DRA[handle].dims[i],
            "ndra_write_sect: d_a dim error");

    /* check if numbers of elements in g_a & d_a sections match */
    gelem = 1;
    delem = 1;
    for (i=0; i<ndim; i++) {
        gelem *= (ghi[i]-glo[i]+1);
        delem *= (dhi[i]-dlo[i]+1);
    }
    if (gelem != delem)
        dai_error("ndra_write_sect: d_a and g_a sections do not match ", 0L);

    dai_assign_request_handle(request);

    /* decompose d_a section into aligned and unaligned subsections
     * -- with respect to underlying array layout on the disk
     */

    Requests[*request].nu=MAX_ALGN;    
    Requests[*request].na=MAX_UNLG;

    nfill_sectionM(d_sect, *d_a, DRA[handle].ndim, dlo, dhi); 
    nfill_sectionM(g_sect, *g_a, ndim, glo, ghi); 

    ndai_decomp_section(d_sect,
            Requests[*request].list_algn, 
            &Requests[*request].na,
            Requests[*request].list_cover, 
            Requests[*request].list_unlgn, 
            &Requests[*request].nu);
    _dra_turn = 0;

    /* process unaligned subsections */
    ndai_transfer_unlgn(DRA_OP_WRITE, (int)*transp, d_sect, g_sect, *request);

    /* process aligned subsections */
    ndai_transfer_algn (DRA_OP_WRITE, (int)*transp, d_sect, g_sect, *request);

    pnga_sync();

    return(ELIO_OK);
}


/**
 * WRITE N-dimensional g_a TO d_a
 *
 * @param g_a[in]      GA handle
 * @param d_a[in]      DRA handle
 * @param request[out] handle to async oper
 */
Integer FATR ndra_write_(Integer *g_a, Integer *d_a, Integer *request)
{
    Integer gdims[MAXDIM], gtype, handle=*d_a+DRA_OFFSET;
    logical transp = FALSE;
    Integer lo[MAXDIM], hi[MAXDIM], ndim, i;

    pnga_sync();

    /* usual argument/type/range checking stuff */

    dai_check_handleM(*d_a,"ndra_write");
    if( !dai_write_allowed(*d_a))
        dai_error("ndra_write: write not allowed to this array",*d_a);

    pnga_inquire(*g_a, &gtype, &ndim, gdims);
    if(DRA[handle].type != (int)gtype)dai_error("ndra_write: type mismatch",gtype);
    if(DRA[handle].ndim != ndim)dai_error("ndra_write: dimension mismatch",ndim);
    for (i=0; i<ndim; i++) {
        if(DRA[handle].dims[i] != gdims[i])
            dai_error("ndra_write: dims mismatch",gdims[i]);
    }

    /* right now, naive implementation just calls ndra_write_section */
    for (i=0; i<ndim; i++) {
        lo[i] = 1;
        hi[i] = DRA[handle].dims[i];
    }

    return(ndra_write_section_(&transp, g_a, lo, hi, d_a, lo, hi, request));
}


/**
 * READ SECTION g_a[glo:ghi] FROM d_a[dlo:dhi]
 *
 * @param transp[in]   transpose operator
 * @param g_a[in]      GA handle
 * @param glo[in]
 * @param ghi[in]
 * @param d_a[in]      DRA handle
 * @param dlo[in]
 * @param dhi[in]
 * @param request[out] request id
 */
Integer FATR ndra_read_section_(logical *transp,
        Integer *g_a, Integer glo[], Integer ghi[],
        Integer *d_a, Integer dlo[], Integer dhi[],
        Integer *request)
{
    Integer gdims[MAXDIM], gtype, handle=*d_a+DRA_OFFSET;
    Integer i, gelem, delem, ndim, me;
    section_t d_sect, g_sect;

    pnga_sync();
    me = pnga_nodeid();
    /* printf("%d: CAME HERE!!!", me); */
    /* usual argument/type/range checking stuff */
    dai_check_handleM(*d_a,"ndra_read_sect");
    if(!dai_read_allowed(*d_a))dai_error("ndra_read_sect: read not allowed",*d_a);
    pnga_inquire(*g_a, &gtype, &ndim, gdims);
    if(DRA[handle].type != (int)gtype)dai_error("ndra_read_sect: type mismatch",gtype);
    if(DRA[handle].ndim != ndim)dai_error("ndra_read_sect: dimension mismatch", ndim);

    for (i=0; i<ndim; i++) dai_check_rangeM(glo[i], ghi[i], gdims[i],
            "ndra_read_sect: g_a dim error");
    for (i=0; i<ndim; i++) dai_check_rangeM(dlo[i], dhi[i], DRA[handle].dims[i],
            "ndra_read_sect: d_a dim error");

    /* check if numbers of elements in g_a & d_a sections match */
    gelem = 1;
    delem = 1;
    for (i=0; i<ndim; i++) {
        gelem *= (ghi[i] - glo[i] + 1);
        delem *= (dhi[i] - dlo[i] + 1);
    }
    if (gelem != delem)
        dai_error("ndra_read_sect: d_a and g_a sections do not match ", 0L);

    dai_assign_request_handle(request);

    /* decompose d_a section into aligned and unaligned subsections
     * -- with respect to underlying array layout on the disk
     */

    Requests[*request].nu=MAX_ALGN;    
    Requests[*request].na=MAX_UNLG;

    if (dra_debug_flag) {
        for (i=0; i<ndim; i++) {
            printf("proc[%ld] ndra_read_section: dlo[%ld] = %ld\n",
                    (long)me, (long)i, (long)dlo[i]);
            printf("proc[%ld] ndra_read_section: dhi[%ld] = %ld\n",
                    (long)me, (long)i, (long)dhi[i]);
        }
        for (i=0; i<ndim; i++) {
            printf("proc[%ld] ndra_read_section: glo[%ld] = %ld\n",
                    (long)me, (long)i, (long)glo[i]);
            printf("proc[%ld] ndra_read_section: ghi[%ld] = %ld\n",
                    (long)me, (long)i, (long)ghi[i]);
        }
    }

    nfill_sectionM(d_sect, *d_a, DRA[handle].ndim, dlo, dhi); 
    nfill_sectionM(g_sect, *g_a, ndim, glo, ghi); 

    ndai_decomp_section(d_sect,
            Requests[*request].list_algn, 
            &Requests[*request].na,
            Requests[*request].list_cover, 
            Requests[*request].list_unlgn, 
            &Requests[*request].nu);

    _dra_turn = 0;
    if (dra_debug_flag && 0) {
        printf("ndra_read_section: Number of aligned sections %d\n",
                Requests[*request].na);
        printf("ndra_read_section: Number of unaligned sections %d\n",
                Requests[*request].nu);
        for (i=0; i<2*ndim; i++) {
            printf("ndra_read_section: list_algn[%ld] =  %ld\n",
                    (long)i,(long)Requests[*request].list_algn[0][i]);
        }
        for (i=0; i<2*ndim; i++) {
            printf("ndra_read_section: list_cover[%ld] =  %ld\n",
                    (long)i,(long)Requests[*request].list_cover[0][i]);
        }
        for (i=0; i<2*ndim; i++) {
            printf("ndra_read_section: list_unlgn[%ld] =  %ld\n",
                    (long)i,(long)Requests[*request].list_unlgn[0][i]);
        } 
    }

    /* process unaligned subsections */
    ndai_transfer_unlgn(DRA_OP_READ, (int)*transp,  d_sect, g_sect, *request);

    /* process aligned subsections */
    ndai_transfer_algn (DRA_OP_READ, (int)*transp,  d_sect, g_sect, *request);

    /* printf(" %d: CAME at the end of ndra_read_section!", me); */
    return(ELIO_OK);
}


/**
 * READ N-dimensional g_a FROM d_a
 */
Integer FATR ndra_read_(Integer* g_a, Integer* d_a, Integer* request)
{
    Integer gdims[MAXDIM], gtype, handle=*d_a+DRA_OFFSET;
    logical transp = FALSE;
    Integer lo[MAXDIM], hi[MAXDIM], ndim, i;

    pnga_sync();
    /* printf("%d: CAME AT ndra_read_!!\n", pnga_nodeid()); */
    /* usual argument/type/range checking stuff */
    dai_check_handleM(*d_a,"ndra_read");
    if(!dai_read_allowed(*d_a))dai_error("ndra_read: read not allowed",*d_a);
    pnga_inquire(*g_a, &gtype, &ndim, gdims);
    /* printf("%d: CAME After pnga_inquire!!\n", pnga_nodeid()); */
    if(DRA[handle].type != (int)gtype)dai_error("ndra_read: type mismatch",gtype);
    if(DRA[handle].ndim != ndim)dai_error("ndra_read: dimension mismatch",ndim);
    for (i=0; i<ndim; i++) {
        if(DRA[handle].dims[i] != gdims[i])
            dai_error("ndra_read: dims mismatch",gdims[i]);
    }

    /* right now, naive implementation just calls ndra_read_section */
    for (i=0; i<ndim; i++) {
        lo[i] = 1;
        hi[i] = DRA[handle].dims[i];
    }
    return(ndra_read_section_(&transp, g_a, lo, hi, d_a, lo, hi, request));
}


/**
 * WRITE SECTION g_a[gilo:gihi, gjlo:gjhi] TO d_a[dilo:dihi, djlo:djhi]
 *
 * @param transp[in]   transpose operator
 * @param g_a[in]      GA handle
 * @param gilo[in]
 * @param gihi[in]
 * @param gjlo[in]
 * @param gjhi[in]
 * @param d_a[in]      DRA handle
 * @param dilo[in]
 * @param dihi[in]
 * @param djlo[in]
 * @param djhi[in]
 * @param request[out] async. request id
 */
Integer FATR dra_write_section_(logical *transp,
        Integer *g_a,Integer *gilo,Integer *gihi,Integer *gjlo,Integer *gjhi,
        Integer *d_a,Integer *dilo,Integer *dihi,Integer *djlo,Integer *djhi,
        Integer *request)
{
    Integer glo[2], ghi[2], dlo[2], dhi[2];

    glo[0] = *gilo;
    glo[1] = *gjlo;
    ghi[0] = *gihi;
    ghi[1] = *gjhi;

    dlo[0] = *dilo;
    dlo[1] = *djlo;
    dhi[0] = *dihi;
    dhi[1] = *djhi;

    return (ndra_write_section_(transp, g_a, glo, ghi, d_a, dlo, dhi, request));

    /*
       Integer gdim1, gdim2, gtype, handle=*d_a+DRA_OFFSET;
       section_t d_sect, g_sect;

       pnga_sync();

       dai_check_handleM(*d_a,"dra_write_sect");
       ga_inquire_internal_(g_a, &gtype, &gdim1, &gdim2);
       if(!dai_write_allowed(*d_a))dai_error("dra_write_sect: write not allowed",*d_a);
       if(DRA[handle].type != (int)gtype)dai_error("dra_write_sect: type mismatch",gtype);
       dai_check_rangeM(*gilo,*gihi, gdim1, "dra_write_sect: g_a dim1 error");
       dai_check_rangeM(*gjlo,*gjhi, gdim2, "dra_write_sect: g_a dim2 error");
       dai_check_rangeM(*dilo,*dihi,DRA[handle].dims[0],"dra_write_sect:d_a dim1 error");
       dai_check_rangeM(*djlo,*djhi,DRA[handle].dims[1],"dra_write_sect:d_a dim2 error");


       if ((*dihi - *dilo + 1) * (*djhi - *djlo + 1) !=
       (*gihi - *gilo + 1) * (*gjhi - *gjlo + 1))
       dai_error("dra_write_sect: d_a and g_a sections do not match ", 0L);

       dai_assign_request_handle(request);


       Requests[*request].nu=MAX_ALGN;    
       Requests[*request].na=MAX_UNLG;

       fill_sectionM(d_sect, *d_a, *dilo, *dihi, *djlo, *djhi); 
       fill_sectionM(g_sect, *g_a, *gilo, *gihi, *gjlo, *gjhi); 

       dai_decomp_section(d_sect,
       Requests[*request].list_algn, 
       &Requests[*request].na,
       Requests[*request].list_cover, 
       Requests[*request].list_unlgn, 
       &Requests[*request].nu);
       _dra_turn = 0;


       dai_transfer_unlgn(DRA_OP_WRITE, (int)*transp, d_sect, g_sect, *request);


       dai_transfer_algn (DRA_OP_WRITE, (int)*transp, d_sect, g_sect, *request);

       pnga_sync();

       return(ELIO_OK);
       */
}


/**
 * WRITE g_a TO d_a
 *
 * @param g_a[in]      GA handle
 * @param d_a[in]      DRA handle
 * @param request[out] handle to async oper.
 */
Integer FATR dra_write_(Integer *g_a, Integer *d_a, Integer *request)
{
    return (ndra_write_(g_a, d_a, request));
    /*
       Integer gdim1, gdim2, gtype, handle=*d_a+DRA_OFFSET;
       logical transp = FALSE;
       Integer ilo, ihi, jlo, jhi;

       pnga_sync();     

       dai_check_handleM(*d_a,"dra_write");
       if( !dai_write_allowed(*d_a))
       dai_error("dra_write: write not allowed to this array",*d_a);

       ga_inquire_internal_(g_a, &gtype, &gdim1, &gdim2);
       if(DRA[handle].type != (int)gtype)dai_error("dra_write: type mismatch",gtype);
       if(DRA[handle].dims[0] != gdim1)dai_error("dra_write: dim1 mismatch",gdim1);
       if(DRA[handle].dims[1] != gdim2)dai_error("dra_write: dim2 mismatch",gdim2);


       ilo = 1; ihi = DRA[handle].dims[0];
       jlo = 1; jhi = DRA[handle].dims[1];
       return(dra_write_section_(&transp, g_a, &ilo, &ihi, &jlo, &jhi,
       d_a, &ilo, &ihi, &jlo, &jhi,request));
       */
}


/**
 * READ SECTION g_a[gilo:gihi, gjlo:gjhi] FROM d_a[dilo:dihi, djlo:djhi]
 *
 * @param transp[in]   transpose operator
 * @param g_a[in]      GA handle
 * @param gilo[in]
 * @param gihi[in]
 * @param gjlo[in]
 * @param gjhi[in]
 * @param d_a[in]      DRA handle
 * @param dilo[in]
 * @param dihi[in]
 * @param djlo[in]
 * @param djhi[in]
 * @param request[out] request id
 */
Integer FATR dra_read_section_(logical *transp,
        Integer *g_a,Integer *gilo,Integer *gihi,Integer *gjlo,Integer *gjhi,
        Integer *d_a,Integer *dilo,Integer *dihi,Integer *djlo,Integer *djhi,
        Integer *request)
{
    Integer glo[2], ghi[2], dlo[2], dhi[2];

    glo[0] = *gilo;
    glo[1] = *gjlo;
    ghi[0] = *gihi;
    ghi[1] = *gjhi;

    dlo[0] = *dilo;
    dlo[1] = *djlo;
    dhi[0] = *dihi;
    dhi[1] = *djhi;

    return (ndra_read_section_(transp, g_a, glo, ghi, d_a, dlo, dhi, request));

    /*
       Integer gdim1, gdim2, gtype, handle=*d_a+DRA_OFFSET;
       section_t d_sect, g_sect;

       pnga_sync();

       dai_check_handleM(*d_a,"dra_read_sect");
       if(!dai_read_allowed(*d_a))dai_error("dra_read_sect: read not allowed",*d_a);
       ga_inquire_internal_(g_a, &gtype, &gdim1, &gdim2);
       if(DRA[handle].type != (int)gtype)dai_error("dra_read_sect: type mismatch",gtype);
       dai_check_rangeM(*gilo, *gihi, gdim1, "dra_read_sect: g_a dim1 error");
       dai_check_rangeM(*gjlo, *gjhi, gdim2, "dra_read_sect: g_a dim2 error");
       dai_check_rangeM(*dilo, *dihi,DRA[handle].dims[0],"dra_read_sect:d_a dim1 error");
       dai_check_rangeM(*djlo, *djhi,DRA[handle].dims[1],"dra_read_sect:d_a dim2 error");


       if ((*dihi - *dilo + 1) * (*djhi - *djlo + 1) !=
       (*gihi - *gilo + 1) * (*gjhi - *gjlo + 1))
       dai_error("dra_read_sect: d_a and g_a sections do not match ", 0L);

       dai_assign_request_handle(request);


       Requests[*request].nu=MAX_ALGN;    
       Requests[*request].na=MAX_UNLG;

       fill_sectionM(d_sect, *d_a, *dilo, *dihi, *djlo, *djhi); 
       fill_sectionM(g_sect, *g_a, *gilo, *gihi, *gjlo, *gjhi); 

       dai_decomp_section(d_sect,
       Requests[*request].list_algn, 
       &Requests[*request].na,
       Requests[*request].list_cover, 
       Requests[*request].list_unlgn, 
       &Requests[*request].nu);

       _dra_turn = 0;


       dai_transfer_unlgn(DRA_OP_READ, (int)*transp,  d_sect, g_sect, *request);


       dai_transfer_algn (DRA_OP_READ, (int)*transp,  d_sect, g_sect, *request);

       return(ELIO_OK);
       */
}


/**
 * READ g_a FROM d_a
 */
Integer FATR dra_read_(Integer* g_a, Integer* d_a, Integer* request)
{

    return (ndra_read_(g_a, d_a, request));

    /*
       Integer gdim1, gdim2, gtype, handle=*d_a+DRA_OFFSET;
       logical transp = FALSE;
       Integer ilo, ihi, jlo, jhi;

       pnga_sync();


       dai_check_handleM(*d_a,"dra_read");
       if(!dai_read_allowed(*d_a))dai_error("dra_read: read not allowed",*d_a);
       ga_inquire_internal_(g_a, &gtype, &gdim1, &gdim2);
       if(DRA[handle].type != (int)gtype)dai_error("dra_read: type mismatch",gtype);
       if(DRA[handle].dims[0] != gdim1)dai_error("dra_read: dim1 mismatch",gdim1);
       if(DRA[handle].dims[1] != gdim2)dai_error("dra_read: dim2 mismatch",gdim2);


       ilo = 1; ihi = DRA[handle].dims[0];
       jlo = 1; jhi = DRA[handle].dims[1];
       return(dra_read_section_(&transp, g_a, &ilo, &ihi, &jlo, &jhi,
       d_a, &ilo, &ihi, &jlo, &jhi, request));
       */
}


/**
 * INQUIRE PARAMETERS OF EXISTING N-DIMENSIONAL DISK ARRAY
 *
 * @param d_a[in]      DRA handle
 * @param type[out]
 * @param ndim[out]
 * @param dims[out]
 * @param name[out]
 * @param filename[out]
 */
Integer ndrai_inquire(Integer *d_a, Integer *type, Integer *ndim,
        Integer dims[], char *name, char *filename)
{
    Integer handle=*d_a+DRA_OFFSET;
    Integer i;

    dai_check_handleM(*d_a,"dra_inquire");

    *type = (Integer)DRA[handle].type;
    *ndim = DRA[handle].ndim;
    for (i=0; i<*ndim; i++) {
      dims[i] = DRA[handle].dims[i];
    }
    strcpy(name, DRA[handle].name);
    strcpy(filename, DRA[handle].fname);

    return(ELIO_OK);
}


/**
 * PRINT OUT INTERNAL PARAMETERS OF DRA
 */
void FATR dra_print_internals_(Integer *d_a)
{
    Integer i;
    Integer *dims, *chunks;
    Integer handle = *d_a + DRA_OFFSET;
    Integer ndim = DRA[handle].ndim;
    Integer me = pnga_nodeid();
    dims = DRA[handle].dims;
    chunks = DRA[handle].chunk;
    if (me == 0) {
        printf("Internal Data for DRA: %s\n",DRA[handle].name);
        printf("  DRA Metafile Name: %s\n",DRA[handle].fname);
        switch(DRA[handle].type){
            case C_DBL:
                printf("  DRA data type is DOUBLE PRECISION\n");
                break;
            case C_FLOAT:
                printf("  DRA data type is SINGLE PRECISION\n");
                break;
            case C_INT:
                printf("  DRA data type is INTEGER\n");
                break;
            case C_DCPL:
                printf("  DRA data type is DOUBLE COMPLEX\n");
                break;
            case C_SCPL:
                printf("  DRA data type is SINGLE COMPLEX\n");
                break;
            case C_LONG:
                printf("  DRA data type is LONG INTEGER\n");
                break;
            default:
                printf("  DRA data type is UNKNOWN\n");
                break;
        }
        switch(DRA[handle].mode) {
            case DRA_RW:
                printf("  DRA access permisions are READ/WRITE\n");
                break;
            case DRA_W:
                printf("  DRA access permisions are WRITE ONLY\n");
                break;
            case DRA_R:
                printf("  DRA access permisions are READ ONLY\n");
                break;
            default:
                printf("  DRA access permisions are UNKNOWN\n");
                break;
        }
        printf("  Dimension of DRA: %d\n",(int)ndim);
        printf("  Dimensions of DRA:\n");
        for (i=0; i<ndim; i++) {
            printf("    Dimension in direction [%d]: %d\n",(int)(i+1),
                    (int)dims[i]);
        }
        printf("  Chunk dimensions of DRA:\n");
        for (i=0; i<ndim; i++) {
            printf("    Chunk dimension in direction [%d]: %d\n",(int)(i+1),
                    (int)chunks[i]);
        }
        if (DRA[handle].actv) {
            printf("  DRA is currently active\n");
        } else {
            printf("  DRA is not currently active\n");
        }
        if (DRA[handle].indep) {
            printf("  DRA is using independent files\n");
        } else {
            printf("  DRA is using shared files\n");
        }
        printf("  Number files used for DRA: %d\n",(int)DRA[handle].numfiles);
        printf("  Number IO processors used for DRA: %d\n",
                (int)DRA[handle].ioprocs);
    }
}


/**
 * SET DEFAULT CONFIGURATION FOR HANDLING DRAs STORED ON OPEN FILE SYSTEMS
 */
void FATR dra_set_default_config_(Integer *numfiles, Integer *numioprocs)
{
    Integer number_of_files, io_procs;
    dai_set_config(*numfiles, *numioprocs, &number_of_files, &io_procs);
    _dra_number_of_files = number_of_files;
    _dra_io_procs = io_procs;
}


/**
 * SET DEBUG FLAG FOR DRA OPERATIONS TO TRUE OR FALSE
 */
void FATR dra_set_debug_(logical *flag)
{
    if (*flag) {
        dra_debug_flag = TRUE;
    } else {
        dra_debug_flag = FALSE;
    }
}
