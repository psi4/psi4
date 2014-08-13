/**
 * @file capi.c
 *
 * Implements the C interface.
 * These calls forward to the (possibly) weak symbols of the internal
 * implementations.
 */
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#include <assert.h>

#include "armci.h"
#include "ga.h"
#include "globalp.h"
#include "ga-papi.h"
#include "matmul.h"

#if ENABLE_PROFILING
#   include "ga-wapi.h"
#else
#   include "ga-wapidefs.h"
#endif

int *_ga_argc=NULL;
char ***_ga_argv=NULL;
int _ga_initialize_args=0;
int _ga_initialize_c=0;

short int _ga_irreg_flag = 0;

static Integer* copy_map(int block[], int block_ndim, int map[]);
static Integer* copy_map64(int64_t block[], int block_ndim, int64_t map[]);

#ifdef USE_FAPI
#  define COPYC2F(carr, farr, n){\
   int i; for(i=0; i< (n); i++)(farr)[i]=(Integer)(carr)[i];} 
#  define COPYF2C(farr, carr, n){\
   int i; for(i=0; i< (n); i++)(carr)[i]=(int)(farr)[i];} 
#  define COPYF2C_64(farr, carr, n){\
   int i; for(i=0; i< (n); i++)(carr)[i]=(int64_t)(farr)[i];} 
#  define COPYINDEX_F2C     COPYF2C
#  define COPYINDEX_F2C_64  COPYF2C_64
#else
#  define COPYC2F(carr, farr, n){\
   int i; for(i=0; i< (n); i++)(farr)[n-i-1]=(Integer)(carr)[i];} 
#  define COPYF2C(farr, carr, n){\
   int i; for(i=0; i< (n); i++)(carr)[n-i-1]=(int)(farr)[i];} 
#  define COPYF2C_64(farr, carr, n){\
   int i; for(i=0; i< (n); i++)(carr)[n-i-1]=(int64_t)(farr)[i];} 
#  define COPYINDEX_C2F(carr, farr, n){\
   int i; for(i=0; i< (n); i++)(farr)[n-i-1]=(Integer)(carr)[i]+1;}
#  define COPYINDEX_F2C(farr, carr, n){\
   int i; for(i=0; i< (n); i++)(carr)[n-i-1]=(int)(farr)[i] -1;}
#  define COPYINDEX_F2C_64(farr, carr, n){\
   int i; for(i=0; i< (n); i++)(carr)[n-i-1]=(int64_t)(farr)[i] -1;}
#define BASE_0
#endif

#define COPY(CAST,src,dst,n) {\
   int i; for(i=0; i< (n); i++)(dst)[i]=(CAST)(src)[i];} 
#define COPY_INC(CAST,src,dst,n) {\
   int i; for(i=0; i< (n); i++)(dst)[i]=(CAST)(src)[i]+1;} 
#define COPY_DEC(CAST,src,dst,n) {\
   int i; for(i=0; i< (n); i++)(dst)[i]=(CAST)(src)[i]-1;} 

int GA_Uses_fapi(void)
{
#ifdef USE_FAPI
return 1;
#else
return 0;
#endif
}


void GA_Initialize_ltd(size_t limit)
{
  Integer lim = (Integer)limit;
  _ga_initialize_c = 1;
  wnga_initialize_ltd(lim);
}

void NGA_Initialize_ltd(size_t limit)
{
  Integer lim = (Integer)limit;
  _ga_initialize_c = 1;
  wnga_initialize_ltd(lim);
}

void GA_Initialize_args(int *argc, char ***argv)
{
  _ga_argc = argc;
  _ga_argv = argv;
  _ga_initialize_c = 1;
  _ga_initialize_args = 1;

  wnga_initialize();
}

void GA_Initialize()
{
  _ga_initialize_c = 1;
  wnga_initialize();
}

void NGA_Initialize()
{
  _ga_initialize_c = 1;
  wnga_initialize();
}

int GA_Initialized()
{
    return wnga_initialized();
}

int NGA_Initialized()
{
    return wnga_initialized();
}

void GA_Terminate() 
{
    wnga_terminate();
    
    _ga_argc = NULL;
    _ga_argv = NULL;
    _ga_initialize_args = 0;
    _ga_initialize_c = 0;
}

void NGA_Terminate() 
{
    wnga_terminate();
    
    _ga_argc = NULL;
    _ga_argv = NULL;
    _ga_initialize_args = 0;
    _ga_initialize_c = 0;
}

int NGA_Create(int type, int ndim, int dims[], char *name, int *chunk)
{
    Integer *ptr, g_a; 
    logical st;
    Integer _ga_work[MAXDIM];
    Integer _ga_dims[MAXDIM];
    if(ndim>MAXDIM)return 0;

    COPYC2F(dims,_ga_dims, ndim);
    if(!chunk)ptr=(Integer*)0;  
    else {
         COPYC2F(chunk,_ga_work, ndim);
         ptr = _ga_work;
    }
    st = wnga_create((Integer)type, (Integer)ndim, _ga_dims, name, ptr, &g_a);
    if(st==TRUE) return (int) g_a;
    else return 0;
}

int NGA_Create64(int type, int ndim, int64_t dims[], char *name, int64_t *chunk)
{
    Integer *ptr, g_a; 
    logical st;
    Integer _ga_dims[MAXDIM];
    Integer _ga_work[MAXDIM];
    if(ndim>MAXDIM)return 0;

    COPYC2F(dims,_ga_dims, ndim);
    if(!chunk)ptr=(Integer*)0;  
    else {
         COPYC2F(chunk,_ga_work, ndim);
         ptr = _ga_work;
    }
    st = wnga_create((Integer)type, (Integer)ndim, _ga_dims, name, ptr, &g_a);
    if(st==TRUE) return (int) g_a;
    else return 0;
}

int NGA_Create_config(int type, int ndim, int dims[], char *name, int chunk[],
                      int p_handle)
{
    Integer *ptr, g_a; 
    logical st;
    Integer _ga_dims[MAXDIM];
    Integer _ga_work[MAXDIM];
    if(ndim>MAXDIM)return 0;

    COPYC2F(dims,_ga_dims, ndim);
    if(!chunk)ptr=(Integer*)0;  
    else {
         COPYC2F(chunk,_ga_work, ndim);
         ptr = _ga_work;
    }
    st = wnga_create_config((Integer)type, (Integer)ndim, _ga_dims, name, ptr,
                    (Integer)p_handle, &g_a);
    if(st==TRUE) return (int) g_a;
    else return 0;
}


int NGA_Create_config64(int type, int ndim, int64_t dims[], char *name, int64_t chunk[], int p_handle)
{
    Integer *ptr, g_a; 
    logical st;
    Integer _ga_dims[MAXDIM];
    Integer _ga_work[MAXDIM];
    if(ndim>MAXDIM)return 0;

    COPYC2F(dims,_ga_dims, ndim);
    if(!chunk)ptr=(Integer*)0;  
    else {
         COPYC2F(chunk,_ga_work, ndim);
         ptr = _ga_work;
    }
    st = wnga_create_config((Integer)type, (Integer)ndim, _ga_dims, name, ptr,
                    (Integer)p_handle, &g_a);
    if(st==TRUE) return (int) g_a;
    else return 0;
}

int NGA_Create_irreg(int type,int ndim,int dims[],char *name,int block[],int map[])
{
    Integer g_a;
    logical st;
    Integer _ga_dims[MAXDIM];
    Integer _ga_work[MAXDIM];
    Integer *_ga_map_capi;
    if(ndim>MAXDIM)return 0;

    COPYC2F(dims,_ga_dims, ndim);
    COPYC2F(block,_ga_work, ndim);
    _ga_map_capi = copy_map(block, ndim, map);

    _ga_irreg_flag = 1; /* set this flag=1, to indicate array is irregular */
    st = wnga_create_irreg(type, (Integer)ndim, _ga_dims, name, _ga_map_capi,
            _ga_work, &g_a);
    _ga_irreg_flag = 0; /* unset it after creating the array */

    free(_ga_map_capi);
    if(st==TRUE) return (int) g_a;
    else return 0;
}

int NGA_Create_irreg64(int type,int ndim,int64_t dims[],char *name,int64_t block[],int64_t map[])
{
    Integer g_a;
    logical st;
    Integer _ga_dims[MAXDIM];
    Integer _ga_work[MAXDIM];
    Integer *_ga_map_capi;
    if(ndim>MAXDIM)return 0;

    COPYC2F(dims,_ga_dims, ndim);
    COPYC2F(block,_ga_work, ndim);
    _ga_map_capi = copy_map64(block, ndim, map);

    _ga_irreg_flag = 1; /* set this flag=1, to indicate array is irregular */
    st = wnga_create_irreg(type, (Integer)ndim, _ga_dims, name, _ga_map_capi,
            _ga_work, &g_a);
    _ga_irreg_flag = 0; /* unset it after creating the array */

    free(_ga_map_capi);
    if(st==TRUE) return (int) g_a;
    else return 0;
}

int NGA_Create_irreg_config(int type,int ndim,int dims[],char *name,int block[],
                            int map[], int p_handle)
{
    Integer g_a;
    logical st;
    Integer _ga_dims[MAXDIM];
    Integer _ga_work[MAXDIM];
    Integer *_ga_map_capi;
    if(ndim>MAXDIM)return 0;

    COPYC2F(dims,_ga_dims, ndim);
    COPYC2F(block,_ga_work, ndim);
    _ga_map_capi = copy_map(block, ndim, map);

    _ga_irreg_flag = 1; /* set this flag=1, to indicate array is irregular */
    st = wnga_create_irreg_config(type, (Integer)ndim, _ga_dims, name,
            _ga_map_capi, _ga_work, (Integer)p_handle, &g_a);
    _ga_irreg_flag = 0; /* unset it, after creating array */

    free(_ga_map_capi);
    if(st==TRUE) return (int) g_a;
    else return 0;
}

int NGA_Create_irreg_config64(int type,int ndim,int64_t dims[],char *name,int64_t block[], int64_t map[], int p_handle)
{
    Integer g_a;
    logical st;
    Integer _ga_dims[MAXDIM];
    Integer _ga_work[MAXDIM];
    Integer *_ga_map_capi;
    if(ndim>MAXDIM)return 0;

    COPYC2F(dims,_ga_dims, ndim);
    COPYC2F(block,_ga_work, ndim);
    _ga_map_capi = copy_map64(block, ndim, map);

    _ga_irreg_flag = 1; /* set this flag=1, to indicate array is irregular */
    st = wnga_create_irreg_config(type, (Integer)ndim, _ga_dims, name,
            _ga_map_capi, _ga_work, (Integer)p_handle, &g_a);
    _ga_irreg_flag = 0; /* unset it, after creating array */

    free(_ga_map_capi);
    if(st==TRUE) return (int) g_a;
    else return 0;
}

int NGA_Create_ghosts_irreg(int type,int ndim,int dims[],int width[],char *name, int block[],int map[])
{
    Integer g_a;
    logical st;
    Integer _ga_dims[MAXDIM];
    Integer _ga_work[MAXDIM];
    Integer _ga_width[MAXDIM];
    Integer *_ga_map_capi;
    if(ndim>MAXDIM)return 0;

    COPYC2F(dims,_ga_dims, ndim);
    COPYC2F(block,_ga_work, ndim);
    COPYC2F(width,_ga_width, ndim);
    _ga_map_capi = copy_map(block, ndim, map);

    _ga_irreg_flag = 1; /* set this flag=1, to indicate array is irregular */
    st = wnga_create_ghosts_irreg(type, (Integer)ndim, _ga_dims, _ga_width,
            name, _ga_map_capi, _ga_work, &g_a);
    _ga_irreg_flag = 0; /* unset it, after creating array */ 

    free(_ga_map_capi);
    if(st==TRUE) return (int) g_a;
    else return 0;
}

int NGA_Create_ghosts_irreg64(int type,int ndim,int64_t dims[],int64_t width[],char *name, int64_t block[],int64_t map[])
{
    Integer g_a;
    logical st;
    Integer _ga_dims[MAXDIM];
    Integer _ga_work[MAXDIM];
    Integer _ga_width[MAXDIM];
    Integer *_ga_map_capi;
    if(ndim>MAXDIM)return 0;

    COPYC2F(dims,_ga_dims, ndim);
    COPYC2F(block,_ga_work, ndim);
    COPYC2F(width,_ga_width, ndim);
    _ga_map_capi = copy_map64(block, ndim, map);

    _ga_irreg_flag = 1; /* set this flag=1, to indicate array is irregular */
    st = wnga_create_ghosts_irreg(type, (Integer)ndim, _ga_dims, _ga_width,
            name, _ga_map_capi, _ga_work, &g_a);
    _ga_irreg_flag = 0; /* unset it, after creating array */ 
     
    free(_ga_map_capi);
    if(st==TRUE) return (int) g_a;
    else return 0;
}

int NGA_Create_ghosts_irreg_config(int type, int ndim, int dims[], int width[], char *name, int block[], int map[], int p_handle)
{
    Integer g_a;
    logical st;
    Integer _ga_dims[MAXDIM];
    Integer _ga_work[MAXDIM];
    Integer _ga_width[MAXDIM];
    Integer *_ga_map_capi;
    if(ndim>MAXDIM)return 0;

    COPYC2F(dims,_ga_dims, ndim);
    COPYC2F(block,_ga_work, ndim);
    COPYC2F(width,_ga_width, ndim);
    _ga_map_capi = copy_map(block, ndim, map);

    _ga_irreg_flag = 1; /* set this flag=1, to indicate array is irregular */
    st = wnga_create_ghosts_irreg_config(type, (Integer)ndim, _ga_dims,
					 _ga_width, name, _ga_map_capi, _ga_work, 
					 (Integer)p_handle, &g_a);
    _ga_irreg_flag = 0; /* unset it, after creating array */ 

    free(_ga_map_capi);
    if(st==TRUE) return (int) g_a;
    else return 0;
}

int NGA_Create_ghosts_irreg_config64(int type, int ndim, int64_t dims[], int64_t width[], char *name, int64_t block[], int64_t map[], int p_handle)
{
    Integer g_a;
    logical st;
    Integer _ga_dims[MAXDIM];
    Integer _ga_work[MAXDIM];
    Integer _ga_width[MAXDIM];
    Integer *_ga_map_capi;
    if(ndim>MAXDIM)return 0;

    COPYC2F(dims,_ga_dims, ndim);
    COPYC2F(block,_ga_work, ndim);
    COPYC2F(width,_ga_width, ndim);
    _ga_map_capi = copy_map64(block, ndim, map);

    _ga_irreg_flag = 1; /* set this flag=1, to indicate array is irregular */
    st = wnga_create_ghosts_irreg_config(type, (Integer)ndim, _ga_dims,
					 _ga_width, name, _ga_map_capi, _ga_work, 
					 (Integer)p_handle, &g_a);
    _ga_irreg_flag = 0; /* unset it, after creating array */ 

    free(_ga_map_capi);
    if(st==TRUE) return (int) g_a;
    else return 0;
}

int NGA_Create_ghosts(int type, int ndim,int dims[], int width[], char *name,
    int chunk[])
{
    Integer *ptr, g_a; 
    logical st;
    Integer _ga_dims[MAXDIM];
    Integer _ga_work[MAXDIM];
    Integer _ga_width[MAXDIM];
    if(ndim>MAXDIM)return 0;

    COPYC2F(dims,_ga_dims, ndim);
    COPYC2F(width,_ga_width, ndim);
    if(!chunk)ptr=(Integer*)0;  
    else {
         COPYC2F(chunk,_ga_work, ndim);
         ptr = _ga_work;
    }
    st = wnga_create_ghosts((Integer)type, (Integer)ndim, _ga_dims,
        _ga_width, name, ptr, &g_a);
    if(st==TRUE) return (int) g_a;
    else return 0;
}


int NGA_Create_ghosts64(int type, int ndim, int64_t dims[], int64_t width[], char *name,
    int64_t chunk[])
{
    Integer *ptr, g_a; 
    logical st;
    Integer _ga_dims[MAXDIM];
    Integer _ga_work[MAXDIM];
    Integer _ga_width[MAXDIM];
    if(ndim>MAXDIM)return 0;

    COPYC2F(dims,_ga_dims, ndim);
    COPYC2F(width,_ga_width, ndim);
    if(!chunk)ptr=(Integer*)0;  
    else {
         COPYC2F(chunk,_ga_work, ndim);
         ptr = _ga_work;
    }
    st = wnga_create_ghosts((Integer)type, (Integer)ndim, _ga_dims,
        _ga_width, name, ptr, &g_a);
    if(st==TRUE) return (int) g_a;
    else return 0;
}


int NGA_Create_ghosts_config(int type, int ndim,int dims[], int width[], char *name,
    int chunk[], int p_handle)
{
    Integer *ptr, g_a; 
    logical st;
    Integer _ga_dims[MAXDIM];
    Integer _ga_work[MAXDIM];
    Integer _ga_width[MAXDIM];
    if(ndim>MAXDIM)return 0;

    COPYC2F(dims,_ga_dims, ndim);
    COPYC2F(width,_ga_width, ndim);
    if(!chunk)ptr=(Integer*)0;  
    else {
         COPYC2F(chunk,_ga_work, ndim);
         ptr = _ga_work;
    }
    st = wnga_create_ghosts_config((Integer)type, (Integer)ndim, _ga_dims,
        _ga_width, name, ptr, (Integer)p_handle, &g_a);
    if(st==TRUE) return (int) g_a;
    else return 0;
}

int NGA_Create_ghosts_config64(int type, int ndim,int64_t  dims[], int64_t width[], char *name,
    int64_t chunk[], int p_handle)
{
    Integer *ptr, g_a; 
    logical st;
    Integer _ga_dims[MAXDIM];
    Integer _ga_work[MAXDIM];
    Integer _ga_width[MAXDIM];
    if(ndim>MAXDIM)return 0;

    COPYC2F(dims,_ga_dims, ndim);
    COPYC2F(width,_ga_width, ndim);
    if(!chunk)ptr=(Integer*)0;  
    else {
         COPYC2F(chunk,_ga_work, ndim);
         ptr = _ga_work;
    }
    st = wnga_create_ghosts_config((Integer)type, (Integer)ndim, _ga_dims,
        _ga_width, name, ptr, (Integer)p_handle, &g_a);
    if(st==TRUE) return (int) g_a;
    else return 0;
}

int GA_Create_handle()
{
    Integer g_a;
    g_a = wnga_create_handle();
    return (int)g_a;
}

int NGA_Create_handle()
{
    Integer g_a;
    g_a = wnga_create_handle();
    return (int)g_a;
}

void GA_Set_data(int g_a, int ndim, int dims[], int type)
{
    Integer aa, nndim, ttype;
    Integer _ga_dims[MAXDIM];
    COPYC2F(dims,_ga_dims, ndim);
    aa = (Integer)g_a;
    nndim = (Integer)ndim;
    ttype = (Integer)type;
    wnga_set_data(aa, nndim, _ga_dims, ttype);
}

void GA_Set_data64(int g_a, int ndim, int64_t dims[], int type)
{
    Integer aa, nndim, ttype;
    Integer _ga_dims[MAXDIM];
    COPYC2F(dims,_ga_dims, ndim);
    aa = (Integer)g_a;
    nndim = (Integer)ndim;
    ttype = (Integer)type;
    wnga_set_data(aa, nndim, _ga_dims, ttype);
}

void NGA_Set_data(int g_a, int ndim, int dims[], int type)
{
    Integer aa, nndim, ttype;
    Integer _ga_dims[MAXDIM];
    COPYC2F(dims,_ga_dims, ndim);
    aa = (Integer)g_a;
    nndim = (Integer)ndim;
    ttype = (Integer)type;
    wnga_set_data(aa, nndim, _ga_dims, ttype);
}

void NGA_Set_data64(int g_a, int ndim, int64_t dims[], int type)
{
    Integer aa, nndim, ttype;
    Integer _ga_dims[MAXDIM];
    COPYC2F(dims,_ga_dims, ndim);
    aa = (Integer)g_a;
    nndim = (Integer)ndim;
    ttype = (Integer)type;
    wnga_set_data(aa, nndim, _ga_dims, ttype);
}

void GA_Set_chunk(int g_a, int chunk[])
{
    Integer aa, *ptr, ndim;
    Integer _ga_work[MAXDIM];
    aa = (Integer)g_a;
    ndim = wnga_get_dimension(aa);
    if(!chunk)ptr=(Integer*)0;  
    else {
      COPYC2F(chunk,_ga_work, ndim);
      ptr = _ga_work;
    }
    wnga_set_chunk(aa, ptr);
}

void GA_Set_chunk64(int g_a, int64_t chunk[])
{
    Integer aa, *ptr, ndim;
    Integer _ga_work[MAXDIM];
    aa = (Integer)g_a;
    ndim = wnga_get_dimension(aa);
    if(!chunk)ptr=(Integer*)0;  
    else {
      COPYC2F(chunk,_ga_work, ndim);
      ptr = _ga_work;
    }
    wnga_set_chunk(aa, ptr);
}

void NGA_Set_chunk(int g_a, int chunk[])
{
    Integer aa, *ptr, ndim;
    Integer _ga_work[MAXDIM];
    aa = (Integer)g_a;
    ndim = wnga_get_dimension(aa);
    if(!chunk)ptr=(Integer*)0;  
    else {
      COPYC2F(chunk,_ga_work, ndim);
      ptr = _ga_work;
    }
    wnga_set_chunk(aa, ptr);
}

void NGA_Set_chunk64(int g_a, int64_t chunk[])
{
    Integer aa, *ptr, ndim;
    Integer _ga_work[MAXDIM];
    aa = (Integer)g_a;
    ndim = wnga_get_dimension(aa);
    if(!chunk)ptr=(Integer*)0;  
    else {
      COPYC2F(chunk,_ga_work, ndim);
      ptr = _ga_work;
    }
    wnga_set_chunk(aa, ptr);
}

void GA_Set_array_name(int g_a, char *name)
{
    Integer aa;
    aa = (Integer)g_a;
    wnga_set_array_name(aa, name);
}

void NGA_Set_array_name(int g_a, char *name)
{
    Integer aa;
    aa = (Integer)g_a;
    wnga_set_array_name(aa, name);
}

void GA_Set_pgroup(int g_a, int p_handle)
{
  Integer aa, pp;
  aa = (Integer)g_a;
  pp = (Integer)p_handle;
  wnga_set_pgroup(aa, pp);
}

void NGA_Set_pgroup(int g_a, int p_handle)
{
  Integer aa, pp;
  aa = (Integer)g_a;
  pp = (Integer)p_handle;
  wnga_set_pgroup(aa, pp);
}

void GA_Set_block_cyclic(int g_a, int dims[])
{
    Integer aa, ndim;
    Integer _ga_dims[MAXDIM];
    aa = (Integer)g_a;
    ndim = wnga_get_dimension(aa);
    COPYC2F(dims,_ga_dims, ndim);
    wnga_set_block_cyclic(aa, _ga_dims);
}

void NGA_Set_block_cyclic(int g_a, int dims[])
{
    Integer aa, ndim;
    Integer _ga_dims[MAXDIM];
    aa = (Integer)g_a;
    ndim = wnga_get_dimension(aa);
    COPYC2F(dims,_ga_dims, ndim);
    wnga_set_block_cyclic(aa, _ga_dims);
}

void GA_Set_restricted(int g_a, int list[], int size)
{
    Integer aa;
    Integer asize = (Integer)size;
    int i;
    Integer *_ga_map_capi;
    aa = (Integer)g_a;
    _ga_map_capi = (Integer*)malloc(size * sizeof(Integer));
    for (i=0; i<size; i++)
       _ga_map_capi[i] = (Integer)list[i];
    wnga_set_restricted(aa,_ga_map_capi,asize);
    free(_ga_map_capi);
}

void NGA_Set_restricted(int g_a, int list[], int size)
{
    Integer aa;
    Integer asize = (Integer)size;
    int i;
    Integer *_ga_map_capi;
    aa = (Integer)g_a;
    _ga_map_capi = (Integer*)malloc(size * sizeof(Integer));
    for (i=0; i<size; i++)
       _ga_map_capi[i] = (Integer)list[i];
    wnga_set_restricted(aa,_ga_map_capi,asize);
    free(_ga_map_capi);
}

void GA_Set_restricted_range(int g_a, int lo_proc, int hi_proc)
{
    Integer aa, lo, hi;
    aa = (Integer)g_a;
    lo = (Integer)lo_proc;
    hi = (Integer)hi_proc;
    wnga_set_restricted_range(aa,lo,hi);
}

void NGA_Set_restricted_range(int g_a, int lo_proc, int hi_proc)
{
    Integer aa, lo, hi;
    aa = (Integer)g_a;
    lo = (Integer)lo_proc;
    hi = (Integer)hi_proc;
    wnga_set_restricted_range(aa,lo,hi);
}

int GA_Total_blocks(int g_a)
{
    Integer aa;
    aa = (Integer)g_a;
    return (int)wnga_total_blocks(aa);
}

int NGA_Total_blocks(int g_a)
{
    Integer aa;
    aa = (Integer)g_a;
    return (int)wnga_total_blocks(aa);
}

void GA_Get_proc_index(int g_a, int iproc, int index[])
{
     Integer aa, proc, ndim;
     Integer _ga_work[MAXDIM];
     aa = (Integer)g_a;
     proc = (Integer)iproc;
     ndim = wnga_get_dimension(aa);
     wnga_get_proc_index(aa, proc, _ga_work);
     COPYF2C(_ga_work,index, ndim);
}

void NGA_Get_proc_index(int g_a, int iproc, int index[])
{
     Integer aa, proc, ndim;
     Integer _ga_work[MAXDIM];
     aa = (Integer)g_a;
     proc = (Integer)iproc;
     ndim = wnga_get_dimension(aa);
     wnga_get_proc_index(aa, proc, _ga_work);
     COPYF2C(_ga_work,index, ndim);
}

void GA_Get_block_info(int g_a, int num_blocks[], int block_dims[])
{
     Integer aa, ndim;
     Integer _ga_work[MAXDIM], _ga_lo[MAXDIM];
     aa = (Integer)g_a;
     ndim = wnga_get_dimension(aa);
     wnga_get_block_info(aa, _ga_work, _ga_lo);
     COPYF2C(_ga_work,num_blocks, ndim);
     COPYF2C(_ga_lo,block_dims, ndim);
}

void NGA_Get_block_info(int g_a, int num_blocks[], int block_dims[])
{
     Integer aa, ndim;
     Integer _ga_work[MAXDIM], _ga_lo[MAXDIM];
     aa = (Integer)g_a;
     ndim = wnga_get_dimension(aa);
     wnga_get_block_info(aa, _ga_work, _ga_lo);
     COPYF2C(_ga_work,num_blocks, ndim);
     COPYF2C(_ga_lo,block_dims, ndim);
}

int GA_Uses_proc_grid(int g_a)
{
     Integer aa = (Integer)g_a;
     return (int)wnga_uses_proc_grid(aa);
}

int NGA_Uses_proc_grid(int g_a)
{
     Integer aa = (Integer)g_a;
     return (int)wnga_uses_proc_grid(aa);
}

int GA_Valid_handle(int g_a)
{
     Integer aa = (Integer)g_a;
     return (int)wnga_valid_handle(aa);
}

int NGA_Valid_handle(int g_a)
{
     Integer aa = (Integer)g_a;
     return (int)wnga_valid_handle(aa);
}

int GA_Verify_handle(int g_a)
{
     Integer aa = (Integer)g_a;
     return (int)wnga_valid_handle(aa);
}

int NGA_Verify_handle(int g_a)
{
     Integer aa = (Integer)g_a;
     return (int)wnga_valid_handle(aa);
}

void GA_Set_block_cyclic_proc_grid(int g_a, int block[], int proc_grid[])
{
    Integer aa, ndim;
    Integer _ga_dims[MAXDIM];
    Integer _ga_lo[MAXDIM];
    aa = (Integer)g_a;
    ndim = wnga_get_dimension(aa);
    COPYC2F(block,_ga_dims, ndim);
    COPYC2F(proc_grid,_ga_lo, ndim);
    wnga_set_block_cyclic_proc_grid(aa, _ga_dims, _ga_lo);
}

void NGA_Set_block_cyclic_proc_grid(int g_a, int block[], int proc_grid[])
{
    Integer aa, ndim;
    Integer _block[MAXDIM];
    Integer _proc_grid[MAXDIM];
    aa = (Integer)g_a;
    ndim = wnga_get_dimension(aa);
    COPYC2F(block,_block, ndim);
    COPYC2F(proc_grid, _proc_grid, ndim);
    wnga_set_block_cyclic_proc_grid(aa, _block, _proc_grid);
}

int GA_Get_pgroup(int g_a)
{
    Integer aa;
    aa = (Integer)g_a;
    return (int)wnga_get_pgroup(aa);
}

int NGA_Get_pgroup(int g_a)
{
    Integer aa;
    aa = (Integer)g_a;
    return (int)wnga_get_pgroup(aa);
}

int GA_Get_pgroup_size(int grp_id)
{
    Integer aa;
    aa = (Integer)grp_id;
    return (int)wnga_get_pgroup_size(aa);
}

int NGA_Get_pgroup_size(int grp_id)
{
    Integer aa;
    aa = (Integer)grp_id;
    return (int)wnga_get_pgroup_size(aa);
}

void GA_Set_ghosts(int g_a, int width[])
{
    Integer aa, *ptr, ndim;
    Integer _ga_work[MAXDIM];
    aa = (Integer)g_a;
    ndim = wnga_get_dimension(aa);
    if(!width)ptr=(Integer*)0;  
    else {
      COPYC2F(width,_ga_work, ndim);
      ptr = _ga_work;
    }
    wnga_set_ghosts(aa, ptr);
}

void GA_Set_ghosts64(int g_a, int64_t width[])
{
    Integer aa, *ptr, ndim;
    Integer _ga_work[MAXDIM];
    aa = (Integer)g_a;
    ndim = wnga_get_dimension(aa);
    if(!width)ptr=(Integer*)0;  
    else {
      COPYC2F(width,_ga_work, ndim);
      ptr = _ga_work;
    }
    wnga_set_ghosts(aa, ptr);
}

void NGA_Set_ghosts(int g_a, int width[])
{
    Integer aa, *ptr, ndim;
    Integer _ga_work[MAXDIM];
    aa = (Integer)g_a;
    ndim = wnga_get_dimension(aa);
    if(!width)ptr=(Integer*)0;  
    else {
      COPYC2F(width,_ga_work, ndim);
      ptr = _ga_work;
    }
    wnga_set_ghosts(aa, ptr);
}

void NGA_Set_ghosts64(int g_a, int64_t width[])
{
    Integer aa, *ptr, ndim;
    Integer _ga_work[MAXDIM];
    aa = (Integer)g_a;
    ndim = wnga_get_dimension(aa);
    if(!width)ptr=(Integer*)0;  
    else {
      COPYC2F(width,_ga_work, ndim);
      ptr = _ga_work;
    }
    wnga_set_ghosts(aa, ptr);
}

void GA_Set_irreg_distr(int g_a, int map[], int block[])
{
    Integer aa, ndim;
    Integer _ga_work[MAXDIM];
    Integer *_ga_map_capi;

    aa = (Integer)g_a;
    ndim = wnga_get_dimension(aa);
    COPYC2F(block,_ga_work, ndim);
    _ga_map_capi = copy_map(block, (int)ndim, map);

    wnga_set_irreg_distr(aa, _ga_map_capi, _ga_work);
    free(_ga_map_capi);
}

void GA_Set_irreg_distr64(int g_a, int64_t map[], int64_t block[])
{
    Integer aa, ndim;
    Integer _ga_work[MAXDIM];
    Integer *_ga_map_capi;

    aa = (Integer)g_a;
    ndim = wnga_get_dimension(aa);
    COPYC2F(block,_ga_work, ndim);
    _ga_map_capi = copy_map64(block, (int)ndim, map);

    wnga_set_irreg_distr(aa, _ga_map_capi, _ga_work);
    free(_ga_map_capi);
}

void NGA_Set_irreg_distr(int g_a, int map[], int block[])
{
    Integer aa, ndim;
    Integer _ga_work[MAXDIM];
    Integer *_ga_map_capi;

    aa = (Integer)g_a;
    ndim = wnga_get_dimension(aa);
    COPYC2F(block,_ga_work, ndim);
    _ga_map_capi = copy_map(block, (int)ndim, map);

    wnga_set_irreg_distr(aa, _ga_map_capi, _ga_work);
    free(_ga_map_capi);
}

void NGA_Set_irreg_distr64(int g_a, int64_t map[], int64_t block[])
{
    Integer aa, ndim;
    Integer _ga_work[MAXDIM];
    Integer *_ga_map_capi;

    aa = (Integer)g_a;
    ndim = wnga_get_dimension(aa);
    COPYC2F(block,_ga_work, ndim);
    _ga_map_capi = copy_map64(block, (int)ndim, map);

    wnga_set_irreg_distr(aa, _ga_map_capi, _ga_work);
    free(_ga_map_capi);
}

void GA_Set_irreg_flag(int g_a, int flag)
{
  Integer aa;
  logical fflag;
  aa = (Integer)g_a;
  fflag = (logical)flag;
  wnga_set_irreg_flag(aa, fflag);
}

void NGA_Set_irreg_flag(int g_a, int flag)
{
  Integer aa;
  logical fflag;
  aa = (Integer)g_a;
  fflag = (logical)flag;
  wnga_set_irreg_flag(aa, fflag);
}

void GA_Set_ghost_corner_flag(int g_a, int flag)
{
  Integer aa;
  logical fflag;
  aa = (Integer)g_a;
  fflag = (logical)flag;
  wnga_set_ghost_corner_flag(aa, fflag);
}

void NGA_Set_ghost_corner_flag(int g_a, int flag)
{
  Integer aa;
  logical fflag;
  aa = (Integer)g_a;
  fflag = (logical)flag;
  wnga_set_ghost_corner_flag(aa, fflag);
}

int GA_Get_dimension(int g_a)
{
  Integer aa;
  aa = (Integer)g_a;
  return (int)wnga_get_dimension(aa);
}

int NGA_Get_dimension(int g_a)
{
  Integer aa;
  aa = (Integer)g_a;
  return (int)wnga_get_dimension(aa);
}

int GA_Allocate(int g_a)
{
  Integer aa;
  aa = (Integer)g_a;
  return (int)wnga_allocate(aa);
}

int NGA_Allocate(int g_a)
{
  Integer aa;
  aa = (Integer)g_a;
  return (int)wnga_allocate(aa);
}

int GA_Pgroup_nodeid(int grp_id)
{
    Integer agrp_id = (Integer)grp_id;
    return (int)wnga_pgroup_nodeid(agrp_id);
}

int NGA_Pgroup_nodeid(int grp_id)
{
    Integer agrp_id = (Integer)grp_id;
    return (int)wnga_pgroup_nodeid(agrp_id);
}

int GA_Pgroup_nnodes(int grp_id)
{
    Integer agrp_id = (Integer)grp_id;
    return (int)wnga_pgroup_nnodes(agrp_id);
}

int NGA_Pgroup_nnodes(int grp_id)
{
    Integer agrp_id = (Integer)grp_id;
    return (int)wnga_pgroup_nnodes(agrp_id);
}

int GA_Pgroup_create(int *list, int count)
{
    Integer acount = (Integer)count;
    int i;
    int grp_id;
    Integer *_ga_map_capi;
    _ga_map_capi = (Integer*)malloc(count * sizeof(Integer));
    for (i=0; i<count; i++)
       _ga_map_capi[i] = (Integer)list[i];
    grp_id = (int)wnga_pgroup_create(_ga_map_capi,acount);
    free(_ga_map_capi);
    return grp_id;
}

int NGA_Pgroup_create(int *list, int count)
{
    Integer acount = (Integer)count;
    int i;
    int grp_id;
    Integer *_ga_map_capi;
    _ga_map_capi = (Integer*)malloc(count * sizeof(Integer));
    for (i=0; i<count; i++)
       _ga_map_capi[i] = (Integer)list[i];
    grp_id = (int)wnga_pgroup_create(_ga_map_capi,acount);
    free(_ga_map_capi);
    return grp_id;
}

int GA_Pgroup_destroy(int grp)
{
    Integer grp_id = (Integer)grp;
    return (int)wnga_pgroup_destroy(grp_id);
}

int NGA_Pgroup_destroy(int grp)
{
    Integer grp_id = (Integer)grp;
    return (int)wnga_pgroup_destroy(grp_id);
}

int GA_Pgroup_split(int grp_id, int num_group)
{
    Integer anum = (Integer)num_group;
    Integer grp  = (Integer)grp_id;
    return (int)wnga_pgroup_split(grp, anum);
}

int NGA_Pgroup_split(int grp_id, int num_group)
{
    Integer anum = (Integer)num_group;
    Integer grp  = (Integer)grp_id;
    return (int)wnga_pgroup_split(grp, anum);
}

int GA_Pgroup_split_irreg(int grp_id, int color)
{
    Integer acolor = (Integer)color;
    Integer grp  = (Integer)grp_id;
    return (int)wnga_pgroup_split_irreg(grp, acolor);
}

int NGA_Pgroup_split_irreg(int grp_id, int color)
{
    Integer acolor = (Integer)color;
    Integer grp  = (Integer)grp_id;
    return (int)wnga_pgroup_split_irreg(grp, acolor);
}

void GA_Update_ghosts(int g_a)
{
    Integer a=(Integer)g_a;
    wnga_update_ghosts(a);
}

void NGA_Update_ghosts(int g_a)
{
    Integer a=(Integer)g_a;
    wnga_update_ghosts(a);
}

void NGA_Update_ghosts_nb(int g_a, ga_nbhdl_t *nbhandle)
{
    Integer a=(Integer)g_a;
    wnga_update_ghosts_nb(a,nbhandle);
}

void GA_Merge_mirrored(int g_a)
{
    Integer a=(Integer)g_a;
    wnga_merge_mirrored(a);
}

void NGA_Merge_mirrored(int g_a)
{
    Integer a=(Integer)g_a;
    wnga_merge_mirrored(a);
}

void GA_Nblock(int g_a, int *nblock)
{
    Integer aa, ndim;
    Integer _ga_work[MAXDIM];
    aa = (Integer)g_a;
    wnga_nblock(aa, _ga_work);
    ndim = wnga_get_dimension(aa);
    COPYF2C(_ga_work,nblock,ndim);
}

void NGA_Nblock(int g_a, int *nblock)
{
    Integer aa, ndim;
    Integer _ga_work[MAXDIM];
    aa = (Integer)g_a;
    wnga_nblock(aa, _ga_work);
    ndim = wnga_get_dimension(aa);
    COPYF2C(_ga_work,nblock,ndim);
}

int GA_Is_mirrored(int g_a)
{
    Integer a=(Integer)g_a;
    return (int)wnga_is_mirrored(a);
}

int NGA_Is_mirrored(int g_a)
{
    Integer a=(Integer)g_a;
    return (int)wnga_is_mirrored(a);
}

void GA_List_nodeid(int *nlist, int nprocs)
{
  Integer i, procs;
  Integer *list;
  procs = (Integer)(nprocs);
  list = malloc(procs*sizeof(int));
  wnga_list_nodeid(list, procs);
  for (i=0; i<procs; i++) {
    nlist[i] = (int)list[i];
  }
  free(list);
}

void NGA_List_nodeid(int *nlist, int nprocs)
{
  Integer i, procs;
  Integer *list;
  procs = (Integer)(nprocs);
  list = malloc(procs*sizeof(int));
  wnga_list_nodeid(list, procs);
  for (i=0; i<procs; i++) {
    nlist[i] = (int)list[i];
  }
  free(list);
}

void NGA_Merge_distr_patch(int g_a, int *alo, int *ahi,
                          int g_b, int *blo, int *bhi)
{
    Integer a = (Integer)g_a;
    Integer andim = wnga_ndim(a);

    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);

    wnga_merge_distr_patch(a, _ga_alo, _ga_ahi, b, _ga_blo, _ga_bhi);
}

void NGA_Merge_distr_patch64(int g_a, int64_t *alo, int64_t *ahi,
                             int g_b, int64_t *blo, int64_t *bhi)
{
    Integer a = (Integer)g_a;
    Integer andim = wnga_ndim(a);

    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);

    wnga_merge_distr_patch(a, _ga_alo, _ga_ahi, b, _ga_blo, _ga_bhi);
}

int NGA_Update_ghost_dir(int g_a, int dimension, int dir, int flag)
{
    Integer a=(Integer)g_a;
    Integer idim = (Integer)dimension;
    Integer idir = (Integer)dir;
    logical iflag = (logical)flag;
    logical st;
    idim++;
    st = wnga_update_ghost_dir(a,idim,idir,iflag);
    return (int)st;
}

void NGA_NbGet_ghost_dir(int g_a, int mask[], ga_nbhdl_t* nbhandle)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM];
    COPYINDEX_C2F(mask,_ga_lo, ndim);
    wnga_nbget_ghost_dir(a, _ga_lo,(Integer *)nbhandle);
}

void NGA_NbGet_ghost_dir64(int g_a, int64_t mask[], ga_nbhdl_t* nbhandle)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM];
    COPYINDEX_C2F(mask,_ga_lo, ndim);
    wnga_nbget_ghost_dir(a, _ga_lo,(Integer *)nbhandle);
}

int GA_Has_ghosts(int g_a)
{
    Integer a=(Integer)g_a;
    logical st;
    st = wnga_has_ghosts(a);
    return (int)st;
}

int NGA_Has_ghosts(int g_a)
{
    Integer a=(Integer)g_a;
    logical st;
    st = wnga_has_ghosts(a);
    return (int)st;
}

void NGA_Get_ghost_block(int g_a, int lo[], int hi[], void *buf, int ld[])
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_get_ghost_block(a, _ga_lo, _ga_hi, buf, _ga_work);
}

void NGA_Get_ghost_block64(int g_a, int64_t lo[], int64_t hi[],
    void *buf, int64_t ld[])
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_get_ghost_block(a, _ga_lo, _ga_hi, buf, _ga_work);
}

void GA_Mask_sync(int first, int last)
{
    Integer ifirst = (Integer)first;
    Integer ilast = (Integer)last;
    wnga_mask_sync(ifirst,ilast);
}

void NGA_Mask_sync(int first, int last)
{
    Integer ifirst = (Integer)first;
    Integer ilast = (Integer)last;
    wnga_mask_sync(ifirst,ilast);
}

int GA_Duplicate(int g_a, char* array_name)
{
    logical st;
    Integer a=(Integer)g_a, b;
    st = wnga_duplicate(a, &b, array_name);
    if(st==TRUE) return (int) b;
    else return 0;
}

int NGA_Duplicate(int g_a, char* array_name)
{
    logical st;
    Integer a=(Integer)g_a, b;
    st = wnga_duplicate(a, &b, array_name);
    if(st==TRUE) return (int) b;
    else return 0;
}

void GA_Destroy(int g_a)
{
    logical st;
    Integer a=(Integer)g_a;
    st = wnga_destroy(a);
    if(st==FALSE)GA_Error("GA (c) destroy failed",g_a);
}

void NGA_Destroy(int g_a)
{
    logical st;
    Integer a=(Integer)g_a;
    st = wnga_destroy(a);
    if(st==FALSE)GA_Error("GA (c) destroy failed",g_a);
}

void GA_Set_memory_limit(size_t limit)
{
Integer lim = (Integer)limit;
     wnga_set_memory_limit(lim);
}

void NGA_Set_memory_limit(size_t limit)
{
Integer lim = (Integer)limit;
     wnga_set_memory_limit(lim);
}

void GA_Zero(int g_a)
{
    Integer a=(Integer)g_a;
    wnga_zero(a);
}

void NGA_Zero(int g_a)
{
    Integer a=(Integer)g_a;
    wnga_zero(a);
}

int GA_Pgroup_get_default()
{
    int value = (int)wnga_pgroup_get_default();
    return value;
}

int NGA_Pgroup_get_default()
{
    int value = (int)wnga_pgroup_get_default();
    return value;
}

void GA_Pgroup_set_default(int p_handle)
{
    Integer grp = (Integer)p_handle;
    wnga_pgroup_set_default(grp);
}

void NGA_Pgroup_set_default(int p_handle)
{
    Integer grp = (Integer)p_handle;
    wnga_pgroup_set_default(grp);
}

int GA_Pgroup_get_mirror()
{
    int value = (int)wnga_pgroup_get_mirror();
    return value;
}

int NGA_Pgroup_get_mirror()
{
    int value = (int)wnga_pgroup_get_mirror();
    return value;
}

int GA_Pgroup_get_world()
{
    int value = (int)wnga_pgroup_get_world();
    return value;
}

int NGA_Pgroup_get_world()
{
    int value = (int)wnga_pgroup_get_world();
    return value;
}

int GA_Idot(int g_a, int g_b)
{
    int value;
    Integer a=(Integer)g_a;
    Integer b=(Integer)g_b;
    wnga_dot(C_INT, a, b, &value);
    return value;
}


long GA_Ldot(int g_a, int g_b)
{
    long value;
    Integer a=(Integer)g_a;
    Integer b=(Integer)g_b;
    wnga_dot(C_LONG, a, b, &value);
    return value;
}


long long GA_Lldot(int g_a, int g_b)
{
    long long value;
    Integer a=(Integer)g_a;
    Integer b=(Integer)g_b;
    wnga_dot(C_LONGLONG, a, b, &value);
    return value;
}

     
double GA_Ddot(int g_a, int g_b)
{
    double value;
    Integer a=(Integer)g_a;
    Integer b=(Integer)g_b;
    wnga_dot(C_DBL, a, b, &value);
    return value;
}


DoubleComplex GA_Zdot(int g_a, int g_b)
{
    DoubleComplex value;
    Integer a=(Integer)g_a;
    Integer b=(Integer)g_b;
    wnga_dot(pnga_type_f2c(MT_F_DCPL),a,b,&value);
    return value;
}


SingleComplex GA_Cdot(int g_a, int g_b)
{
    SingleComplex value;
    Integer a=(Integer)g_a;
    Integer b=(Integer)g_b;
    wnga_dot(pnga_type_f2c(MT_F_SCPL),a,b,&value);
    return value;
}


float GA_Fdot(int g_a, int g_b)
{
    float sum;
    Integer a=(Integer)g_a;
    Integer b=(Integer)g_b;
    wnga_dot(C_FLOAT, a, b, &sum);
    return sum;
}    


void GA_Randomize(int g_a, void *value)
{
    Integer a=(Integer)g_a;
    wnga_randomize(a, value);
}

void NGA_Randomize(int g_a, void *value)
{
    Integer a=(Integer)g_a;
    wnga_randomize(a, value);
}

void GA_Fill(int g_a, void *value)
{
    Integer a=(Integer)g_a;
    wnga_fill(a, value);
}

void NGA_Fill(int g_a, void *value)
{
    Integer a=(Integer)g_a;
    wnga_fill(a, value);
}

void GA_Scale(int g_a, void *value)
{
    Integer a=(Integer)g_a;
    wnga_scale(a,value);
}


void GA_Add(void *alpha, int g_a, void* beta, int g_b, int g_c)
{
    Integer a=(Integer)g_a;
    Integer b=(Integer)g_b;
    Integer c=(Integer)g_c;
    wnga_add(alpha, a, beta, b, c);
}


void GA_Copy(int g_a, int g_b)
{
    Integer a=(Integer)g_a;
    Integer b=(Integer)g_b;
    wnga_copy(a, b);
}


void NGA_Get(int g_a, int lo[], int hi[], void* buf, int ld[])
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_get(a, _ga_lo, _ga_hi, buf, _ga_work);
}

void NGA_Get64(int g_a, int64_t lo[], int64_t hi[], void* buf, int64_t ld[])
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_get(a, _ga_lo, _ga_hi, buf, _ga_work);
}

void NGA_NbGet(int g_a, int lo[], int hi[], void* buf, int ld[],
               ga_nbhdl_t* nbhandle)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_nbget(a, _ga_lo, _ga_hi, buf, _ga_work,(Integer *)nbhandle);
}

void NGA_NbGet64(int g_a, int64_t lo[], int64_t hi[], void* buf, int64_t ld[],
               ga_nbhdl_t* nbhandle)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_nbget(a, _ga_lo, _ga_hi, buf, _ga_work,(Integer *)nbhandle);
}

void NGA_Put(int g_a, int lo[], int hi[], void* buf, int ld[])
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_put(a, _ga_lo, _ga_hi, buf, _ga_work);
}    

void NGA_Put64(int g_a, int64_t lo[], int64_t hi[], void* buf, int64_t ld[])
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_put(a, _ga_lo, _ga_hi, buf, _ga_work);
}    

void NGA_NbPut(int g_a, int lo[], int hi[], void* buf, int ld[],
               ga_nbhdl_t* nbhandle)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_nbput(a, _ga_lo, _ga_hi, buf, _ga_work,(Integer *)nbhandle);
}

void NGA_NbPut64(int g_a, int64_t lo[], int64_t hi[], void* buf, int64_t ld[],
                 ga_nbhdl_t* nbhandle)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_nbput(a, _ga_lo, _ga_hi, buf, _ga_work,(Integer *)nbhandle);
}

int NGA_NbTest(ga_nbhdl_t* nbhandle)
{
    return(wnga_nbtest((Integer *)nbhandle));
}

void NGA_NbWait(ga_nbhdl_t* nbhandle)
{
    wnga_nbwait((Integer *)nbhandle);
}

void NGA_Strided_acc(int g_a, int lo[], int hi[], int skip[],
                     void* buf, int ld[], void *alpha)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM], _ga_skip[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    COPYC2F(skip, _ga_skip, ndim);
    wnga_strided_acc(a, _ga_lo, _ga_hi, _ga_skip, buf, _ga_work, alpha);
}    

void NGA_Strided_acc64(int g_a, int64_t lo[], int64_t hi[], int64_t skip[],
                     void* buf, int64_t ld[], void *alpha)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_skip[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    COPYC2F(skip, _ga_skip, ndim);
    wnga_strided_acc(a, _ga_lo, _ga_hi, _ga_skip, buf, _ga_work, alpha);
}

void NGA_Strided_get(int g_a, int lo[], int hi[], int skip[],
                     void* buf, int ld[])
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_skip[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    COPYC2F(skip, _ga_skip, ndim);
    wnga_strided_get(a, _ga_lo, _ga_hi, _ga_skip, buf, _ga_work);
}    

void NGA_Strided_get64(int g_a, int64_t lo[], int64_t hi[], int64_t skip[],
                     void* buf, int64_t ld[])
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_skip[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    COPYC2F(skip, _ga_skip, ndim);
    wnga_strided_get(a, _ga_lo, _ga_hi, _ga_skip, buf, _ga_work);
}

void NGA_Strided_put(int g_a, int lo[], int hi[], int skip[],
                     void* buf, int ld[])
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_skip[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    COPYC2F(skip, _ga_skip, ndim);
    wnga_strided_put(a, _ga_lo, _ga_hi, _ga_skip, buf, _ga_work);
}    

void NGA_Strided_put64(int g_a, int64_t lo[], int64_t hi[], int64_t skip[],
                     void* buf, int64_t ld[])
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_skip[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    COPYC2F(skip, _ga_skip, ndim);
    wnga_strided_put(a, _ga_lo, _ga_hi, _ga_skip, buf, _ga_work);
}

void NGA_Acc(int g_a, int lo[], int hi[], void* buf,int ld[], void* alpha)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_acc(a, _ga_lo, _ga_hi, buf, _ga_work, alpha);
}    

void NGA_Acc64(int g_a, int64_t lo[], int64_t hi[], void* buf, int64_t ld[], void* alpha)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_acc(a, _ga_lo, _ga_hi, buf, _ga_work, alpha);
}    

void NGA_NbAcc(int g_a, int lo[], int hi[], void* buf,int ld[], void* alpha,
               ga_nbhdl_t* nbhandle)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_nbacc(a, _ga_lo,_ga_hi,buf,_ga_work,alpha,(Integer *)nbhandle);
}

void NGA_NbAcc64(int g_a, int64_t lo[], int64_t hi[], void* buf,int64_t ld[], void* alpha,
               ga_nbhdl_t* nbhandle)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_nbacc(a, _ga_lo,_ga_hi,buf,_ga_work,alpha,(Integer *)nbhandle);
}

void NGA_Periodic_get(int g_a, int lo[], int hi[], void* buf, int ld[])
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_periodic(a, _ga_lo, _ga_hi, buf, _ga_work, NULL, PERIODIC_GET);
}

void NGA_Periodic_get64(int g_a, int64_t lo[], int64_t hi[], void* buf, int64_t ld[])
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_periodic(a, _ga_lo, _ga_hi, buf, _ga_work, NULL, PERIODIC_GET);
}

void NGA_Periodic_put(int g_a, int lo[], int hi[], void* buf, int ld[])
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_periodic(a, _ga_lo, _ga_hi, buf, _ga_work, NULL, PERIODIC_PUT);
}    

void NGA_Periodic_put64(int g_a, int64_t lo[], int64_t hi[], void* buf, int64_t ld[])
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_periodic(a, _ga_lo, _ga_hi, buf, _ga_work, NULL, PERIODIC_PUT);
}

void NGA_Periodic_acc(int g_a, int lo[], int hi[], void* buf,int ld[], void* alpha)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_periodic(a, _ga_lo, _ga_hi, buf, _ga_work, alpha, PERIODIC_ACC);
}

void NGA_Periodic_acc64(int g_a, int64_t lo[], int64_t hi[], void* buf,int64_t ld[], void* alpha)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    Integer _ga_work[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    COPYC2F(ld,_ga_work, ndim-1);
    wnga_periodic(a, _ga_lo, _ga_hi, buf, _ga_work, alpha, PERIODIC_ACC);
}

long NGA_Read_inc(int g_a, int subscript[], long inc)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer in=(Integer)inc;
    Integer _ga_lo[MAXDIM];
    COPYINDEX_C2F(subscript, _ga_lo, ndim);
    return (long)wnga_read_inc(a, _ga_lo, in);
}

long NGA_Read_inc64(int g_a, int64_t subscript[], long inc)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer in=(Integer)inc;
    Integer _ga_lo[MAXDIM];
    COPYINDEX_C2F(subscript, _ga_lo, ndim);
    return (long)wnga_read_inc(a, _ga_lo, in);
}

void NGA_Distribution(int g_a, int iproc, int lo[], int hi[])
{
     Integer a=(Integer)g_a;
     Integer p=(Integer)iproc;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
     wnga_distribution(a, p, _ga_lo, _ga_hi);
     COPYINDEX_F2C(_ga_lo,lo, ndim);
     COPYINDEX_F2C(_ga_hi,hi, ndim);
}

void NGA_Distribution64(int g_a, int iproc, int64_t lo[], int64_t hi[])
{
     Integer a=(Integer)g_a;
     Integer p=(Integer)iproc;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
     wnga_distribution(a, p, _ga_lo, _ga_hi);
     COPYINDEX_F2C_64(_ga_lo,lo, ndim);
     COPYINDEX_F2C_64(_ga_hi,hi, ndim);
}

void NGA_Select_elem(int g_a, char* op, void* val, int* index)
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM];
     wnga_select_elem(a, op, val, _ga_lo);
     COPYINDEX_F2C(_ga_lo,index,ndim);
}

void NGA_Select_elem64(int g_a, char* op, void* val, int64_t* index)
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM];
     wnga_select_elem(a, op, val, _ga_lo);
     COPYINDEX_F2C_64(_ga_lo,index,ndim);
}

void GA_Scan_add(int g_a, int g_b, int g_sbit, int lo,
                 int hi, int excl)
{
     Integer a = (Integer)g_a;
     Integer b = (Integer)g_b;
     Integer s = (Integer)g_sbit;
     Integer x = (Integer)excl;
     Integer alo = lo+1;
     Integer ahi = hi+1;
     wnga_scan_add(a, b, s, alo, ahi, x);

}

void GA_Scan_add64(int g_a, int g_b, int g_sbit, int64_t lo,
                 int64_t hi, int excl)
{
     Integer a = (Integer)g_a;
     Integer b = (Integer)g_b;
     Integer s = (Integer)g_sbit;
     Integer x = (Integer)excl;
     Integer alo = lo+1;
     Integer ahi = hi+1;
     wnga_scan_add(a, b, s, alo, ahi, x);

}

void GA_Scan_copy(int g_a, int g_b, int g_sbit, int lo,
                  int hi)
{
     Integer a = (Integer)g_a;
     Integer b = (Integer)g_b;
     Integer s = (Integer)g_sbit;
     Integer alo = lo+1;
     Integer ahi = hi+1;
     wnga_scan_copy(a, b, s, alo, ahi);

}

void GA_Scan_copy64(int g_a, int g_b, int g_sbit, int64_t lo,
                    int64_t hi)
{
     Integer a = (Integer)g_a;
     Integer b = (Integer)g_b;
     Integer s = (Integer)g_sbit;
     Integer alo = lo+1;
     Integer ahi = hi+1;
     wnga_scan_copy(a, b, s, alo, ahi);

}

void GA_Patch_enum(int g_a, int lo, int hi, void *start, void *inc)
{
     Integer a = (Integer)g_a;
     Integer alo = lo+1;
     Integer ahi = hi+1;
     wnga_patch_enum(a, alo, ahi, start, inc);
}

void GA_Patch_enum64(int g_a, int64_t lo, int64_t hi, void *start, void *inc)
{
     Integer a = (Integer)g_a;
     Integer alo = lo+1;
     Integer ahi = hi+1;
     wnga_patch_enum(a, alo, ahi, start, inc);
}

void GA_Pack(int g_src, int g_dest, int g_mask, int lo, int hi, int *icount)
{
     Integer a = (Integer)g_src;
     Integer b = (Integer)g_dest;
     Integer s = (Integer)g_mask;
     Integer icnt;
     Integer alo = lo+1;
     Integer ahi = hi+1;
     wnga_pack(a, b, s, alo, ahi, &icnt); 
     *icount = icnt;
}

void GA_Pack64(int g_src, int g_dest, int g_mask, int64_t lo, int64_t hi, int64_t *icount)
{
     Integer a = (Integer)g_src;
     Integer b = (Integer)g_dest;
     Integer s = (Integer)g_mask;
     Integer icnt;
     Integer alo = lo+1;
     Integer ahi = hi+1;
     wnga_pack(a, b, s, alo, ahi, &icnt); 
     *icount = icnt;
}

void GA_Unpack(int g_src, int g_dest, int g_mask, int lo, int hi, int *icount)
{
     Integer a = (Integer)g_src;
     Integer b = (Integer)g_dest;
     Integer s = (Integer)g_mask;
     Integer icnt;
     Integer alo = lo+1;
     Integer ahi = hi+1;
     wnga_unpack(a, b, s, alo, ahi, &icnt); 
     *icount = icnt;
}

void GA_Unpack64(int g_src, int g_dest, int g_mask, int64_t lo, int64_t hi, int64_t *icount)
{
     Integer a = (Integer)g_src;
     Integer b = (Integer)g_dest;
     Integer s = (Integer)g_mask;
     Integer icnt;
     Integer alo = lo+1;
     Integer ahi = hi+1;
     wnga_unpack(a, b, s, alo, ahi, &icnt); 
     *icount = icnt;
}

int GA_Compare_distr(int g_a, int g_b)
{
    logical st;
    Integer a=(Integer)g_a;
    Integer b=(Integer)g_b;
    st = wnga_compare_distr(a,b);
    if(st == TRUE) return 0;
    else return 1;
}

int NGA_Compare_distr(int g_a, int g_b)
{
    logical st;
    Integer a=(Integer)g_a;
    Integer b=(Integer)g_b;
    st = wnga_compare_distr(a,b);
    if(st == TRUE) return 0;
    else return 1;
}

void NGA_Access(int g_a, int lo[], int hi[], void *ptr, int ld[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
     Integer _ga_work[MAXDIM];
     COPYINDEX_C2F(lo,_ga_lo,ndim);
     COPYINDEX_C2F(hi,_ga_hi,ndim);

     wnga_access_ptr(a,_ga_lo, _ga_hi, ptr, _ga_work);
     COPYF2C(_ga_work,ld, ndim-1);
}

void NGA_Access64(int g_a, int64_t lo[], int64_t hi[], void *ptr, int64_t ld[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
     Integer _ga_work[MAXDIM];
     COPYINDEX_C2F(lo,_ga_lo,ndim);
     COPYINDEX_C2F(hi,_ga_hi,ndim);

     wnga_access_ptr(a,_ga_lo, _ga_hi, ptr, _ga_work);
     COPYF2C_64(_ga_work,ld, ndim-1);
}

void NGA_Access_block(int g_a, int idx, void *ptr, int ld[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer iblock = (Integer)idx;
     Integer _ga_work[MAXDIM];
     wnga_access_block_ptr(a,iblock,ptr,_ga_work);
     COPYF2C(_ga_work,ld, ndim-1);
}

void NGA_Access_block64(int g_a, int64_t idx, void *ptr, int64_t ld[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer iblock = (Integer)idx;
     Integer _ga_work[MAXDIM];
     wnga_access_block_ptr(a,iblock,ptr,_ga_work);
     COPYF2C_64(_ga_work,ld, ndim-1);
}

void NGA_Access_block_grid(int g_a, int index[], void *ptr, int ld[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_work[MAXDIM], _ga_lo[MAXDIM];
     COPYC2F(index, _ga_lo, ndim);
     wnga_access_block_grid_ptr(a,_ga_lo,ptr,_ga_work);
     COPYF2C(_ga_work,ld, ndim-1);
}

void NGA_Access_block_grid64(int g_a, int64_t index[], void *ptr, int64_t ld[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM], _ga_work[MAXDIM];
     COPYC2F(index, _ga_lo, ndim);
     wnga_access_block_grid_ptr(a,_ga_lo,ptr,_ga_work);
     COPYF2C_64(_ga_work,ld, ndim-1);
}

void NGA_Access_block_segment(int g_a, int proc, void *ptr, int *len)
{
     Integer a=(Integer)g_a;
     Integer iblock = (Integer)proc;
     Integer ilen = (Integer)(*len);
     wnga_access_block_segment_ptr(a,iblock,ptr,&ilen);
     *len = (int)ilen;
}

void NGA_Access_block_segment64(int g_a, int proc, void *ptr, int64_t *len)
{
     Integer a=(Integer)g_a;
     Integer iblock = (Integer)proc;
     Integer ilen = (Integer)(*len);
     wnga_access_block_segment_ptr(a,iblock,ptr,&ilen);
     *len = (int64_t)ilen;
}

void NGA_Access_ghosts(int g_a, int dims[], void *ptr, int ld[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM];
     Integer _ga_work[MAXDIM];
     wnga_access_ghost_ptr(a, _ga_lo, ptr, _ga_work);
     COPYF2C(_ga_work,ld, ndim-1);
     COPYF2C(_ga_lo, dims, ndim);
}

void NGA_Access_ghosts64(int g_a, int64_t dims[], void *ptr, int64_t ld[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM];
     Integer _ga_work[MAXDIM];

     wnga_access_ghost_ptr(a, _ga_lo, ptr, _ga_work);
     COPYF2C_64(_ga_work,ld, ndim-1);
     COPYF2C_64(_ga_lo, dims, ndim);
}

void NGA_Access_ghost_element(int g_a, void *ptr, int subscript[], int ld[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _subscript[MAXDIM];
     Integer _ld[MAXDIM];

     COPYINDEX_C2F(subscript, _subscript, ndim);
     COPYINDEX_C2F(ld, _ld, ndim-1);
     wnga_access_ghost_element_ptr(a, ptr, _subscript, _ld);
}

void NGA_Access_ghost_element64(int g_a, void *ptr, int64_t subscript[], int64_t ld[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _subscript[MAXDIM];
     Integer _ld[MAXDIM];

     COPYINDEX_C2F(subscript, _subscript, ndim);
     COPYINDEX_C2F(ld, _ld, ndim-1);
     wnga_access_ghost_element_ptr(a, ptr, _subscript, _ld);
}

void NGA_Release(int g_a, int lo[], int hi[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
     COPYINDEX_C2F(lo,_ga_lo,ndim);
     COPYINDEX_C2F(hi,_ga_hi,ndim);

     wnga_release(a,_ga_lo, _ga_hi);
}

void NGA_Release64(int g_a, int64_t lo[], int64_t hi[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
     COPYINDEX_C2F(lo,_ga_lo,ndim);
     COPYINDEX_C2F(hi,_ga_hi,ndim);

     wnga_release(a,_ga_lo, _ga_hi);
}

void NGA_Release_block(int g_a, int idx)
{
     Integer a=(Integer)g_a;
     Integer iblock = (Integer)idx;

     wnga_release_block(a, iblock);
}

void NGA_Release_block_grid(int g_a, int index[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM];
     COPYINDEX_C2F(index,_ga_lo,ndim);

     wnga_release_block_grid(a, _ga_lo);
}

void NGA_Release_block_segment(int g_a, int idx)
{
     Integer a=(Integer)g_a;
     Integer iproc = (Integer)idx;

     wnga_release_block_segment(a, iproc);
}

void NGA_Release_ghost_element(int g_a, int index[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM];
     COPYINDEX_C2F(index,_ga_lo,ndim);

     wnga_release_ghost_element(a, _ga_lo);
}

void NGA_Release_ghost_element64(int g_a, int64_t index[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM];
     COPYINDEX_C2F(index,_ga_lo,ndim);

     wnga_release_ghost_element(a, _ga_lo);
}


void NGA_Release_ghosts(int g_a)
{
     Integer a=(Integer)g_a;

     wnga_release_ghosts(a);
}

int NGA_Locate(int g_a, int subscript[])
{
    logical st;
    Integer a=(Integer)g_a, owner;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM];
    COPYINDEX_C2F(subscript,_ga_lo,ndim);

    st = wnga_locate(a,_ga_lo,&owner);
    if(st == TRUE) return (int)owner;
    else return -1;
}

int NGA_Locate64(int g_a, int64_t subscript[])
{
    logical st;
    Integer a=(Integer)g_a, owner;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM];
    COPYINDEX_C2F(subscript,_ga_lo,ndim);

    st = wnga_locate(a,_ga_lo,&owner);
    if(st == TRUE) return (int)owner;
    else return -1;
}


int NGA_Locate_nnodes(int g_a, int lo[], int hi[])
{
     /* logical st; */
     Integer a=(Integer)g_a, np;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];

     COPYINDEX_C2F(lo,_ga_lo,ndim);
     COPYINDEX_C2F(hi,_ga_hi,ndim);
     /* st = wnga_locate_nnodes(a, _ga_lo, _ga_hi, &np); */
     (void)wnga_locate_nnodes(a, _ga_lo, _ga_hi, &np);
     return (int)np;
}


int NGA_Locate_nnodes64(int g_a, int64_t lo[], int64_t hi[])
{
     /* logical st; */
     Integer a=(Integer)g_a, np;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];

     COPYINDEX_C2F(lo,_ga_lo,ndim);
     COPYINDEX_C2F(hi,_ga_hi,ndim);
     /* st = wnga_locate_nnodes(a, _ga_lo, _ga_hi, &np); */
     (void)wnga_locate_nnodes(a, _ga_lo, _ga_hi, &np);
     return (int)np;
}


int NGA_Locate_region(int g_a,int lo[],int hi[],int map[],int procs[])
{
     logical st;
     Integer a=(Integer)g_a, np_guess, np_actual;
     Integer ndim = wnga_ndim(a);
     Integer *tmap;
     int i;
     Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
     Integer *_ga_map_capi;
     COPYINDEX_C2F(lo,_ga_lo,ndim);
     COPYINDEX_C2F(hi,_ga_hi,ndim);
     st = wnga_locate_nnodes(a, _ga_lo, _ga_hi, &np_guess);
     tmap = (Integer *)malloc( (int)(np_guess*2*ndim *sizeof(Integer)));
     if(!map)GA_Error("NGA_Locate_region: unable to allocate memory",g_a);
     _ga_map_capi = (Integer*)malloc(np_guess*sizeof(Integer));

     st = wnga_locate_region(a,_ga_lo, _ga_hi, tmap, _ga_map_capi, &np_actual);
     assert(np_guess == np_actual);
     if(st==FALSE){
       free(tmap);
       free(_ga_map_capi);
       return 0;
     }

     COPY(int,_ga_map_capi,procs, np_actual);

        /* might have to swap lo/hi when copying */

     for(i=0; i< np_actual*2; i++){
        Integer *ptmap = tmap+i*ndim;
        int *pmap = map +i*ndim;
        COPYINDEX_F2C(ptmap, pmap, ndim);  
     }
     free(tmap);
     free(_ga_map_capi);
     return (int)np_actual;
}

int NGA_Locate_region64(int g_a,int64_t lo[],int64_t hi[],int64_t map[],int procs[])
{
     logical st;
     Integer a=(Integer)g_a, np_guess, np_actual;
     Integer ndim = wnga_ndim(a);
     Integer *tmap;
     int i;
     Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
     Integer *_ga_map_capi;
     COPYINDEX_C2F(lo,_ga_lo,ndim);
     COPYINDEX_C2F(hi,_ga_hi,ndim);
     st = wnga_locate_nnodes(a, _ga_lo, _ga_hi, &np_guess);
     tmap = (Integer *)malloc( (int)(np_guess*2*ndim *sizeof(Integer)));
     if(!map)GA_Error("NGA_Locate_region: unable to allocate memory",g_a);
     _ga_map_capi = (Integer*)malloc(np_guess*sizeof(Integer));

     st = wnga_locate_region(a,_ga_lo, _ga_hi, tmap, _ga_map_capi, &np_actual);
     assert(np_guess == np_actual);
     if(st==FALSE){
       free(tmap);
       free(_ga_map_capi);
       return 0;
     }

     COPY(int,_ga_map_capi,procs, np_actual);

        /* might have to swap lo/hi when copying */

     for(i=0; i< np_actual*2; i++){
        Integer *ptmap = tmap+i*ndim;
        int64_t *pmap = map +i*ndim;
        COPYINDEX_F2C_64(ptmap, pmap, ndim);  
     }
     free(tmap);
     free(_ga_map_capi);
     return (int)np_actual;
}

int NGA_Locate_num_blocks(int g_a, int *lo, int *hi)
{
  Integer ret;
  Integer a = (Integer)g_a;
  Integer ndim = wnga_ndim(a);
  Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
  COPYINDEX_C2F(lo,_ga_lo,ndim);
  COPYINDEX_C2F(hi,_ga_hi,ndim);
  ret = wnga_locate_num_blocks(a, _ga_lo, _ga_hi);
  return (int)ret;
}


void NGA_Inquire(int g_a, int *type, int *ndim, int dims[])
{
     Integer a=(Integer)g_a;
     Integer ttype, nndim;
     Integer _ga_dims[MAXDIM];
     wnga_inquire(a,&ttype, &nndim, _ga_dims);
     COPYF2C(_ga_dims, dims,nndim);  
     *ndim = (int)nndim;
     *type = (int)ttype;
}

void NGA_Inquire64(int g_a, int *type, int *ndim, int64_t dims[])
{
     Integer a=(Integer)g_a;
     Integer ttype, nndim;
     Integer _ga_dims[MAXDIM];
     wnga_inquire(a,&ttype, &nndim, _ga_dims);
     COPYF2C_64(_ga_dims, dims,nndim);  
     *ndim = (int)nndim;
     *type = (int)ttype;
}

char* GA_Inquire_name(int g_a)
{
     Integer a=(Integer)g_a;
     char *ptr;
     wnga_inquire_name(a, &ptr);
     return(ptr);
}

char* NGA_Inquire_name(int g_a)
{
     Integer a=(Integer)g_a;
     char *ptr;
     wnga_inquire_name(a, &ptr);
     return(ptr);
}

size_t GA_Memory_avail(void)
{
    return (size_t)wnga_memory_avail();
}

size_t NGA_Memory_avail(void)
{
    return (size_t)wnga_memory_avail();
}

int GA_Memory_limited(void)
{
    if(wnga_memory_limited() ==TRUE) return 1;
    else return 0;
}

int NGA_Memory_limited(void)
{
    if(wnga_memory_limited() ==TRUE) return 1;
    else return 0;
}

void NGA_Proc_topology(int g_a, int proc, int coord[])
{
     Integer a=(Integer)g_a;
     Integer p=(Integer)proc;
     Integer ndim = wnga_ndim(a);
     Integer _ga_work[MAXDIM];
     wnga_proc_topology(a, p, _ga_work);
     COPY(int,_ga_work, coord,ndim);  
}

void GA_Get_proc_grid(int g_a, int dims[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_work[MAXDIM];
     wnga_get_proc_grid(a, _ga_work);
     COPY(int,_ga_work, dims ,ndim);  
}

void NGA_Get_proc_grid(int g_a, int dims[])
{
     Integer a=(Integer)g_a;
     Integer _ga_work[MAXDIM];
     Integer ndim = wnga_ndim(a);
     wnga_get_proc_grid(a, _ga_work);
     COPY(int,_ga_work, dims ,ndim);  
}

void GA_Check_handle(int g_a, char* string)
{
     Integer a=(Integer)g_a;
     wnga_check_handle(a,string);
}

int GA_Create_mutexes(int number)
{
     Integer n = (Integer)number;
     if(wnga_create_mutexes(n) == TRUE)return 1;
     else return 0;
}

int NGA_Create_mutexes(int number)
{
     Integer n = (Integer)number;
     if(wnga_create_mutexes(n) == TRUE)return 1;
     else return 0;
}

int GA_Destroy_mutexes(void) {
  if(wnga_destroy_mutexes() == TRUE) return 1;
  else return 0;
}

int NGA_Destroy_mutexes(void) {
  if(wnga_destroy_mutexes() == TRUE) return 1;
  else return 0;
}

void GA_Lock(int mutex)
{
     Integer m = (Integer)mutex;
     wnga_lock(m);
}

void NGA_Lock(int mutex)
{
     Integer m = (Integer)mutex;
     wnga_lock(m);
}

void GA_Unlock(int mutex)
{
     Integer m = (Integer)mutex;
     wnga_unlock(m);
}

void NGA_Unlock(int mutex)
{
     Integer m = (Integer)mutex;
     wnga_unlock(m);
}

void GA_Brdcst(void *buf, int lenbuf, int root)
{
  Integer type=GA_TYPE_BRD;
  Integer len = (Integer)lenbuf;
  Integer orig = (Integer)root;
  wnga_msg_brdcst(type, buf, len, orig);
}
   
void GA_Pgroup_brdcst(int grp_id, void *buf, int lenbuf, int root)
{
    Integer type=GA_TYPE_BRD;
    Integer len = (Integer)lenbuf;
    Integer orig = (Integer)root;
    Integer grp = (Integer)grp_id;
    wnga_pgroup_brdcst(grp, type, buf, len, orig);
}

void GA_Pgroup_sync(int grp_id)
{
    Integer grp = (Integer)grp_id;
    wnga_pgroup_sync(grp);
}

void NGA_Pgroup_sync(int grp_id)
{
    Integer grp = (Integer)grp_id;
    wnga_pgroup_sync(grp);
}

void GA_Gop(int type, void *x, int n, char *op)
{ wnga_gop(type, x, n, op); }

void GA_Igop(int x[], int n, char *op)
{ wnga_gop(C_INT, x, n, op); }

void GA_Lgop(long x[], int n, char *op)
{ wnga_gop(C_LONG, x, n, op); }

void GA_Llgop(long long x[], int n, char *op)
{ wnga_gop(C_LONGLONG, x, n, op); }

void GA_Fgop(float x[], int n, char *op)
{ wnga_gop(C_FLOAT, x, n, op); }       

void GA_Dgop(double x[], int n, char *op)
{ wnga_gop(C_DBL, x, n, op); }

void GA_Cgop(SingleComplex x[], int n, char *op)
{ wnga_gop(C_SCPL, x, n, op); }

void GA_Zgop(DoubleComplex x[], int n, char *op)
{ wnga_gop(C_DCPL, x, n, op); }

void NGA_Gop(int type, void *x, int n, char *op)
{ wnga_gop(type, x, n, op); }

void NGA_Igop(int x[], int n, char *op)
{ wnga_gop(C_INT, x, n, op); }

void NGA_Lgop(long x[], int n, char *op)
{ wnga_gop(C_LONG, x, n, op); }

void NGA_Llgop(long long x[], int n, char *op)
{ wnga_gop(C_LONGLONG, x, n, op); }

void NGA_Fgop(float x[], int n, char *op)
{ wnga_gop(C_FLOAT, x, n, op); }       

void NGA_Dgop(double x[], int n, char *op)
{ wnga_gop(C_DBL, x, n, op); }

void NGA_Cgop(SingleComplex x[], int n, char *op)
{ wnga_gop(C_SCPL, x, n, op); }

void NGA_Zgop(DoubleComplex x[], int n, char *op)
{ wnga_gop(C_DCPL, x, n, op); }

void GA_Pgroup_gop(int grp_id, int type, double x[], int n, char *op)
{ wnga_pgroup_gop(grp_id, type, x, n, op); }

void GA_Pgroup_igop(int grp_id, int x[], int n, char *op)
{ wnga_pgroup_gop(grp_id, C_INT, x, n, op); }
 
void GA_Pgroup_lgop(int grp_id, long x[], int n, char *op)
{ wnga_pgroup_gop(grp_id, C_LONG, x, n, op); }
 
void GA_Pgroup_llgop(int grp_id, long long x[], int n, char *op)
{ wnga_pgroup_gop(grp_id, C_LONGLONG, x, n, op); }
 
void GA_Pgroup_fgop(int grp_id, float x[], int n, char *op)
{ wnga_pgroup_gop(grp_id, C_FLOAT, x, n, op); }
 
void GA_Pgroup_dgop(int grp_id, double x[], int n, char *op)
{ wnga_pgroup_gop(grp_id, C_DBL, x, n, op); }
 
void GA_Pgroup_cgop(int grp_id, SingleComplex x[], int n, char *op)
{ wnga_pgroup_gop(grp_id, C_SCPL, x, n, op); }
 
void GA_Pgroup_zgop(int grp_id, DoubleComplex x[], int n, char *op)
{ wnga_pgroup_gop(grp_id, C_DCPL, x, n, op); }
 
void NGA_Pgroup_gop(int grp_id, int type, double x[], int n, char *op)
{ wnga_pgroup_gop(grp_id, type, x, n, op); }

void NGA_Pgroup_igop(int grp_id, int x[], int n, char *op)
{ wnga_pgroup_gop(grp_id, C_INT, x, n, op); }
 
void NGA_Pgroup_lgop(int grp_id, long x[], int n, char *op)
{ wnga_pgroup_gop(grp_id, C_LONG, x, n, op); }
 
void NGA_Pgroup_llgop(int grp_id, long long x[], int n, char *op)
{ wnga_pgroup_gop(grp_id, C_LONGLONG, x, n, op); }
 
void NGA_Pgroup_fgop(int grp_id, float x[], int n, char *op)
{ wnga_pgroup_gop(grp_id, C_FLOAT, x, n, op); }
 
void NGA_Pgroup_dgop(int grp_id, double x[], int n, char *op)
{ wnga_pgroup_gop(grp_id, C_DBL, x, n, op); }
 
void NGA_Pgroup_cgop(int grp_id, SingleComplex x[], int n, char *op)
{ wnga_pgroup_gop(grp_id, C_SCPL, x, n, op); }
 
void NGA_Pgroup_zgop(int grp_id, DoubleComplex x[], int n, char *op)
{ wnga_pgroup_gop(grp_id, C_DCPL, x, n, op); }
 
void NGA_Alloc_gatscat_buf(int nelems)
{
  Integer elems = (Integer)nelems;
  wnga_alloc_gatscat_buf(elems);
}

void NGA_Free_gatscat_buf()
{
  wnga_free_gatscat_buf();
}

void NGA_Scatter(int g_a, void *v, int* subsArray[], int n)
{
    int idx, i;
    Integer a = (Integer)g_a;
    Integer nv = (Integer)n;
#ifndef USE_GATSCAT_NEW
    Integer ndim = wnga_ndim(a);
    Integer *_subs_array;
    _subs_array = (Integer *)malloc((int)ndim* n * sizeof(Integer));
    if(_subs_array == NULL) GA_Error("Memory allocation failed.", 0);
    for(idx=0; idx<n; idx++)
        for(i=0; i<ndim; i++)
            _subs_array[idx*ndim+(ndim-i-1)] = subsArray[idx][i] + 1;
    wnga_scatter(a, v, _subs_array, 0, nv);
    free(_subs_array);
#else
    wnga_scatter(a, v, subsArray, 1, nv);
#endif
}

void NGA_Scatter_flat(int g_a, void *v, int subsArray[], int n)
{
    int idx, i;
    Integer a = (Integer)g_a;
    Integer nv = (Integer)n;
    Integer ndim = wnga_ndim(a);
    Integer *_subs_array;
    _subs_array = (Integer *)malloc((int)ndim* n * sizeof(Integer));
    if(_subs_array == NULL) GA_Error("Memory allocation failed.", 0);

    /* adjust the indices for fortran interface */
    for(idx=0; idx<n; idx++)
        for(i=0; i<ndim; i++)
            _subs_array[idx*ndim+(ndim-i-1)] = subsArray[idx*ndim+i] + 1;
    
    wnga_scatter(a, v, _subs_array, 0, nv);
    
    free(_subs_array);
}


void NGA_Scatter64(int g_a, void *v, int64_t* subsArray[], int64_t n)
{
    int64_t idx;
    int i;
    Integer a = (Integer)g_a;
    Integer nv = (Integer)n;
#ifndef USE_GATSCAT_NEW
    Integer ndim = wnga_ndim(a);
    Integer *_subs_array;
    _subs_array = (Integer *)malloc((int)ndim* n * sizeof(Integer));
    if(_subs_array == NULL) GA_Error("Memory allocation failed.", 0);
    for(idx=0; idx<n; idx++)
        for(i=0; i<ndim; i++)
            _subs_array[idx*ndim+(ndim-i-1)] = subsArray[idx][i] + 1;
    wnga_scatter(a, v, _subs_array, 0, nv);
    free(_subs_array);
#else
    wnga_scatter(a, v, subsArray, 1, nv);
#endif
}

void NGA_Scatter_flat64(int g_a, void *v, int64_t subsArray[], int64_t n)
{
    int idx, i;
    Integer a = (Integer)g_a;
    Integer nv = (Integer)n;
    Integer ndim = wnga_ndim(a);
    Integer *_subs_array;
    _subs_array = (Integer *)malloc((int)ndim* n * sizeof(Integer));
    if(_subs_array == NULL) GA_Error("Memory allocation failed.", 0);

    /* adjust the indices for fortran interface */
    for(idx=0; idx<n; idx++)
        for(i=0; i<ndim; i++)
            _subs_array[idx*ndim+(ndim-i-1)] = subsArray[idx*ndim+i] + 1;
    
    wnga_scatter(a, v, _subs_array, 0, nv);
    
    free(_subs_array);
}


void NGA_Scatter_acc(int g_a, void *v, int* subsArray[], int n, void *alpha)
{
    int idx, i;
    Integer a = (Integer)g_a;
    Integer nv = (Integer)n;
#ifndef USE_GATSCAT_NEW
    Integer ndim = wnga_ndim(a);
    Integer *_subs_array;
    _subs_array = (Integer *)malloc((int)ndim* n * sizeof(Integer));
    if(_subs_array == NULL) GA_Error("Memory allocation failed.", 0);
    for(idx=0; idx<n; idx++)
        for(i=0; i<ndim; i++)
            _subs_array[idx*ndim+(ndim-i-1)] = subsArray[idx][i] + 1;
    wnga_scatter_acc(a, v, _subs_array, 0, nv, alpha);
    free(_subs_array);
#else
    wnga_scatter_acc(a, v, subsArray, 1, nv, alpha);
#endif
}

void NGA_Scatter_acc_flat(int g_a, void *v, int subsArray[], int n, void *alpha)
{
    int idx, i;
    Integer a = (Integer)g_a;
    Integer nv = (Integer)n;
    Integer ndim = wnga_ndim(a);
    Integer *_subs_array;
    _subs_array = (Integer *)malloc((int)ndim* n * sizeof(Integer));
    if(_subs_array == NULL) GA_Error("Memory allocation failed.", 0);

    /* adjust the indices for fortran interface */
    for(idx=0; idx<n; idx++)
        for(i=0; i<ndim; i++)
            _subs_array[idx*ndim+(ndim-i-1)] = subsArray[idx*ndim+i] + 1;
    
    wnga_scatter_acc(a, v, _subs_array, 0, nv, alpha);
    
    free(_subs_array);
}

void NGA_Scatter_acc64(int g_a, void *v, int64_t* subsArray[], int64_t n, void *alpha)
{
    int64_t idx;
    int i;
    Integer a = (Integer)g_a;
    Integer nv = (Integer)n;
#ifndef USE_GATSCAT_NEW
    Integer ndim = wnga_ndim(a);
    Integer *_subs_array;
    _subs_array = (Integer *)malloc((int)ndim* n * sizeof(Integer));
    if(_subs_array == NULL) GA_Error("Memory allocation failed.", 0);
    for(idx=0; idx<n; idx++)
        for(i=0; i<ndim; i++)
            _subs_array[idx*ndim+(ndim-i-1)] = subsArray[idx][i] + 1;
    wnga_scatter_acc(a, v, _subs_array, 0, nv, alpha);
    free(_subs_array);
#else
    wnga_scatter_acc(a, v, subsArray, 1, nv, alpha);
#endif
}

void NGA_Scatter_acc_flat64(int g_a, void *v, int64_t subsArray[], int64_t n, void *alpha)
{
    int idx, i;
    Integer a = (Integer)g_a;
    Integer nv = (Integer)n;
    Integer ndim = wnga_ndim(a);
    Integer *_subs_array;
    _subs_array = (Integer *)malloc((int)ndim* n * sizeof(Integer));
    if(_subs_array == NULL) GA_Error("Memory allocation failed.", 0);

    /* adjust the indices for fortran interface */
    for(idx=0; idx<n; idx++)
        for(i=0; i<ndim; i++)
            _subs_array[idx*ndim+(ndim-i-1)] = subsArray[idx*ndim+i] + 1;
    
    wnga_scatter_acc(a, v, _subs_array, 0, nv, alpha);
    
    free(_subs_array);
}

void NGA_Gather(int g_a, void *v, int* subsArray[], int n)
{
    int idx, i;
    Integer a = (Integer)g_a;
    Integer nv = (Integer)n;
#ifndef USE_GATSCAT_NEW
    Integer ndim = wnga_ndim(a);
    Integer *_subs_array;
    _subs_array = (Integer *)malloc((int)ndim* n * sizeof(Integer));
    if(_subs_array == NULL) GA_Error("Memory allocation failed.", 0);

    /* adjust the indices for fortran interface */
    for(idx=0; idx<n; idx++)
        for(i=0; i<ndim; i++)
            _subs_array[idx*ndim+(ndim-i-1)] = subsArray[idx][i] + 1;
    wnga_gather(a, v, _subs_array, 0, nv);
    free(_subs_array);
#else
    wnga_gather(a, v, subsArray, 1, nv);
#endif
}


void NGA_Gather_flat(int g_a, void *v, int subsArray[], int n)
{
    int idx, i;
    Integer a = (Integer)g_a;
    Integer nv = (Integer)n;
    Integer ndim = wnga_ndim(a);
    Integer *_subs_array;
    _subs_array = (Integer *)malloc((int)ndim* n * sizeof(Integer));
    if(_subs_array == NULL) GA_Error("Memory allocation failed.", 0);

    /* adjust the indices for fortran interface */
    for(idx=0; idx<n; idx++)
        for(i=0; i<ndim; i++)
            _subs_array[idx*ndim+(ndim-i-1)] = subsArray[idx*ndim+i] + 1;
    
    wnga_gather(a, v, _subs_array, 0, nv);
    
    free(_subs_array);
}


void NGA_Gather64(int g_a, void *v, int64_t* subsArray[], int64_t n)
{
    int64_t idx;
    int i;
    Integer a = (Integer)g_a;
    Integer nv = (Integer)n;
#ifndef USE_GATSCAT_NEW
    Integer ndim = wnga_ndim(a);
    Integer *_subs_array;
    _subs_array = (Integer *)malloc((int)ndim* n * sizeof(Integer));
    if(_subs_array == NULL) GA_Error("Memory allocation failed.", 0);

    /* adjust the indices for fortran interface */
    for(idx=0; idx<n; idx++)
        for(i=0; i<ndim; i++)
            _subs_array[idx*ndim+(ndim-i-1)] = subsArray[idx][i] + 1;
    wnga_gather(a, v, _subs_array, 0, nv);
    free(_subs_array);
#else
    wnga_gather(a, v, subsArray, 1, nv);
#endif
}


void NGA_Gather_flat64(int g_a, void *v, int64_t subsArray[], int64_t n)
{
    int idx, i;
    Integer a = (Integer)g_a;
    Integer nv = (Integer)n;
    Integer ndim = wnga_ndim(a);
    Integer *_subs_array;
    _subs_array = (Integer *)malloc((int)ndim* n * sizeof(Integer));
    if(_subs_array == NULL) GA_Error("Memory allocation failed.", 0);

    /* adjust the indices for fortran interface */
    for(idx=0; idx<n; idx++)
        for(i=0; i<ndim; i++)
            _subs_array[idx*ndim+(ndim-i-1)] = subsArray[idx*ndim+i] + 1;
    
    wnga_gather(a, v, _subs_array, 0, nv);
    
    free(_subs_array);
}


void GA_Dgemm_c(char ta, char tb, int m, int n, int k,
              double alpha, int g_a, int g_b, double beta, int g_c )
{
  Integer G_a = g_a;
  Integer G_b = g_b;
  Integer G_c = g_c;

  Integer ailo = 1;
  Integer aihi = m;
  Integer ajlo = 1;
  Integer ajhi = k;
  
  Integer bilo = 1;
  Integer bihi = k;
  Integer bjlo = 1;
  Integer bjhi = n;
  
  Integer cilo = 1;
  Integer cihi = m;
  Integer cjlo = 1;
  Integer cjhi = n;
  
  pnga_matmul(&ta, &tb, (DoublePrecision *)&alpha,(DoublePrecision *)&beta,
	    G_a, ailo, aihi, ajlo, ajhi,
	    G_b, bilo, bihi, bjlo, bjhi,
	    G_c, cilo, cihi, cjlo, cjhi);
}

void GA_Zgemm_c(char ta, char tb, int m, int n, int k,
              DoubleComplex alpha, int g_a, int g_b, 
	      DoubleComplex beta, int g_c )
{
  Integer G_a = g_a;
  Integer G_b = g_b;
  Integer G_c = g_c;

  Integer ailo = 1;
  Integer aihi = m;
  Integer ajlo = 1;
  Integer ajhi = k;
  
  Integer bilo = 1;
  Integer bihi = k;
  Integer bjlo = 1;
  Integer bjhi = n;
  
  Integer cilo = 1;
  Integer cihi = m;
  Integer cjlo = 1;
  Integer cjhi = n;
  
  pnga_matmul(&ta, &tb, (DoublePrecision *)&alpha,(DoublePrecision *)&beta,
	    G_a, ailo, aihi, ajlo, ajhi,
	    G_b, bilo, bihi, bjlo, bjhi,
	    G_c, cilo, cihi, cjlo, cjhi);
}

void GA_Cgemm_c(char ta, char tb, int m, int n, int k,
              SingleComplex alpha, int g_a, int g_b, 
	      SingleComplex beta, int g_c )
{
  Integer G_a = g_a;
  Integer G_b = g_b;
  Integer G_c = g_c;

  Integer ailo = 1;
  Integer aihi = m;
  Integer ajlo = 1;
  Integer ajhi = k;
  
  Integer bilo = 1;
  Integer bihi = k;
  Integer bjlo = 1;
  Integer bjhi = n;
  
  Integer cilo = 1;
  Integer cihi = m;
  Integer cjlo = 1;
  Integer cjhi = n;
  
  pnga_matmul(&ta, &tb, (float *)&alpha,(float *)&beta,
	    G_a, ailo, aihi, ajlo, ajhi,
	    G_b, bilo, bihi, bjlo, bjhi,
	    G_c, cilo, cihi, cjlo, cjhi);
}

void GA_Sgemm_c(char ta, char tb, int m, int n, int k,
              float alpha, int g_a, int g_b, 
	      float beta,  int g_c )
{
  Integer G_a = g_a;
  Integer G_b = g_b;
  Integer G_c = g_c;

  Integer ailo = 1;
  Integer aihi = m;
  Integer ajlo = 1;
  Integer ajhi = k;
  
  Integer bilo = 1;
  Integer bihi = k;
  Integer bjlo = 1;
  Integer bjhi = n;
  
  Integer cilo = 1;
  Integer cihi = m;
  Integer cjlo = 1;
  Integer cjhi = n;
  
  pnga_matmul(&ta, &tb, (float*)&alpha, (float*)&beta,
	    G_a, ailo, aihi, ajlo, ajhi,
	    G_b, bilo, bihi, bjlo, bjhi,
	    G_c, cilo, cihi, cjlo, cjhi);
}

void GA_Dgemm64_c(char ta, char tb, int64_t m, int64_t n, int64_t k,
                double alpha, int g_a, int g_b, double beta, int g_c )
{
  Integer G_a = g_a;
  Integer G_b = g_b;
  Integer G_c = g_c;

  Integer ailo = 1;
  Integer aihi = m;
  Integer ajlo = 1;
  Integer ajhi = k;
  
  Integer bilo = 1;
  Integer bihi = k;
  Integer bjlo = 1;
  Integer bjhi = n;
  
  Integer cilo = 1;
  Integer cihi = m;
  Integer cjlo = 1;
  Integer cjhi = n;
  
  pnga_matmul(&ta, &tb, (DoublePrecision *)&alpha,(DoublePrecision *)&beta,
	    G_a, ailo, aihi, ajlo, ajhi,
	    G_b, bilo, bihi, bjlo, bjhi,
	    G_c, cilo, cihi, cjlo, cjhi);
}

void GA_Zgemm64_c(char ta, char tb, int64_t m, int64_t n, int64_t k,
                DoubleComplex alpha, int g_a, int g_b, 
                DoubleComplex beta, int g_c )
{
  Integer G_a = g_a;
  Integer G_b = g_b;
  Integer G_c = g_c;

  Integer ailo = 1;
  Integer aihi = m;
  Integer ajlo = 1;
  Integer ajhi = k;
  
  Integer bilo = 1;
  Integer bihi = k;
  Integer bjlo = 1;
  Integer bjhi = n;
  
  Integer cilo = 1;
  Integer cihi = m;
  Integer cjlo = 1;
  Integer cjhi = n;
  
  pnga_matmul(&ta, &tb, (DoublePrecision *)&alpha,(DoublePrecision *)&beta,
	    G_a, ailo, aihi, ajlo, ajhi,
	    G_b, bilo, bihi, bjlo, bjhi,
	    G_c, cilo, cihi, cjlo, cjhi);
}

void GA_Cgemm64_c(char ta, char tb, int64_t m, int64_t n, int64_t k,
                SingleComplex alpha, int g_a, int g_b, 
                SingleComplex beta, int g_c )
{
  Integer G_a = g_a;
  Integer G_b = g_b;
  Integer G_c = g_c;

  Integer ailo = 1;
  Integer aihi = m;
  Integer ajlo = 1;
  Integer ajhi = k;
  
  Integer bilo = 1;
  Integer bihi = k;
  Integer bjlo = 1;
  Integer bjhi = n;
  
  Integer cilo = 1;
  Integer cihi = m;
  Integer cjlo = 1;
  Integer cjhi = n;
  
  pnga_matmul(&ta, &tb, (float *)&alpha,(float *)&beta,
	    G_a, ailo, aihi, ajlo, ajhi,
	    G_b, bilo, bihi, bjlo, bjhi,
	    G_c, cilo, cihi, cjlo, cjhi);
}

void GA_Sgemm64_c(char ta, char tb, int64_t m, int64_t n, int64_t k,
                float alpha, int g_a, int g_b, 
                float beta,  int g_c )
{
  Integer G_a = g_a;
  Integer G_b = g_b;
  Integer G_c = g_c;

  Integer ailo = 1;
  Integer aihi = m;
  Integer ajlo = 1;
  Integer ajhi = k;
  
  Integer bilo = 1;
  Integer bihi = k;
  Integer bjlo = 1;
  Integer bjhi = n;
  
  Integer cilo = 1;
  Integer cihi = m;
  Integer cjlo = 1;
  Integer cjhi = n;
  
  pnga_matmul(&ta, &tb, (float*)&alpha, (float*)&beta,
	    G_a, ailo, aihi, ajlo, ajhi,
	    G_b, bilo, bihi, bjlo, bjhi,
	    G_c, cilo, cihi, cjlo, cjhi);
}

/**
 * When calling GA _dgemm from the C, it is represented as follows (since the
 * underlying GA dgemm implementation pnga_matmul is in Fortran style)
 *   C(m,n) = A(m,k) * B(k,n)
 * Since GA internally creates GAs in Fortran style, we remap the above
 * expression to: C(n,m) = B(n,k) * A(k,m) and pass it to fortran. This
 * produces the output C in (n,m) format, which translates to C as (m,n),and
 * ultimately giving us the correct result.
 */
void GA_Dgemm(char ta, char tb, int m, int n, int k,
              double alpha, int g_a, int g_b, double beta, int g_c )
{
   GA_Dgemm_c(tb, ta, n, m, k, alpha, g_b, g_a, beta, g_c);
}

void GA_Zgemm(char ta, char tb, int m, int n, int k,
              DoubleComplex alpha, int g_a, int g_b, 
	      DoubleComplex beta, int g_c )
{
    GA_Zgemm_c(tb, ta, n, m, k, alpha, g_b, g_a, beta, g_c);
}

void GA_Cgemm(char ta, char tb, int m, int n, int k,
              SingleComplex alpha, int g_a, int g_b, 
	      SingleComplex beta, int g_c )
{
    GA_Cgemm_c(tb, ta, n, m, k, alpha, g_b, g_a, beta, g_c);
}

void GA_Sgemm(char ta, char tb, int m, int n, int k,
              float alpha, int g_a, int g_b, 
	      float beta,  int g_c )
{
    GA_Sgemm_c(tb, ta, n, m, k, alpha, g_b, g_a, beta, g_c);
}

void GA_Dgemm64(char ta, char tb, int64_t m, int64_t n, int64_t k,
                double alpha, int g_a, int g_b, double beta, int g_c )
{
    GA_Dgemm64_c(tb, ta, n, m, k, alpha, g_b, g_a, beta, g_c);
}

void GA_Zgemm64(char ta, char tb, int64_t m, int64_t n, int64_t k,
                DoubleComplex alpha, int g_a, int g_b, 
                DoubleComplex beta, int g_c )
{
    GA_Zgemm64_c(tb, ta, n, m, k, alpha, g_b, g_a, beta, g_c);
}

void GA_Cgemm64(char ta, char tb, int64_t m, int64_t n, int64_t k,
                SingleComplex alpha, int g_a, int g_b, 
                SingleComplex beta, int g_c )
{
    GA_Cgemm64_c(tb, ta, n, m, k, alpha, g_b, g_a, beta, g_c);
}

void GA_Sgemm64(char ta, char tb, int64_t m, int64_t n, int64_t k,
                float alpha, int g_a, int g_b, 
                float beta,  int g_c )
{
    GA_Sgemm64_c(tb, ta, n, m, k, alpha, g_b, g_a, beta, g_c);
}

/* Patch related */

void NGA_Copy_patch(char trans, int g_a, int alo[], int ahi[],
                                int g_b, int blo[], int bhi[])
{
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);

    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
    
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);

    wnga_copy_patch(&trans, a, _ga_alo, _ga_ahi, b, _ga_blo, _ga_bhi);
}

void NGA_Copy_patch64(char trans, int g_a, int64_t alo[], int64_t ahi[],
                                  int g_b, int64_t blo[], int64_t bhi[])
{
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);

    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
    
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);

    wnga_copy_patch(&trans, a, _ga_alo, _ga_ahi, b, _ga_blo, _ga_bhi);
}

void GA_Matmul_patch(char transa, char transb, void* alpha, void *beta,
		     int g_a, int ailo, int aihi, int ajlo, int ajhi,
		     int g_b, int bilo, int bihi, int bjlo, int bjhi,
		     int g_c, int cilo, int cihi, int cjlo, int cjhi)

{
    Integer Ig_a = g_a;
    Integer Ig_b = g_b;
    Integer Ig_c = g_c;
#if 0
    Integer alo[2], ahi[2]; 
    Integer blo[2], bhi[2];
    Integer clo[2], chi[2];
    alo[0]=ailo+1; ahi[0]=aihi+1; alo[1]=ajlo+1; ahi[1]=ajhi+1;
    blo[0]=bilo+1; bhi[0]=bihi+1; blo[1]=bjlo+1; bhi[1]=bjhi+1;
    clo[0]=cilo+1; chi[0]=cihi+1; clo[1]=cjlo+1; chi[1]=cjhi+1;
    pnga_matmul_patch(transa, transb, alpha, beta, g_a, alo, ahi,
                         g_b, blo, bhi, g_c, clo, chi);
#else
    Integer Iailo=ailo+1, Iaihi=aihi+1, Iajlo=ajlo+1, Iajhi=ajhi+1;
    Integer Ibilo=bilo+1, Ibihi=bihi+1, Ibjlo=bjlo+1, Ibjhi=bjhi+1;
    Integer Icilo=cilo+1, Icihi=cihi+1, Icjlo=cjlo+1, Icjhi=cjhi+1;
    if(pnga_is_mirrored(Ig_a))
       wnga_matmul_mirrored(&transa, &transb, (void*)alpha, (void*)beta,
                  Ig_a, Iailo, Iaihi, Iajlo, Iajhi,
                  Ig_b, Ibilo, Ibihi, Ibjlo, Ibjhi,
                  Ig_c, Icilo, Icihi, Icjlo, Icjhi);
    else {
       gai_matmul_patch_flag(SET);
       wnga_matmul(&transa, &transb, (void*)alpha, (void*)beta,
             Ig_a, Iailo, Iaihi, Iajlo, Iajhi,
             Ig_b, Ibilo, Ibihi, Ibjlo, Ibjhi,
             Ig_c, Icilo, Icihi, Icjlo, Icjhi);
       gai_matmul_patch_flag(UNSET);
    }
#endif
}

void GA_Matmul_patch64(char transa, char transb, void* alpha, void *beta,
                       int g_a, int64_t ailo, int64_t aihi, int64_t ajlo, int64_t ajhi,
                       int g_b, int64_t bilo, int64_t bihi, int64_t bjlo, int64_t bjhi,
                       int g_c, int64_t cilo, int64_t cihi, int64_t cjlo, int64_t cjhi)

{
    Integer Ig_a = g_a;
    Integer Ig_b = g_b;
    Integer Ig_c = g_c;
#if 0
    Integer alo[2], ahi[2]; 
    Integer blo[2], bhi[2];
    Integer clo[2], chi[2];
    alo[0]=ailo+1; ahi[0]=aihi+1; alo[1]=ajlo+1; ahi[1]=ajhi+1;
    blo[0]=bilo+1; bhi[0]=bihi+1; blo[1]=bjlo+1; bhi[1]=bjhi+1;
    clo[0]=cilo+1; chi[0]=cihi+1; clo[1]=cjlo+1; chi[1]=cjhi+1;
    pnga_matmul_patch(transa, transb, alpha, beta, g_a, alo, ahi,
                         g_b, blo, bhi, g_c, clo, chi);
#else
    Integer Iailo=ailo+1, Iaihi=aihi+1, Iajlo=ajlo+1, Iajhi=ajhi+1;
    Integer Ibilo=bilo+1, Ibihi=bihi+1, Ibjlo=bjlo+1, Ibjhi=bjhi+1;
    Integer Icilo=cilo+1, Icihi=cihi+1, Icjlo=cjlo+1, Icjhi=cjhi+1;
    if(pnga_is_mirrored(Ig_a))
       wnga_matmul_mirrored(&transa, &transb, (void*)alpha, (void*)beta,
                  Ig_a, Iailo, Iaihi, Iajlo, Iajhi,
                  Ig_b, Ibilo, Ibihi, Ibjlo, Ibjhi,
                  Ig_c, Icilo, Icihi, Icjlo, Icjhi);
    else {
       gai_matmul_patch_flag(SET);
       wnga_matmul(&transa, &transb, (void*)alpha, (void*)beta,
             Ig_a, Iailo, Iaihi, Iajlo, Iajhi,
             Ig_b, Ibilo, Ibihi, Ibjlo, Ibjhi,
             Ig_c, Icilo, Icihi, Icjlo, Icjhi);
       gai_matmul_patch_flag(UNSET);
    }
#endif
}

void NGA_Matmul_patch(char transa, char transb, void* alpha, void *beta,
		      int g_a, int alo[], int ahi[], 
		      int g_b, int blo[], int bhi[], 
		      int g_c, int clo[], int chi[]) 

{
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);
    
    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
    
    Integer c=(Integer)g_c;
    Integer cndim = wnga_ndim(c);
    
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    Integer _ga_clo[MAXDIM], _ga_chi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
    
    COPYINDEX_C2F(clo,_ga_clo, cndim);
    COPYINDEX_C2F(chi,_ga_chi, cndim);
    
    pnga_matmul_patch(&transb, &transa, alpha, beta,
		     b, _ga_blo, _ga_bhi,
		     a, _ga_alo, _ga_ahi,
		     c, _ga_clo, _ga_chi);
}

void NGA_Matmul_patch64(char transa, char transb, void* alpha, void *beta,
                        int g_a, int64_t alo[], int64_t ahi[], 
                        int g_b, int64_t blo[], int64_t bhi[], 
                        int g_c, int64_t clo[], int64_t chi[]) 

{
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);
    
    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
    
    Integer c=(Integer)g_c;
    Integer cndim = wnga_ndim(c);
    
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    Integer _ga_clo[MAXDIM], _ga_chi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
    
    COPYINDEX_C2F(clo,_ga_clo, cndim);
    COPYINDEX_C2F(chi,_ga_chi, cndim);
    
    pnga_matmul_patch(&transb, &transa, alpha, beta,
		     b, _ga_blo, _ga_bhi,
		     a, _ga_alo, _ga_ahi,
		     c, _ga_clo, _ga_chi);
}

int NGA_Idot_patch(int g_a, char t_a, int alo[], int ahi[],
                   int g_b, char t_b, int blo[], int bhi[])
{
    int res;
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);

    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
    
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);

    wnga_dot_patch(a, &t_a, _ga_alo, _ga_ahi,
                   b, &t_b, _ga_blo, _ga_bhi, &res);

    return res;
}

long NGA_Ldot_patch(int g_a, char t_a, int alo[], int ahi[],
                    int g_b, char t_b, int blo[], int bhi[])
{
    long res;
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);

    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
    
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);

    wnga_dot_patch(a, &t_a, _ga_alo, _ga_ahi,
                   b, &t_b, _ga_blo, _ga_bhi, &res);

    return res;
}

long long NGA_Lldot_patch(int g_a, char t_a, int alo[], int ahi[],
                    int g_b, char t_b, int blo[], int bhi[])
{
    long res;
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);

    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
    
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);

    wnga_dot_patch(a, &t_a, _ga_alo, _ga_ahi,
                   b, &t_b, _ga_blo, _ga_bhi, &res);

    return res;
}

double NGA_Ddot_patch(int g_a, char t_a, int alo[], int ahi[],
                   int g_b, char t_b, int blo[], int bhi[])
{
    double res;
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);

    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);

    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);

    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);

    wnga_dot_patch(a, &t_a, _ga_alo, _ga_ahi,
                   b, &t_b, _ga_blo, _ga_bhi, &res);

    return res;
}

DoubleComplex NGA_Zdot_patch(int g_a, char t_a, int alo[], int ahi[],
                             int g_b, char t_b, int blo[], int bhi[])
{
    DoubleComplex res;
    
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);
    
    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
    
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
    
    wnga_dot_patch(a, &t_a, _ga_alo, _ga_ahi,
                   b, &t_b, _ga_blo, _ga_bhi, &res);
    
    return res;
}

SingleComplex NGA_Cdot_patch(int g_a, char t_a, int alo[], int ahi[],
                             int g_b, char t_b, int blo[], int bhi[])
{
    SingleComplex res;
    
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);
    
    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
    
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
    
    wnga_dot_patch(a, &t_a, _ga_alo, _ga_ahi,
                   b, &t_b, _ga_blo, _ga_bhi, &res);
    
    return res;
}

float NGA_Fdot_patch(int g_a, char t_a, int alo[], int ahi[],
                   int g_b, char t_b, int blo[], int bhi[])
{
    float res;
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);
 
    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
 
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
 
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
 
    wnga_dot_patch(a, &t_a, _ga_alo, _ga_ahi,
                   b, &t_b, _ga_blo, _ga_bhi, &res);
 
    return res;
}                                           

int NGA_Idot_patch64(int g_a, char t_a, int64_t alo[], int64_t ahi[],
                     int g_b, char t_b, int64_t blo[], int64_t bhi[])
{
    int res;
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);

    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
    
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);

    wnga_dot_patch(a, &t_a, _ga_alo, _ga_ahi,
                   b, &t_b, _ga_blo, _ga_bhi, &res);

    return res;
}

long NGA_Ldot_patch64(int g_a, char t_a, int64_t alo[], int64_t ahi[],
                     int g_b, char t_b, int64_t blo[], int64_t bhi[])
{
    long res;
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);

    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
    
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);

    wnga_dot_patch(a, &t_a, _ga_alo, _ga_ahi,
                   b, &t_b, _ga_blo, _ga_bhi, &res);

    return res;
}

long long NGA_Lldot_patch64(int g_a, char t_a, int64_t alo[], int64_t ahi[],
                     int g_b, char t_b, int64_t blo[], int64_t bhi[])
{
    long long res;
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);

    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
    
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);

    wnga_dot_patch(a, &t_a, _ga_alo, _ga_ahi,
                   b, &t_b, _ga_blo, _ga_bhi, &res);

    return res;
}

double NGA_Ddot_patch64(int g_a, char t_a, int64_t alo[], int64_t ahi[],
                        int g_b, char t_b, int64_t blo[], int64_t bhi[])
{
    double res;
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);

    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);

    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);

    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);

    wnga_dot_patch(a, &t_a, _ga_alo, _ga_ahi,
                   b, &t_b, _ga_blo, _ga_bhi, &res);

    return res;
}

DoubleComplex NGA_Zdot_patch64(int g_a, char t_a, int64_t alo[], int64_t ahi[],
                               int g_b, char t_b, int64_t blo[], int64_t bhi[])
{
    DoubleComplex res;
    
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);
    
    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
    
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
    
    wnga_dot_patch(a, &t_a, _ga_alo, _ga_ahi,
                   b, &t_b, _ga_blo, _ga_bhi, &res);
    
    return res;
}

SingleComplex NGA_Cdot_patch64(int g_a, char t_a, int64_t alo[], int64_t ahi[],
                               int g_b, char t_b, int64_t blo[], int64_t bhi[])
{
    SingleComplex res;
    
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);
    
    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
    
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
    
    wnga_dot_patch(a, &t_a, _ga_alo, _ga_ahi,
                   b, &t_b, _ga_blo, _ga_bhi, &res);
    
    return res;
}

float NGA_Fdot_patch64(int g_a, char t_a, int64_t alo[], int64_t ahi[],
                       int g_b, char t_b, int64_t blo[], int64_t bhi[])
{
    float res;
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);
 
    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);
 
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
 
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
 
    wnga_dot_patch(a, &t_a, _ga_alo, _ga_ahi,
                   b, &t_b, _ga_blo, _ga_bhi, &res);
 
    return res;
}


void NGA_Fill_patch(int g_a, int lo[], int hi[], void *val)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);

    wnga_fill_patch(a, _ga_lo, _ga_hi, val);
}

void NGA_Fill_patch64(int g_a, int64_t lo[], int64_t hi[], void *val)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);

    wnga_fill_patch(a, _ga_lo, _ga_hi, val);
}

void NGA_Zero_patch(int g_a, int lo[], int hi[])
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);

    wnga_zero_patch(a, _ga_lo, _ga_hi);
}

void NGA_Zero_patch64(int g_a, int64_t lo[], int64_t  hi[])
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);

    wnga_zero_patch(a, _ga_lo, _ga_hi);
}

void NGA_Scale_patch(int g_a, int lo[], int hi[], void *alpha)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);

    wnga_scale_patch(a, _ga_lo, _ga_hi, alpha);
}

void NGA_Scale_patch64(int g_a, int64_t lo[], int64_t hi[], void *alpha)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);

    wnga_scale_patch(a, _ga_lo, _ga_hi, alpha);
}

void NGA_Add_patch(void * alpha, int g_a, int alo[], int ahi[],
                   void * beta,  int g_b, int blo[], int bhi[],
                   int g_c, int clo[], int chi[])
{
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);

    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);

    Integer c=(Integer)g_c;
    Integer cndim = wnga_ndim(c);
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    Integer _ga_clo[MAXDIM], _ga_chi[MAXDIM];

    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);

    COPYINDEX_C2F(clo,_ga_clo, cndim);
    COPYINDEX_C2F(chi,_ga_chi, cndim);

    wnga_add_patch(alpha, a, _ga_alo, _ga_ahi, beta, b, _ga_blo, _ga_bhi,
                   c, _ga_clo, _ga_chi);
}

void NGA_Add_patch64(void * alpha, int g_a, int64_t alo[], int64_t ahi[],
                     void * beta,  int g_b, int64_t blo[], int64_t bhi[],
                                   int g_c, int64_t clo[], int64_t chi[])
{
    Integer a=(Integer)g_a;
    Integer andim = wnga_ndim(a);

    Integer b=(Integer)g_b;
    Integer bndim = wnga_ndim(b);

    Integer c=(Integer)g_c;
    Integer cndim = wnga_ndim(c);

    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    Integer _ga_clo[MAXDIM], _ga_chi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);

    COPYINDEX_C2F(clo,_ga_clo, cndim);
    COPYINDEX_C2F(chi,_ga_chi, cndim);

    wnga_add_patch(alpha, a, _ga_alo, _ga_ahi, beta, b, _ga_blo, _ga_bhi,
                   c, _ga_clo, _ga_chi);
}


void GA_Print_patch(int g_a,int ilo,int ihi,int jlo,int jhi,int pretty)
{
    Integer a = (Integer)g_a;
    Integer lo[2];
    Integer hi[2];
    Integer p = (Integer) pretty;
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];

    lo[0] = ilo; lo[1] = jlo;
    hi[0] = ihi; lo[1] = jhi;
    COPYINDEX_C2F(lo,_ga_lo,2);
    COPYINDEX_C2F(hi,_ga_hi,2);
    wnga_print_patch(a, _ga_lo, _ga_hi, p);
}


void NGA_Print_patch(int g_a, int lo[], int hi[], int pretty)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer p = (Integer)pretty;
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];

    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    wnga_print_patch(a, _ga_lo, _ga_hi, p);
}

void NGA_Print_patch64(int g_a, int64_t lo[], int64_t hi[], int pretty)
{
    Integer a=(Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer p = (Integer)pretty;
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
  
    wnga_print_patch(a, _ga_lo, _ga_hi, p);
}

void GA_Print(int g_a)
{
    Integer a=(Integer)g_a;
    wnga_print(a);
}

void GA_Print_file(FILE *file, int g_a)
{
  Integer G_a = g_a;
  wnga_print_file(file, G_a);
}

void GA_Diag_seq(int g_a, int g_s, int g_v, void *eval)
{
    Integer a = (Integer)g_a;
    Integer s = (Integer)g_s;
    Integer v = (Integer)g_v;

    wnga_diag_seq(a, s, v, eval);
}

void GA_Diag_std_seq(int g_a, int g_v, void *eval)
{
    Integer a = (Integer)g_a;
    Integer v = (Integer)g_v;

    wnga_diag_std_seq(a, v, eval);
}

void GA_Diag(int g_a, int g_s, int g_v, void *eval)
{
    Integer a = (Integer)g_a;
    Integer s = (Integer)g_s;
    Integer v = (Integer)g_v;

    wnga_diag(a, s, v, eval);
}

void GA_Diag_std(int g_a, int g_v, void *eval)
{
    Integer a = (Integer)g_a;
    Integer v = (Integer)g_v;

    wnga_diag_std(a, v, eval);
}

void GA_Diag_reuse(int reuse, int g_a, int g_s, int g_v, void *eval)
{
    Integer r = (Integer)reuse;
    Integer a = (Integer)g_a;
    Integer s = (Integer)g_s;
    Integer v = (Integer)g_v;

    wnga_diag_reuse(r, a, s, v, eval);
}

void GA_Lu_solve(char tran, int g_a, int g_b)
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;

    Integer t;

    if(tran == 't' || tran == 'T') t = 1;
    else t = 0;

    wnga_lu_solve_alt(t, a, b);
}

int GA_Llt_solve(int g_a, int g_b)
{
    Integer res;
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;

    res = wnga_llt_solve(a, b);

    return((int)res);
}

int GA_Solve(int g_a, int g_b)
{
    Integer res;
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;

    res = wnga_solve(a, b);

    return((int)res);
}

int GA_Spd_invert(int g_a)
{
    Integer res;
    Integer a = (Integer)g_a;

    res = wnga_spd_invert(a);

    return((int)res);
}

void GA_Summarize(int verbose)
{
    Integer v = (Integer)verbose;

    wnga_summarize(v);
}

void GA_Symmetrize(int g_a)
{
    Integer a = (Integer)g_a;

    wnga_symmetrize(a);
}

void NGA_Symmetrize(int g_a)
{
    Integer a = (Integer)g_a;

    wnga_symmetrize(a);
}

void GA_Transpose(int g_a, int g_b)
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;

    wnga_transpose(a, b);
}


void GA_Print_distribution(int g_a)
{
#ifdef USE_FAPI
    wnga_print_distribution(1,(Integer)g_a);
#else
    wnga_print_distribution(0,(Integer)g_a);
#endif
}


void NGA_Release_update(int g_a, int lo[], int hi[])
{
  Integer a = (Integer)g_a;
  Integer ndim = wnga_ndim(a);
  Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
  COPYINDEX_C2F(lo,_ga_lo,ndim);
  COPYINDEX_C2F(hi,_ga_hi,ndim);

  wnga_release_update(a,_ga_lo, _ga_hi);
}

void NGA_Release_update64(int g_a, int64_t lo[], int64_t hi[])
{
  Integer a = (Integer)g_a;
  Integer ndim = wnga_ndim(a);
  Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
  COPYINDEX_C2F(lo,_ga_lo,ndim);
  COPYINDEX_C2F(hi,_ga_hi,ndim);

  wnga_release_update(a,_ga_lo, _ga_hi);
}

void NGA_Release_update_block(int g_a, int idx)
{
     Integer a=(Integer)g_a;
     Integer iblock = (Integer)idx;

     wnga_release_update_block(a, iblock);
}

void NGA_Release_update_block_grid(int g_a, int index[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM];
     COPYINDEX_C2F(index,_ga_lo,ndim);
     wnga_release_update_block_grid(a, _ga_lo);
}

void NGA_Release_update_block_segment(int g_a, int idx)
{
     Integer a=(Integer)g_a;
     Integer iproc = (Integer)idx;

     wnga_release_update_block_segment(a, iproc);
}

void NGA_Release_update_ghost_element(int g_a, int index[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM];
     COPYINDEX_C2F(index,_ga_lo,ndim);
     wnga_release_update_ghost_element(a, _ga_lo);
}

void NGA_Release_update_ghost_element64(int g_a, int64_t index[])
{
     Integer a=(Integer)g_a;
     Integer ndim = wnga_ndim(a);
     Integer _ga_lo[MAXDIM];
     COPYINDEX_C2F(index,_ga_lo,ndim);
     wnga_release_update_ghost_element(a, _ga_lo);
}

void NGA_Release_update_ghosts(int g_a)
{
     Integer a=(Integer)g_a;
     wnga_release_update_ghosts(a);
}

int GA_Ndim(int g_a)
{
    Integer a = (Integer)g_a;
    return (int)wnga_ndim(a);
}

int NGA_Ndim(int g_a)
{
    Integer a = (Integer)g_a;
    return (int)wnga_ndim(a);
}

/*Limin's functions */


void GA_Step_bound_info(int g_xx, int g_vv, int g_xxll, int g_xxuu,  void *boundmin, void *wolfemin, void *boundmax)
{
    Integer xx = (Integer)g_xx;
    Integer vv = (Integer)g_vv;
    Integer xxll = (Integer)g_xxll;
    Integer xxuu = (Integer)g_xxuu;
    wnga_step_bound_info(xx, vv, xxll, xxuu, boundmin, wolfemin, boundmax);
}

void GA_Step_bound_info_patch(int g_xx, int xxlo[], int xxhi[],  int g_vv, int vvlo[], int vvhi[], int g_xxll, int xxlllo[], int xxllhi[], int g_xxuu,  int xxuulo[], int xxuuhi[], void *boundmin, void *wolfemin, void *boundmax)
{
    Integer xx = (Integer)g_xx;
    Integer vv = (Integer)g_vv;
    Integer xxll = (Integer)g_xxll;
    Integer xxuu = (Integer)g_xxuu;
    Integer ndim = wnga_ndim(xx);
    Integer _ga_xxlo[MAXDIM], _ga_xxhi[MAXDIM];
    Integer _ga_vvlo[MAXDIM], _ga_vvhi[MAXDIM];
    Integer _ga_xxlllo[MAXDIM], _ga_xxllhi[MAXDIM];
    Integer _ga_xxuulo[MAXDIM], _ga_xxuuhi[MAXDIM];
    COPYINDEX_C2F(xxlo,_ga_xxlo, ndim);
    COPYINDEX_C2F(xxhi,_ga_xxhi, ndim);
    COPYINDEX_C2F(vvlo,_ga_vvlo, ndim);
    COPYINDEX_C2F(vvhi,_ga_vvhi, ndim);
    COPYINDEX_C2F(xxlllo,_ga_xxlllo, ndim);
    COPYINDEX_C2F(xxllhi,_ga_xxllhi, ndim);
    COPYINDEX_C2F(xxuulo,_ga_xxuulo, ndim);
    COPYINDEX_C2F(xxuuhi,_ga_xxuuhi, ndim);
    wnga_step_bound_info_patch(xx, _ga_xxlo, _ga_xxhi, vv, _ga_vvlo, _ga_vvhi, xxll, _ga_xxlllo, _ga_xxllhi, xxuu, _ga_xxuulo, _ga_xxuuhi , 
    boundmin,wolfemin,boundmax);
}

void GA_Step_bound_info_patch64(int g_xx, int64_t xxlo[], int64_t xxhi[],
                                int g_vv, int64_t vvlo[], int64_t vvhi[],
                                int g_xxll, int64_t xxlllo[], int64_t xxllhi[],
                                int64_t g_xxuu,  int64_t xxuulo[], int64_t xxuuhi[],
                                void *boundmin, void *wolfemin, void *boundmax)
{
    Integer xx = (Integer)g_xx;
    Integer vv = (Integer)g_vv;
    Integer xxll = (Integer)g_xxll;
    Integer xxuu = (Integer)g_xxuu;
    Integer ndim = wnga_ndim(xx);
    Integer _ga_xxlo[MAXDIM], _ga_xxhi[MAXDIM];
    Integer _ga_vvlo[MAXDIM], _ga_vvhi[MAXDIM];
    Integer _ga_xxlllo[MAXDIM], _ga_xxllhi[MAXDIM];
    Integer _ga_xxuulo[MAXDIM], _ga_xxuuhi[MAXDIM];
    COPYINDEX_C2F(xxlo,_ga_xxlo, ndim);
    COPYINDEX_C2F(xxhi,_ga_xxhi, ndim);
    COPYINDEX_C2F(vvlo,_ga_vvlo, ndim);
    COPYINDEX_C2F(vvhi,_ga_vvhi, ndim);
    COPYINDEX_C2F(xxlllo,_ga_xxlllo, ndim);
    COPYINDEX_C2F(xxllhi,_ga_xxllhi, ndim);
    COPYINDEX_C2F(xxuulo,_ga_xxuulo, ndim);
    COPYINDEX_C2F(xxuuhi,_ga_xxuuhi, ndim);
    wnga_step_bound_info_patch(xx, _ga_xxlo, _ga_xxhi, vv, _ga_vvlo, _ga_vvhi, xxll, _ga_xxlllo, _ga_xxllhi, xxuu, _ga_xxuulo, _ga_xxuuhi , 
    boundmin,wolfemin,boundmax);
}

void GA_Step_max(int g_a, int g_b, void *step)
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    wnga_step_max(a, b, step);
}

void GA_Step_max_patch(int g_a, int alo[], int ahi[], int g_b, int blo[], int bhi[], void *step)
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    Integer ndim = wnga_ndim(a);
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, ndim);
    COPYINDEX_C2F(ahi,_ga_ahi, ndim);
    COPYINDEX_C2F(blo,_ga_blo, ndim);
    COPYINDEX_C2F(bhi,_ga_bhi, ndim);
    wnga_step_max_patch(a, _ga_alo, _ga_ahi, b, _ga_blo, _ga_bhi, step);
}

void GA_Step_max_patch64(int g_a, int64_t alo[], int64_t  ahi[],
                         int g_b, int64_t blo[], int64_t  bhi[], void *step)
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    Integer ndim = wnga_ndim(a);
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, ndim);
    COPYINDEX_C2F(ahi,_ga_ahi, ndim);
    COPYINDEX_C2F(blo,_ga_blo, ndim);
    COPYINDEX_C2F(bhi,_ga_bhi, ndim);
    wnga_step_max_patch(a, _ga_alo, _ga_ahi, b, _ga_blo, _ga_bhi, step);
}

void GA_Abs_value(int g_a)
{
    Integer a = (Integer)g_a;
    wnga_abs_value(a);
}

void GA_Add_constant(int g_a, void *alpha)
{
    Integer a = (Integer)g_a;
    wnga_add_constant(a, alpha);
}


void GA_Recip(int g_a)
{
    Integer a = (Integer)g_a;
    wnga_recip(a);
}

void GA_Elem_multiply(int g_a, int g_b, int g_c)
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    Integer c = (Integer)g_c;
    wnga_elem_multiply(a, b, c);
}

void GA_Elem_divide(int g_a, int g_b, int g_c)
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    Integer c = (Integer)g_c;
    wnga_elem_divide(a, b, c);
}

void GA_Elem_maximum(int g_a, int g_b, int g_c)
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    Integer c = (Integer)g_c;
    wnga_elem_maximum(a, b, c);
}


void GA_Elem_minimum(int g_a, int g_b, int g_c)
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    Integer c = (Integer)g_c;
    wnga_elem_minimum(a, b, c);
}


void GA_Abs_value_patch(int g_a, int *lo, int *hi)
{
    Integer a = (Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    wnga_abs_value_patch(a,_ga_lo, _ga_hi);
}

void GA_Abs_value_patch64(int g_a, int64_t *lo, int64_t *hi)
{
    Integer a = (Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    wnga_abs_value_patch(a,_ga_lo, _ga_hi);
}

void GA_Add_constant_patch(int g_a, int *lo, int* hi, void *alpha)
{
    Integer a = (Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    wnga_add_constant_patch(a, _ga_lo, _ga_hi, alpha);
}

void GA_Add_constant_patch64(int g_a, int64_t *lo, int64_t *hi, void *alpha)
{
    Integer a = (Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    wnga_add_constant_patch(a, _ga_lo, _ga_hi, alpha);
}

void GA_Recip_patch(int g_a, int *lo, int *hi)
{
    Integer a = (Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    wnga_recip_patch(a,_ga_lo, _ga_hi);
}

void GA_Recip_patch64(int g_a, int64_t *lo,  int64_t *hi)
{
    Integer a = (Integer)g_a;
    Integer ndim = wnga_ndim(a);
    Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM];
    COPYINDEX_C2F(lo,_ga_lo, ndim);
    COPYINDEX_C2F(hi,_ga_hi, ndim);
    wnga_recip_patch(a,_ga_lo, _ga_hi);
}

void GA_Elem_multiply_patch(int g_a, int alo[], int ahi[],
                            int g_b, int blo[], int bhi[],
                            int g_c, int clo[], int chi[])
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    Integer c = (Integer)g_c;
    Integer andim = wnga_ndim(a);
    Integer bndim = wnga_ndim(b);
    Integer cndim = wnga_ndim(c);
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    Integer _ga_clo[MAXDIM], _ga_chi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
    COPYINDEX_C2F(clo,_ga_clo, cndim);
    COPYINDEX_C2F(chi,_ga_chi, cndim);
    wnga_elem_multiply_patch(a, _ga_alo, _ga_ahi, b, _ga_blo, _ga_bhi, c, _ga_clo, _ga_chi);
}

void GA_Elem_multiply_patch64(int g_a, int64_t alo[], int64_t ahi[],
                              int g_b, int64_t blo[], int64_t bhi[],
                              int g_c, int64_t clo[], int64_t chi[])
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    Integer c = (Integer)g_c;
    Integer andim = wnga_ndim(a);
    Integer bndim = wnga_ndim(b);
    Integer cndim = wnga_ndim(c);
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    Integer _ga_clo[MAXDIM], _ga_chi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
    COPYINDEX_C2F(clo,_ga_clo, cndim);
    COPYINDEX_C2F(chi,_ga_chi, cndim);
    wnga_elem_multiply_patch(a, _ga_alo, _ga_ahi, b, _ga_blo, _ga_bhi, c, _ga_clo, _ga_chi);
}

void GA_Elem_divide_patch(int g_a, int alo[], int ahi[],
                          int g_b, int blo[], int bhi[],
                          int g_c, int clo[], int chi[])
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    Integer c = (Integer)g_c;
    Integer andim = wnga_ndim(a);
    Integer bndim = wnga_ndim(b);
    Integer cndim = wnga_ndim(c);
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    Integer _ga_clo[MAXDIM], _ga_chi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
    COPYINDEX_C2F(clo,_ga_clo, cndim);
    COPYINDEX_C2F(chi,_ga_chi, cndim);
    wnga_elem_divide_patch(a, _ga_alo, _ga_ahi, b, _ga_blo, _ga_bhi, c, _ga_clo, _ga_chi);
}

void GA_Elem_divide_patch64(int g_a, int64_t alo[], int64_t ahi[],
                            int g_b, int64_t blo[], int64_t bhi[],
                            int g_c, int64_t clo[], int64_t chi[])
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    Integer c = (Integer)g_c;
    Integer andim = wnga_ndim(a);
    Integer bndim = wnga_ndim(b);
    Integer cndim = wnga_ndim(c);
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    Integer _ga_clo[MAXDIM], _ga_chi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
    COPYINDEX_C2F(clo,_ga_clo, cndim);
    COPYINDEX_C2F(chi,_ga_chi, cndim);
    wnga_elem_divide_patch(a, _ga_alo, _ga_ahi, b, _ga_blo, _ga_bhi, c, _ga_clo, _ga_chi);
}


void GA_Elem_maximum_patch(int g_a, int alo[], int ahi[],
                           int g_b, int blo[], int bhi[],
                           int g_c, int clo[], int chi[])
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    Integer c = (Integer)g_c;
    Integer andim = wnga_ndim(a);
    Integer bndim = wnga_ndim(b);
    Integer cndim = wnga_ndim(c);
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    Integer _ga_clo[MAXDIM], _ga_chi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
    COPYINDEX_C2F(clo,_ga_clo, cndim);
    COPYINDEX_C2F(chi,_ga_chi, cndim);
    wnga_elem_maximum_patch(a, _ga_alo, _ga_ahi, b, _ga_blo, _ga_bhi, c, _ga_clo, _ga_chi);
}

void GA_Elem_maximum_patch64(int g_a, int64_t alo[], int64_t ahi[],
                             int g_b, int64_t blo[], int64_t bhi[],
                             int g_c, int64_t clo[], int64_t chi[])
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    Integer c = (Integer)g_c;
    Integer andim = wnga_ndim(a);
    Integer bndim = wnga_ndim(b);
    Integer cndim = wnga_ndim(c);
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    Integer _ga_clo[MAXDIM], _ga_chi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
    COPYINDEX_C2F(clo,_ga_clo, cndim);
    COPYINDEX_C2F(chi,_ga_chi, cndim);
    wnga_elem_maximum_patch(a, _ga_alo, _ga_ahi, b, _ga_blo, _ga_bhi, c, _ga_clo, _ga_chi);
}

void GA_Elem_minimum_patch(int g_a, int alo[], int ahi[],
                           int g_b, int blo[], int bhi[],
                           int g_c, int clo[], int chi[])
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    Integer c = (Integer)g_c;
    Integer andim = wnga_ndim(a);
    Integer bndim = wnga_ndim(b);
    Integer cndim = wnga_ndim(c);
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    Integer _ga_clo[MAXDIM], _ga_chi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
    COPYINDEX_C2F(clo,_ga_clo, cndim);
    COPYINDEX_C2F(chi,_ga_chi, cndim);
    wnga_elem_minimum_patch(a, _ga_alo, _ga_ahi, b, _ga_blo, _ga_bhi, c, _ga_clo, _ga_chi);
}

void GA_Elem_minimum_patch64(int g_a, int64_t alo[], int64_t ahi[],
                             int g_b, int64_t blo[], int64_t bhi[],
                             int g_c, int64_t clo[], int64_t chi[])
{
    Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    Integer c = (Integer)g_c;
    Integer andim = wnga_ndim(a);
    Integer bndim = wnga_ndim(b);
    Integer cndim = wnga_ndim(c);
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    Integer _ga_clo[MAXDIM], _ga_chi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
    COPYINDEX_C2F(clo,_ga_clo, cndim);
    COPYINDEX_C2F(chi,_ga_chi, cndim);
    wnga_elem_minimum_patch(a, _ga_alo, _ga_ahi, b, _ga_blo, _ga_bhi, c, _ga_clo, _ga_chi);
}

void GA_Shift_diagonal(int g_a, void *c){
 Integer a = (Integer )g_a;
 wnga_shift_diagonal(a, c);
}

void GA_Set_diagonal(int g_a, int g_v){
 Integer a = (Integer )g_a;
 Integer v = (Integer )g_v;
 wnga_set_diagonal(a, v);
}

void GA_Zero_diagonal(int g_a){
 Integer a = (Integer )g_a;
 wnga_zero_diagonal(a);
}
void GA_Add_diagonal(int g_a, int g_v){
 Integer a = (Integer )g_a;
 Integer v = (Integer )g_v;
 wnga_add_diagonal(a, v);
}

void GA_Get_diag(int g_a, int g_v){
 Integer a = (Integer )g_a;
 Integer v = (Integer )g_v;
 wnga_get_diag(a, v);
}

void GA_Scale_rows(int g_a, int g_v){
 Integer a = (Integer )g_a;
 Integer v = (Integer )g_v;
 wnga_scale_rows(a, v);
}

void GA_Scale_cols(int g_a, int g_v){
 Integer a = (Integer )g_a;
 Integer v = (Integer )g_v;
 wnga_scale_cols(a, v);
}
void GA_Norm1(int g_a, double *nm){
 Integer a = (Integer )g_a;
 wnga_norm1(a, nm);
}

void GA_Norm_infinity(int g_a, double *nm){
 Integer a = (Integer )g_a;
 wnga_norm_infinity(a, nm);
}


void GA_Median(int g_a, int g_b, int g_c, int g_m){
 Integer a = (Integer )g_a;
 Integer b = (Integer )g_b;
 Integer c= (Integer )g_c;
 Integer m = (Integer )g_m;
 wnga_median(a, b, c, m);
}

void GA_Median_patch(int g_a, int *alo, int *ahi,
                     int g_b, int *blo, int *bhi,
                     int g_c, int *clo, int *chi,
                     int g_m, int *mlo, int *mhi){

   Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    Integer c = (Integer)g_c;
    Integer m = (Integer )g_m;
    Integer andim = wnga_ndim(a);
    Integer bndim = wnga_ndim(b);
    Integer cndim = wnga_ndim(c);
    Integer mndim = wnga_ndim(m);
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    Integer _ga_clo[MAXDIM], _ga_chi[MAXDIM];
    Integer _ga_mlo[MAXDIM], _ga_mhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
    COPYINDEX_C2F(clo,_ga_clo, cndim);
    COPYINDEX_C2F(chi,_ga_chi, cndim);
    COPYINDEX_C2F(mlo,_ga_mlo, mndim);
    COPYINDEX_C2F(mhi,_ga_mhi, mndim);
    wnga_median_patch(a, _ga_alo, _ga_ahi, b, _ga_blo, _ga_bhi, c, _ga_clo, _ga_chi, m, _ga_mlo, _ga_mhi);
}

void GA_Median_patch64(int g_a, int64_t *alo, int64_t *ahi,
                       int g_b, int64_t *blo, int64_t *bhi,
                       int g_c, int64_t *clo, int64_t *chi,
                       int g_m, int64_t *mlo, int64_t *mhi){

   Integer a = (Integer)g_a;
    Integer b = (Integer)g_b;
    Integer c = (Integer)g_c;
    Integer m = (Integer )g_m;
    Integer andim = wnga_ndim(a);
    Integer bndim = wnga_ndim(b);
    Integer cndim = wnga_ndim(c);
    Integer mndim = wnga_ndim(m);
    Integer _ga_alo[MAXDIM], _ga_ahi[MAXDIM];
    Integer _ga_blo[MAXDIM], _ga_bhi[MAXDIM];
    Integer _ga_clo[MAXDIM], _ga_chi[MAXDIM];
    Integer _ga_mlo[MAXDIM], _ga_mhi[MAXDIM];
    COPYINDEX_C2F(alo,_ga_alo, andim);
    COPYINDEX_C2F(ahi,_ga_ahi, andim);
    COPYINDEX_C2F(blo,_ga_blo, bndim);
    COPYINDEX_C2F(bhi,_ga_bhi, bndim);
    COPYINDEX_C2F(clo,_ga_clo, cndim);
    COPYINDEX_C2F(chi,_ga_chi, cndim);
    COPYINDEX_C2F(mlo,_ga_mlo, mndim);
    COPYINDEX_C2F(mhi,_ga_mhi, mndim);
    wnga_median_patch(a, _ga_alo, _ga_ahi, b, _ga_blo, _ga_bhi, c, _ga_clo, _ga_chi, m, _ga_mlo, _ga_mhi);
}

/* return number of nodes being used in a cluster */
int GA_Cluster_nnodes()
{
    return wnga_cluster_nnodes();
} 

/* returns ClusterNode id of the calling process */
int GA_Cluster_nodeid() 
{
    return wnga_cluster_nodeid();
}

/* returns ClusterNode id of the specified process */
int GA_Cluster_proc_nodeid(int proc)
{
    Integer aproc = proc;
    return wnga_cluster_proc_nodeid(aproc);
}

/* return number of processes being used on the specified node */
int GA_Cluster_nprocs(int x) 
{
    Integer ax = x;
    return wnga_cluster_nprocs(ax);
}

/* global id of the calling process */
int GA_Cluster_procid(int node, int loc_proc)
{
    Integer anode = node;
    Integer aloc_proc = loc_proc;
    return wnga_cluster_procid(anode, aloc_proc);
}

/* wrapper for timer routines */
double GA_Wtime()
{
    return (double)wnga_wtime();
}

double NGA_Wtime()
{
    return (double)wnga_wtime();
}

void GA_Set_debug(int flag)
{
    Integer aa;
    aa = (Integer)flag;
    wnga_set_debug(aa);
}

void NGA_Set_debug(int flag)
{
    Integer aa;
    aa = (Integer)flag;
    wnga_set_debug(aa);
}

int GA_Get_debug()
{
    return (int)wnga_get_debug();
}

int NGA_Get_debug()
{
    return (int)wnga_get_debug();
}

#ifdef ENABLE_CHECKPOINT
void GA_Checkpoint(int* gas, int num)
{
    wnga_checkpoint_arrays(gas,num);
}
#endif

int GA_Pgroup_absolute_id(int grp_id, int pid) {
  Integer agrp_id = (Integer)grp_id;
  Integer apid = (Integer) pid;
  return (int)wnga_pgroup_absolute_id(agrp_id, apid);
}

int NGA_Pgroup_absolute_id(int grp_id, int pid) {
  Integer agrp_id = (Integer)grp_id;
  Integer apid = (Integer) pid;
  return (int)wnga_pgroup_absolute_id(agrp_id, apid);
}

void GA_Error(char *str, int code)
{
    Integer icode = code;
    wnga_error(str, icode);
}

void NGA_Error(char *str, int code)
{
    Integer icode = code;
    wnga_error(str, icode);
}

size_t GA_Inquire_memory()
{
    return (size_t)wnga_inquire_memory();
}

size_t NGA_Inquire_memory()
{
    return (size_t)wnga_inquire_memory();
}

void GA_Sync()
{
    wnga_sync();
}

void NGA_Sync()
{
    wnga_sync();
}

int GA_Uses_ma()
{
    return wnga_uses_ma();
}

int NGA_Uses_ma()
{
    return wnga_uses_ma();
}

void GA_Print_stats()
{
    wnga_print_stats();
}

void GA_Init_fence()
{
    wnga_init_fence();
}

void NGA_Init_fence()
{
    wnga_init_fence();
}

void GA_Fence()
{
    wnga_fence();
}

void NGA_Fence()
{
    wnga_fence();
}

int GA_Nodeid()
{
    return wnga_nodeid();
}

int NGA_Nodeid()
{
    return wnga_nodeid();
}

int GA_Nnodes()
{
    return wnga_nnodes();
}

int NGA_Nnodes()
{
    return wnga_nnodes();
}

static Integer* copy_map(int block[], int block_ndim, int map[])
{
    int d;
    int i,sum=0,capi_offset=0,map_offset=0;
    Integer *_ga_map_capi;

    for (d=0; d<block_ndim; d++) {
        sum += block[d];
    }

    _ga_map_capi = (Integer*)malloc(sum * sizeof(Integer));

    capi_offset = sum;
    for (d=0; d<block_ndim; d++) {
        capi_offset -= block[d];
        for (i=0; i<block[d]; i++) {
            _ga_map_capi[capi_offset+i] = map[map_offset+i] + 1;
        }
        map_offset += block[d];
    }

    return _ga_map_capi;
}

static Integer* copy_map64(int64_t block[], int block_ndim, int64_t map[])
{
    int d;
    int64_t i,sum=0,capi_offset=0,map_offset=0;
    Integer *_ga_map_capi;

    for (d=0; d<block_ndim; d++) {
        sum += block[d];
    }

    _ga_map_capi = (Integer*)malloc(sum * sizeof(Integer));

    capi_offset = sum;
    for (d=0; d<block_ndim; d++) {
        capi_offset -= block[d];
        for (i=0; i<block[d]; i++) {
            _ga_map_capi[capi_offset+i] = map[map_offset+i] + 1;
        }
        map_offset += block[d];
    }

    return _ga_map_capi;
}

int NGA_Register_type(size_t bytes) {
  return wnga_register_type(bytes);
}

int NGA_Deregister_type(int type) {
  return wnga_deregister_type(type);
}


void NGA_Get_field(int g_a, int *lo, int *hi, int foff, int fsize,
		   void *buf, int *ld) {
  Integer a = (Integer)g_a;
  Integer andim = wnga_ndim(a);
  Integer _alo[MAXDIM], _ahi[MAXDIM];
  Integer _ld[MAXDIM];
  COPYINDEX_C2F(lo,_alo, andim);
  COPYINDEX_C2F(hi,_ahi, andim);
  COPYINDEX_C2F(ld, _ld, andim-1);

  wnga_get_field(a, _alo, _ahi, foff, fsize, buf, _ld);
}

void NGA_Nbget_field(int g_a, int *lo, int *hi, int foff, int fsize,
		     void *buf, int *ld, ga_nbhdl_t *nbhandle) {
  Integer a = (Integer)g_a;
  Integer andim = wnga_ndim(a);
  Integer _alo[MAXDIM], _ahi[MAXDIM];
  Integer _ld[MAXDIM];
  COPYINDEX_C2F(lo,_alo, andim);
  COPYINDEX_C2F(hi,_ahi, andim);
  COPYINDEX_C2F(ld, _ld, andim-1);

  wnga_nbget_field(a, _alo, _ahi, foff, fsize, buf, _ld, (Integer*)nbhandle);
}

void NGA_Nbput_field(int g_a, int *lo, int *hi, int foff, int fsize,
		     void *buf, int *ld, ga_nbhdl_t *nbhandle) {
  Integer a = (Integer)g_a;
  Integer andim = wnga_ndim(a);
  Integer _alo[MAXDIM], _ahi[MAXDIM];
  Integer _ld[MAXDIM];
  COPYINDEX_C2F(lo,_alo, andim);
  COPYINDEX_C2F(hi,_ahi, andim);
  COPYINDEX_C2F(ld, _ld, andim-1);

  wnga_nbput_field(a, _alo, _ahi, foff, fsize, buf, _ld, (Integer*)nbhandle);
}

void NGA_Put_field(int g_a, int *lo, int *hi, int foff, int fsize,
		   void *buf, int *ld) {
  Integer a = (Integer)g_a;
  Integer andim = wnga_ndim(a);
  Integer _alo[MAXDIM], _ahi[MAXDIM];
  Integer _ld[MAXDIM];
  COPYINDEX_C2F(lo,_alo, andim);
  COPYINDEX_C2F(hi,_ahi, andim);
  COPYINDEX_C2F(ld, _ld, andim-1);

  wnga_put_field(a, _alo, _ahi, foff, fsize, buf, _ld);
}



