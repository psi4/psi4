/******************* header file for Disk Arrays *****************/
#ifndef _DRA_H_
#define _DRA_H_

#include "chemio.h"
#include "typesf2c.h"

typedef long dra_size_t;
#define DRA_RW ELIO_RW
#define DRA_R  ELIO_R
#define DRA_W  ELIO_W
#define DRA_REQ_INVALID -333

#ifdef __cplusplus
extern "C" {
#endif

/* C-interface prototypes */

extern int NDRA_Create(       int type,
                              int ndim,
                              dra_size_t dims[],
                              char *name,
                              char* filename,
                              int mode,
                              dra_size_t reqdims[],
                              int *d_a);

extern int NDRA_Inquire(      int d_a,
                              int *type,
                              int *ndim,
                              dra_size_t dims[],
                              char *name,
                              char* filename);

extern int NDRA_Write(        int g_a,
                              int d_a,
                              int *request);

extern int NDRA_Read(         int g_a,
                              int d_a,
                              int *request);

extern int NDRA_Write_section(logical transp,
                              int g_a,
                              int glo[],
                              int ghi[],
                              int d_a,
                              dra_size_t dlo[],
                              dra_size_t dhi[],
                              int *request);

extern int NDRA_Read_section( logical transp,
                              int g_a,
                              int glo[],
                              int ghi[],
                              int d_a,
                              dra_size_t dlo[],
                              dra_size_t dhi[],
                              int *request);

extern int DRA_Init(          int max_arrays,
                              double max_array_size,
                              double total_disk_space,
                              double max_memory);

extern int DRA_Terminate();

extern int DRA_Open(          char* filename,
                              int mode,
                              int *d_a);

extern int DRA_Probe(         int request,
                              int *compl_status);

extern void DRA_Set_debug(    logical flag);

extern void DRA_Print_internals(    int d_a);

extern void DRA_Set_default_config(    int numfiles, int numioprocs);

extern int DRA_Wait(          int request);

extern int DRA_Delete(        int d_a);

extern int DRA_Close(         int d_a);

extern void DRA_Flick();

#ifdef __cplusplus
}       
#endif

#endif /* _DRA_H_ */
