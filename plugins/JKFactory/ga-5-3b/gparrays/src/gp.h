#include <stdlib.h>

extern void GP_Access_element(int g_p, int *subscript, void *ptr, int *size);
extern int GP_Allocate(int g_p);
extern void GP_Assign_local_element(int g_p, int *subscript, void *ptr, int size);
extern int   GP_Create_handle();
extern void  GP_Debug(int g_p);
extern int   GP_Destroy(int g_p);
extern void  GP_Distribution(int g_p, int proc, int *lo, int *hi);
extern void  GP_Free(void* ptr);
extern void* GP_Free_local_element(int g_p, int *subscript);
extern int   GP_Get_dimension(int g_p);
extern void GP_Gather_size(int g_p, int nv, int *subscript, int *size);
extern void GP_Gather(int g_p, int nv, int *subscript, void *buf, void **buf_ptr, int *buf_size, int *size, int setbuf);
extern void  GP_Get_size(int g_p, int *lo, int *hi, int *size);
extern void  GP_Get(int g_p, int *lo, int *hi, void *buf, void **buf_ptr, int *ld, void *buf_size, int *ld_sz, int *size, int setbuf);
extern void  GP_Initialize();
extern void* GP_Malloc(size_t size);
extern void GP_Memzero(int g_p);
extern void GP_Put(int g_p, int *lo, int *hi, void **buf_ptr, int *ld, void *buf_size, int *ld_sz, int *size, int checksize);
extern void GP_Release_element(int g_p, int *subscript);
extern void GP_Release_update_element(int g_p, int *subscript);
extern void GP_Scatter(int g_p, int nv, int *subscript, void **buf_ptr, int *buf_size, int *size, int checksize);
extern void GP_Set_chunk(int g_p, int *chunk);
extern void GP_Set_dimensions(int g_p, int ndim, int *dims);
extern void GP_Set_irreg_distr(int g_p, int *mapc, int *blocks);
extern void GP_Sync();
extern void GP_Terminate();
