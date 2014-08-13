#ifndef __GPBASE_H__
#define __GPBASE_H__

#include "gacommon.h"
#include "typesf2c.h"

#if SIZEOF_VOIDP == SIZEOF_INT
#   define GP_Int int
#elif SIZEOF_VOIDP == SIZEOF_LONG
#   define GP_Int long
#else
#   error sizeof(void*) is not sizeof(int) nor sizeof(long)
#endif

/* Set maximum number of Global Pointer Arrays */
#define GP_MAX_ARRAYS 1024

/* Set maximum dimension of Global Pointer ARRAYS */
#define GP_MAX_DIM 7

/* Define handle numbering offset for GP Arrays */ 
#define GP_OFFSET 1000

typedef struct{
  Integer g_size_array;       /* Handle to Global Array holding sizes    */
  Integer g_ptr_array;        /* Handle to Global Array holding pointers */
  Integer active;             /* Handle is currently active              */
  Integer ndim;               /* Dimension of GP                         */
  Integer dims[GP_MAX_DIM];   /* Axes dimensions of GP                   */
  Integer lo[GP_MAX_DIM];     /* Lower indices of local block            */
  Integer hi[GP_MAX_DIM];     /* Upper indices of local block            */
  Integer ld[GP_MAX_DIM-1];   /* Stride of local block                   */
} gp_array_t;

extern gp_array_t *GP;
#endif /* __GPBASE_H__ */
