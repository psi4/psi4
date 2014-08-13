#if HAVE_CONFIG_H
#   include "config.h"
#endif


/** GA Memory Allocation Routines: uses either MA or external allocator */

#include "globalp.h"
#include "ga-papi.h"
#include "ga-wapi.h"
#define GA_MAXMEM_AVAIL ( ( (long)1 << (8*sizeof(Integer)-2) ) -1)
#define CHECK           0
#define ALIGNMENT       sizeof(DoubleComplex)

static void * (*ga_ext_alloc)(size_t, int, char *);
static void (*ga_ext_free)(void *);
short int ga_usesMA = 1; 

void GA_Register_stack_memory(
        void * (*ext_alloc)(size_t, int, char *),
        void   (*ext_free)(void *))
{
    if(ext_alloc == NULL || ext_free == NULL)
      pnga_error("GA_Register_stack_memory():Invalid pointer(s) passed\n", 0);
    ga_ext_alloc = ext_alloc; ga_ext_free  = ext_free; ga_usesMA=0;
}

void* ga_malloc(Integer nelem, int type, char *name)
{
    void *ptr;  
    unsigned long addr;
    Integer handle, adjust=0, bytes, item_size=GAsizeofM(pnga_type_f2c(type));
    Integer extra;

#if NOFORT
    type = pnga_type_f2c(type);
#endif

    /* extra space for 1.ALIGNMENT and 2.storing handle */
    if(ALIGNMENT%item_size)
       pnga_error("ga_malloc: GA datatype cannot be aligned.Adjust ALIGNMENT",0);
    extra = 2*ALIGNMENT/item_size;
    nelem += extra;

    if(ga_usesMA) { /* Uses Memory Allocator (MA) */
       if(MA_push_stack(type,nelem,name,&handle))  MA_get_pointer(handle,&ptr);
       else pnga_error("ga_malloc: MA_push_stack failed",0);
       addr = (unsigned long)ptr;
    }
    else { /* else, using external memory allocator */
       bytes = nelem*item_size;
       addr  = (unsigned long)(*ga_ext_alloc)(
               (size_t)bytes, (int)item_size, name);
    }

    /* Address Alignment */
    adjust = (Integer) (addr%ALIGNMENT);
    if(adjust != 0) { adjust=ALIGNMENT-adjust; addr+=adjust; }
    ptr = (void *)addr; 
    if(!ga_usesMA) handle = adjust;

    if(ptr == NULL) pnga_error("ga_malloc failed", 0L);
    *((Integer*)ptr)=handle;/*store handle or adjustment-value in this buffer*/
    ptr = ((char*)ptr) + ALIGNMENT;

    return ptr;
}

void ga_free(void *ptr)
{
    Integer handle;
    ptr = ((char*)ptr)-ALIGNMENT;
    handle= *((Integer*)ptr); /* retreive handle */

    if(ga_usesMA) {
      if(!MA_pop_stack(handle)) pnga_error("ga_free: MA_pop_stack failed",0);}
    else /*make sure to free original(before address alignment) pointer*/
      (*ga_ext_free)((char *)ptr - handle);
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_memory_avail_type = pnga_memory_avail_type
#endif
Integer pnga_memory_avail_type(Integer datatype)
{
    if(ga_usesMA)  return MA_inquire_avail(datatype);
    else return GA_MAXMEM_AVAIL;
}
