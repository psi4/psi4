/** @file
 * Public header file for a portable dynamic memory allocator.
 *
 * This file may be included by internal and external C files.
 */
#ifndef _macdecls_h
#define _macdecls_h

#ifdef __cplusplus
extern "C" {
#endif

#include "macommon.h"
#include "matypes.h"

/**
 ** constants
 **/

/* datatypes */
#define MT_CHAR     MT_C_CHAR     /**< char */
#define MT_INT      MT_C_INT      /**< int */
#define MT_LONGINT  MT_C_LONGINT  /**< long int */
#define MT_LONGLONG MT_C_LONGLONG /**< long long */
#define MT_REAL     MT_C_FLOAT    /**< float */
#define MT_DBL      MT_C_DBL      /**< double */
#define MT_LDBL     MT_C_LDBL     /**< long double */
#define MT_SCPL     MT_C_SCPL     /**< single precision complex */
#define MT_DCPL     MT_C_DCPL     /**< double precision complex */
#define MT_LDCPL    MT_C_LDCPL    /**< long double precision complex */
#define MT_C_FIRST  MT_CHAR       /**< first type */
#define MT_C_LAST   MT_LDCPL      /**< last type */

/**
 ** function types
 **/

extern Boolean MA_alloc_get(
    Integer     datatype,       /**< of elements in this block */
    Integer     nelem,          /**< # of elements in this block */
    const char  *name,          /**< assigned to this block by client */
    Integer     *memhandle,     /**< RETURN: handle for this block */
    MA_AccessIndex *index       /**< RETURN: index for this block */   );
extern Boolean MA_allocate_heap(
    Integer     datatype,       /**< of elements in this block */
    Integer     nelem,          /**< # of elements in this block */
    const char  *name,          /**< assigned to this block by client */
    Integer     *memhandle      /**< RETURN: handle for this block */ );
extern Boolean MA_chop_stack(Integer memhandle);
extern Boolean MA_free_heap(Integer memhandle);
extern Boolean MA_free_heap_piece(
    Integer     memhandle,      /**< the block to deallocate a piece of */
    Integer     nelem           /**< # of elements to deallocate */);
extern Boolean MA_get_index(
    Integer     memhandle,      /**< block to get index for */
    MA_AccessIndex *index       /**< RETURN: base index */);
extern Pointer MA_get_mbase(Integer datatype);   /**< to get base address of */
extern Boolean MA_get_next_memhandle(
    Integer     *ithandle,      /**< handle for this iterator */
    Integer     *memhandle      /**< RETURN: handle for the next block */);
extern Boolean MA_get_numalign(Integer *value);
extern Boolean MA_get_pointer(
    Integer     memhandle,      /**< block to get pointer for */
    void       *pointer         /**< JN: void** = void*    */ );
extern Boolean MA_init(
    Integer     datatype,      /**< for computing storage requirement */
    Integer     nominal_stack, /**< # of datatype elements desired for stack */
    Integer     nominal_heap   /**< # of datatype elements desired for heap */);
extern Boolean MA_initialized();
extern Boolean MA_init_memhandle_iterator( Integer *ithandle);
extern Integer MA_inquire_avail(Integer datatype);
extern Integer MA_inquire_heap(Integer datatype);
extern Integer MA_inquire_heap_check_stack(Integer datatype);
extern Integer MA_inquire_heap_no_partition(Integer datatype);
extern Integer MA_inquire_stack(Integer datatype);
extern Integer MA_inquire_stack_check_heap(Integer datatype);
extern Integer MA_inquire_stack_no_partition(Integer datatype);
extern Boolean MA_pop_stack(Integer memhandle);
extern void MA_print_stats(Boolean printroutines);
extern Boolean MA_push_get(
    Integer     datatype,       /**< of elements in this block */
    Integer     nelem,          /**< # of elements in this block */
    const char  *name,          /**< assigned to this block by client */
    Integer     *memhandle,     /**< RETURN: handle for this block */
    MA_AccessIndex *index       /**< RETURN: index for this block */);
extern Boolean MA_push_stack(
    Integer     datatype,       /**< of elements in this block */
    Integer     nelem,          /**< # of elements in this block */
    const char  *name,          /**< assigned to this block by client */
    Integer     *memhandle      /**< RETURN: handle for this block */);
extern Boolean MA_set_auto_verify(Boolean  value /* to set flag to */);
extern Boolean MA_set_error_print(Boolean value /* to set flag to */);
extern Boolean MA_set_hard_fail( Boolean value /* to set flag to */);
extern Boolean MA_set_numalign(Integer  value);
extern Integer MA_sizeof(
    Integer     datatype1,      /**< of source elements */
    Integer     nelem1,         /**< # of source elements */
    Integer     datatype2       /**< of target elements */);
extern Integer MA_sizeof_overhead(Integer datatype);
extern void MA_summarize_allocated_blocks();
extern void MA_trace(Boolean value);
extern Boolean MA_verify_allocator_stuff();
extern void MA_set_error_callback(void(*func)());

extern void ma_set_error_callback();

/**
 ** variables
 **/

/* base arrays for the C datatypes */
extern char                 ma_cb_char[];     /**< MT_C_CHAR */
extern int                  ma_cb_int[];      /**< MT_C_INT */
extern long                 ma_cb_long[];     /**< MT_C_LONGINT */
extern long long            ma_cb_longlong[]; /**< MT_C_LONGLONG */
extern float                ma_cb_float[];    /**< MT_C_FLOAT */
extern double               ma_cb_dbl[];      /**< MT_C_DBL */
extern MA_LongDouble        ma_cb_ldbl[];     /**< MT_C_LDBL */
extern MA_SingleComplex     ma_cb_scpl[];     /**< MT_C_SCPL */
extern MA_DoubleComplex     ma_cb_dcpl[];     /**< MT_C_DCPL */
extern MA_LongDoubleComplex ma_cb_ldcpl[];    /**< MT_C_LDCPL */

#ifdef __cplusplus
}
#endif

#endif /* _macdecls_h */
