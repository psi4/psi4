/** @file
 * Private header file containing symbolic constants and type declarations
 * for internal C routines.
 *
 * This file should only be included by internal C files.
 */
#ifndef _ma_h
#define _ma_h

#include "macdecls.h"

/**
 ** function types
 **/

extern Boolean MAi_inform_base(Integer datatype, Pointer address1, Pointer address2);
extern void MAi_summarize_allocated_blocks(int index_base);

#define f2c_alloc_get_                  F77_FUNC_(f2c_alloc_get,F2C_ALLOC_GET)
#define f2c_allocate_heap_              F77_FUNC_(f2c_allocate_heap,F2C_ALLOCATE_HEAP)
#define f2c_chop_stack_                 F77_FUNC_(f2c_chop_stack,F2C_CHOP_STACK)
#define f2c_free_heap_                  F77_FUNC_(f2c_free_heap,F2C_FREE_HEAP)
#define f2c_free_heap_piece_            F77_FUNC_(f2c_free_heap_piece,F2C_FREE_HEAP_PIECE)
#define f2c_get_index_                  F77_FUNC_(f2c_get_index,F2C_GET_INDEX)
#define f2c_get_next_memhandle_         F77_FUNC_(f2c_get_next_memhandle,F2C_GET_NEXT_MEMHANDLE)
#define f2c_get_numalign_               F77_FUNC_(f2c_get_numalign,F2C_GET_NUMALIGN)
#define f2c_inform_base_                F77_FUNC_(f2c_inform_base,F2C_INFORM_BASE)
#define f2c_inform_base_fcd_            F77_FUNC_(f2c_inform_base_fcd,F2C_INFORM_BASE_FCD)
#define f2c_init_                       F77_FUNC_(f2c_init,F2C_INIT)
#define f2c_init_memhandle_iterator_    F77_FUNC_(f2c_init_memhandle_iterator,F2C_INIT_MEMHANDLE_ITERATOR)
#define f2c_initialized_                F77_FUNC_(f2c_initialized,F2C_INITIALIZED)
#define f2c_inquire_avail_              F77_FUNC_(f2c_inquire_avail,F2C_INQUIRE_AVAIL)
#define f2c_inquire_heap_               F77_FUNC_(f2c_inquire_heap,F2C_INQUIRE_HEAP)
#define f2c_inquire_stack_              F77_FUNC_(f2c_inquire_stack,F2C_INQUIRE_STACK)
#define f2c_inquire_heap_check_stack_   F77_FUNC_(f2c_inquire_heap_check_stack,F2C_INQUIRE_HEAP_CHECK_STACK)
#define f2c_inquire_heap_no_partition_  F77_FUNC_(f2c_inquire_heap_no_partition,F2C_INQUIRE_HEAP_NO_PARTITION)
#define f2c_inquire_stack_check_heap_   F77_FUNC_(f2c_inquire_stack_check_heap,F2C_INQUIRE_STACK_CHECK_HEAP)
#define f2c_inquire_stack_no_partition_ F77_FUNC_(f2c_inquire_stack_no_partition,F2C_INQUIRE_STACK_NO_PARTITION)
#define f2c_pop_stack_                  F77_FUNC_(f2c_pop_stack,F2C_POP_STACK)
#define f2c_print_stats_                F77_FUNC_(f2c_print_stats,F2C_PRINT_STATS)
#define f2c_push_get_                   F77_FUNC_(f2c_push_get,F2C_PUSH_GET)
#define f2c_push_stack_                 F77_FUNC_(f2c_push_stack,F2C_PUSH_STACK)
#define f2c_set_auto_verify_            F77_FUNC_(f2c_set_auto_verify,F2C_SET_AUTO_VERIFY)
#define f2c_set_error_print_            F77_FUNC_(f2c_set_error_print,F2C_SET_ERROR_PRINT)
#define f2c_set_hard_fail_              F77_FUNC_(f2c_set_hard_fail,F2C_SET_HARD_FAIL)
#define f2c_set_numalign_               F77_FUNC_(f2c_set_numalign,F2C_SET_NUMALIGN)
#define f2c_sizeof_                     F77_FUNC_(f2c_sizeof,F2C_SIZEOF)
#define f2c_sizeof_overhead_            F77_FUNC_(f2c_sizeof_overhead,F2C_SIZEOF_OVERHEAD)
#define f2c_summarize_allocated_blocks_ F77_FUNC_(f2c_summarize_allocated_blocks,F2C_SUMMARIZE_ALLOCATED_BLOCKS)
#define f2c_trace_                      F77_FUNC_(f2c_trace,F2C_TRACE)
#define f2c_verify_allocator_stuff_     F77_FUNC_(f2c_verify_allocator_stuff,F2C_VERIFY_ALLOCATOR_STUFF)
#define ma_set_sizes_                   F77_FUNC_(ma_set_sizes,MA_SET_SIZES)

#endif /* _ma_h */
