#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*
 * Portable dynamic memory allocator.
 */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_MALLOC_H
#   include <malloc.h>
#endif
#include "error.h"
#include "farg.h"
#include "ma.h"
#include "memcpy.h"
#include "scope.h"
#include "table.h"

#ifdef ENABLE_ARMCI_MEM_OPTION
extern void* ARMCI_Malloc_local(long bytes);
#endif

/*
 * Memory layout:
 *
 *    segment = heap_region stack_region
 *    region = block block block ...
 *    block = AD gap1 guard1 client_space guard2 gap2
 *
 * A segment of memory is obtained from the OS upon initialization.
 * The low end of the segment is managed as a heap; the heap region
 * grows from low addresses to high addresses.  The high end of the
 * segment is managed as a stack; the stack region grows from high
 * addresses to low addresses.
 *
 * Each region consists of a series of contiguous blocks, one per
 * allocation request, and possibly some unused space.  Blocks in
 * the heap region are either in use by the client (allocated and
 * not yet deallocated) or not in use by the client (allocated and
 * already deallocated).  A block on the rightmost end of the heap
 * region becomes part of the unused space upon deallocation.
 * Blocks in the stack region are always in use by the client,
 * because when a stack block is deallocated, it becomes part of
 * the unused space.
 *
 * A block consists of the client space, i.e., the range of memory
 * available for use by the application; guard words adjacent to
 * each end of the client space to help detect improper memory access
 * by the client; bookkeeping info (in an "allocation descriptor,"
 * AD); and two gaps, each zero or more bytes long, to satisfy
 * alignment constraints (specifically, to ensure that AD and
 * client_space are aligned properly).
 */

/**
 ** constants
 **/

/* return value for returns that should never execute */
#define DONTCARE (Integer)0

/* default total # of bytes */
#define DEFAULT_TOTAL_HEAP  524288 /* 2^19 */
#define DEFAULT_TOTAL_STACK 524288 /* 2^19 */

/* estimate of max # of outstanding allocation requests */
#define DEFAULT_REQUESTS_HEAP  1
#define DEFAULT_REQUESTS_STACK 1

/* bytes per address */
#define BPA    1

/* per-allocation storage overhead, excluding alignment gaps */
#define BLOCK_OVERHEAD_FIXED (sizeof(AD) + (2 * sizeof(Guard)))

/* block lengths are integral multiples of this */
/*
 * Note that for machines such as the KSR on which sizeof(pointer)
 * and sizeof(long) are different than sizeof(int), alignment issues
 * can be tricky.  For example, the fields of a struct (e.g.,
 * client_space of AD) can be improperly aligned if the struct is
 * dynamically placed (by MA) in such a way that the first field is
 * properly aligned but sizes of subsequent fields accumulate to cause
 * a later field to be misaligned.  By defining the unit of alignment
 * to be the biggest of the integer and pointer types, part of the
 * problem is solved, but the sum of sizes of preceding fields can
 * still potentially cause difficulty.
 */
#if defined(BGP) || defined(BGQ)
#define ALIGNMENT	32
#else
#define ALIGNMENT	sizeof(long)
#endif

/* min size of block split and placed on free list */
#define MINBLOCKSIZE mai_round((long)(ALIGNMENT + BLOCK_OVERHEAD_FIXED), \
        (ulongi)ALIGNMENT)

/* signatures for guard words */
#define GUARD1 (Guard)0xaaaaaaaa /* start signature */
#define GUARD2 (Guard)0x55555555 /* stop signature */

/**
 ** types
 **/

typedef unsigned int Guard;   /* for detection of memory trashing */
typedef unsigned long ulongi; /* for brevity */

/* allocation request for a block */
typedef struct _AR
{
    Integer    datatype; /* of elements */
    Integer    nelem;    /* # of elements */
} AR;

/* allocation descriptor for a block */
typedef struct _AD
{
    Integer     datatype;          /* of elements */
    Integer     nelem;             /* # of elements */
    char        name[MA_NAMESIZE]; /* given by client */
    Pointer     client_space;      /* start of client space */
    ulongi      nbytes;            /* total # of bytes */
    struct _AD *next;              /* AD in linked list */
    ulongi      checksum;          /* of AD */
} AD;

/* block location for mh2ad */
typedef enum
{
    BL_HeapOrStack,
    BL_Heap,
    BL_Stack,
    BL_StackTop
} BlockLocation;

/**
 ** function types
 **/

private Boolean ad_big_enough(AD *ad, Pointer ar);
private Boolean ad_eq(AD *ad, Pointer ad_target);
private Boolean ad_gt(AD *ad, Pointer ad_target);
private Boolean ad_le(AD *ad, Pointer ad_target);
private Boolean ad_lt(AD *ad, Pointer ad_target);
private void ad_print(AD *ad, char *block_type);
private void balloc_after(AR *ar, Pointer address, Pointer *client_space, ulongi *nbytes);
private void balloc_before(AR *ar, Pointer address, Pointer *client_space, ulongi *nbytes);
private void block_free_heap(AD *ad);
private AD *block_split(AD *ad, ulongi bytes_needed, Boolean insert_free);
private ulongi checksum(AD *ad);

#ifdef DEBUG
private void debug_ad_print(AD *ad);
#endif /* DEBUG */

private Boolean guard_check(AD *ad);
private void guard_set(AD *ad);
private void list_coalesce(AD *list);
private AD *list_delete(AD *ad, AD **list);
private int list_delete_many(AD **list, Boolean (*pred)(), Pointer closure, void (*action)());
private AD *list_delete_one(AD **list, Boolean (*pred)(), Pointer closure);
private void list_insert(AD *ad, AD **list);
private void list_insert_ordered(AD *ad, AD **list, Boolean (*pred)());
private Boolean list_member(AD *ad, AD *list);
private int list_print(AD *list, char *block_type, int index_base);
private void list_verify(AD *list, char *block_type, char *preamble, int *blocks, int *bad_blocks, int *bad_checksums, int *bad_lguards, int *bad_rguards);
private Integer ma_max_heap_frag_nelem(Integer datatype, Integer min_nelem);
private Integer ma_nelem(Pointer address, ulongi length, Integer datatype);
private void ma_preinitialize(char *caller);
private Boolean mh2ad(Integer memhandle, AD **adout, BlockLocation location, char *caller);
private void mh_free(AD *ad);
private long mai_round(long value, ulongi unit);
private void str_ncopy(char *to, char *from, int maxchars);

/* foreign routines */

extern Integer ma_set_sizes_();    /* from the MA FORTRAN interface */

/**
 ** variables
 **/

/* base addresses of the datatypes */
private Pointer ma_base[] =
{
    (Pointer)ma_cb_char,    /* MT_C_CHAR */
    (Pointer)ma_cb_int,     /* MT_C_INT */
    (Pointer)ma_cb_long,    /* MT_C_LONGINT */
    (Pointer)ma_cb_float,   /* MT_C_FLOAT */
    (Pointer)ma_cb_dbl,     /* MT_C_DBL */
    (Pointer)ma_cb_ldbl,    /* MT_C_LDBL */
    (Pointer)ma_cb_scpl,    /* MT_C_SCPL */
    (Pointer)ma_cb_dcpl,    /* MT_C_DCPL */
    (Pointer)ma_cb_ldcpl,   /* MT_C_LDCPL */
    0,                      /* MT_F_BYTE */
    0,                      /* MT_F_INT */
    0,                      /* MT_F_LOG */
    0,                      /* MT_F_REAL */
    0,                      /* MT_F_DBL */
    0,                      /* MT_F_SCPL */
    0,                      /* MT_F_DCPL */
    (Pointer)ma_cb_longlong /* MT_C_LONGLONG */
};

/* names of the datatypes */
private char *ma_datatype[] =
{
    "char",
    "int",
    "long int",
    "float",
    "double",
    "long double",
    "single precision complex",
    "double precision complex",
    "long double precision complex",
    "byte",
    "integer",
    "logical",
    "real",
    "double precision",
    "single precision complex",
    "double precision complex",
    "long long"
};

/* numbers of bytes in the datatypes */
private int ma_sizeof[] =
{
    sizeof(char),                 /* MT_C_CHAR */
    sizeof(int),                  /* MT_C_INT */
    sizeof(long int),             /* MT_C_LONGINT */
    sizeof(float),                /* MT_C_FLOAT */
    sizeof(double),               /* MT_C_DBL */
    sizeof(MA_LongDouble),        /* MT_C_LDBL */
    sizeof(MA_SingleComplex),     /* MT_C_SCPL */
    sizeof(MA_DoubleComplex),     /* MT_C_DCPL */
    sizeof(MA_LongDoubleComplex), /* MT_C_LDCPL */
    0,                            /* MT_F_BYTE */
    0,                            /* MT_F_INT */
    0,                            /* MT_F_LOG */
    0,                            /* MT_F_REAL */
    0,                            /* MT_F_DBL */
    0,                            /* MT_F_SCPL */
    0,                            /* MT_F_DCPL */
    sizeof(long long)             /* MT_C_LONGLONG */
};

/*
 * Initially, ma_hp points to the start of the segment, and ma_sp
 * points to the first address past the end of the segment.  The
 * start of the segment is always pointed to by ma_segment, and
 * the first address past the end of the segment is always pointed
 * to by ma_eos.  The (unenforced) boundary between the heap region
 * and the stack region, defined at initialization, is always pointed
 * to by ma_partition.
 *
 *    ................................................
 *    ^                      ^                        ^
 *    ma_segment, ma_hp      ma_partition             ma_eos, ma_sp
 *
 * Later, ma_hp points to the first address past the end of the
 * rightmost heap block, and ma_sp points to the leftmost stack block.
 *
 *    hhhhhhhhhhhhhhhh.....................sssssssssss
 *    ^               ^      ^             ^          ^
 *    ma_segment      ma_hp  ma_partition  ma_sp      ma_eos
 */

private Pointer ma_segment;   /* memory from OS */
private Pointer ma_partition; /* boundary between heap and stack */
private Pointer ma_eos;       /* end of segment */
private Pointer ma_hp;        /* heap pointer */
private Pointer ma_sp;        /* stack pointer */

private AD *ma_hfree;         /* free list for heap */
private AD *ma_hused;         /* used list for heap */
private AD *ma_sused;         /* used list for stack */

/* toggled when ma_preinitialize succeeds */
private Boolean ma_preinitialized = MA_FALSE;

/* toggled when MA_init succeeds */
private Boolean ma_initialized = MA_FALSE;

/* invoke MA_verify_allocator_stuff in each public routine? */
private Boolean ma_auto_verify = MA_FALSE;

/* print push/pop/alloc/free? */
private Boolean ma_trace = MA_FALSE;

/* base arrays for the C datatypes */
public char                 ma_cb_char[2];    /* MT_C_CHAR */
public int                  ma_cb_int[2];     /* MT_C_INT */
public long                 ma_cb_long[2];    /* MT_C_LONGINT */
public long long            ma_cb_longlong[2];/* MT_C_LONGLONG */
public float                ma_cb_float[2];   /* MT_C_FLOAT */
public double               ma_cb_dbl[2];     /* MT_C_DBL */
public MA_LongDouble        ma_cb_ldbl[2];    /* MT_C_LDBL */
public MA_SingleComplex     ma_cb_scpl[2];    /* MT_C_SCPL */
public MA_DoubleComplex     ma_cb_dcpl[2];    /* MT_C_DCPL */
public MA_LongDoubleComplex ma_cb_ldcpl[2];   /* MT_C_LDCPL */

/* requested power-of-two alignment */
private Integer ma_numalign = 0;

/**
 ** macros
 **/

/* minimum of two values */
#ifdef min
#undef min
#endif
#define min(a, b) (((b) < (a)) ? (b) : (a))

/* maximum of two values */
#ifdef max
#undef max
#endif
#define max(a, b) (((b) > (a)) ? (b) : (a))

/* proper word ending corresponding to n */
#define plural(n) (((n) == 1) ? "" : "s")

/* convert between internal and external datatype values */
#define mt_import(d) ((d) - MT_BASE)
#define mt_export(d) ((d) + MT_BASE)

/* return nonzero if d is a valid (external) datatype */
#define mt_valid(d) (((d) >= MT_FIRST) && ((d) <= MT_LAST))

/* convert between pointer (address) and equivalent byte address */
#define p2b(p) ((ulongi)(p) * BPA)
#define b2p(b) ((Pointer)((b) / BPA))

/* return nonzero if a is a potentially valid address */
#define reasonable_address(a) (((a) >= ma_segment) && ((a) < ma_eos))

/* worst case bytes of overhead for any block of elements of datatype d */
#define max_block_overhead(d) \
    (BLOCK_OVERHEAD_FIXED + (ma_sizeof[d] - 1) + (ALIGNMENT - 1))

/* compute 0-based index for client_space from AD */
#define client_space_index(ad) \
    ((MA_AccessIndex)((long)((ad)->client_space - ma_base[(ad)->datatype]) / \
                 ma_sizeof[(ad)->datatype]))

/* compute address of guard from AD */
#define guard1(ad) ((Pointer)((ad)->client_space - sizeof(Guard)))
#define guard2(ad) ((Pointer)((ad)->client_space \
                + ((ad)->nelem * ma_sizeof[(ad)->datatype])))

/*
 * When reading or writing guard values, it is necessary to do an
 * explicit byte copy to avoid bus errors caused by guards that
 * are not suitably aligned.
 */

/* copy from guard to value */
#define guard_read(guard, value) bytecopy((guard), (value), sizeof(Guard))

/* copy from value to guard */
#define guard_write(guard, value) bytecopy((value), (guard), sizeof(Guard))

/**
 ** statistics stuff
 **/

#ifdef STATS

/* the number of routines for which calls are counted */
#define NUMROUTINES ((int)FID_MA_verify_allocator_stuff + 1)

/* function identifiers */
typedef enum
{
    FID_MA_alloc_get = 0,
    FID_MA_allocate_heap,
    FID_MA_chop_stack,
    FID_MA_free_heap,
    FID_MA_free_heap_piece,
    FID_MA_get_index,
    FID_MA_get_mbase,
    FID_MA_get_next_memhandle,
    FID_MA_get_numalign,
    FID_MA_get_pointer,
    FID_MA_init,
    FID_MA_initialized,
    FID_MA_init_memhandle_iterator,
    FID_MA_inquire_avail,
    FID_MA_inquire_heap,
    FID_MA_inquire_heap_check_stack,
    FID_MA_inquire_heap_no_partition,
    FID_MA_inquire_stack,
    FID_MA_inquire_stack_check_heap,
    FID_MA_inquire_stack_no_partition,
    FID_MA_pop_stack,
    FID_MA_print_stats,
    FID_MA_push_get,
    FID_MA_push_stack,
    FID_MA_set_auto_verify,
    FID_MA_set_error_print,
    FID_MA_set_hard_fail,
    FID_MA_set_numalign,
    FID_MA_sizeof,
    FID_MA_sizeof_overhead,
    FID_MA_summarize_allocated_blocks,
    FID_MA_trace,
    FID_MA_verify_allocator_stuff
} FID;

/* MA usage statistics */
typedef struct
{
    ulongi    hblocks;           /* current # of heap blocks */
    ulongi    hblocks_max;       /* max # of heap blocks */
    ulongi    hbytes;            /* current # of heap bytes */
    ulongi    hbytes_max;        /* max # of heap bytes */
    ulongi    sblocks;           /* current # of stack blocks */
    ulongi    sblocks_max;       /* max # of stack blocks */
    ulongi    sbytes;            /* current # of stack bytes */
    ulongi    sbytes_max;        /* max # of stack bytes */
    ulongi    calls[NUMROUTINES];/* # of calls to each routine */
} Stats;

/* names of the routines */
private char *ma_routines[] =
{
    "MA_alloc_get",
    "MA_allocate_heap",
    "MA_chop_stack",
    "MA_free_heap",
    "MA_free_heap_piece",
    "MA_get_index",
    "MA_get_mbase",
    "MA_get_next_memhandle",
    "MA_get_numalign",
    "MA_get_pointer",
    "MA_init",
    "MA_initialized",
    "MA_init_memhandle_iterator",
    "MA_inquire_avail",
    "MA_inquire_heap",
    "MA_inquire_heap_check_stack",
    "MA_inquire_heap_no_partition",
    "MA_inquire_stack",
    "MA_inquire_stack_check_heap",
    "MA_inquire_stack_no_partition",
    "MA_pop_stack",
    "MA_print_stats",
    "MA_push_get",
    "MA_push_stack",
    "MA_set_auto_verify",
    "MA_set_error_print",
    "MA_set_hard_fail",
    "MA_set_numalign",
    "MA_sizeof",
    "MA_sizeof_overhead",
    "MA_summarize_allocated_blocks",
    "MA_trace",
    "MA_verify_allocator_stuff"
};

/* MA usage statistics */
private Stats ma_stats;

#endif /* STATS */

/**
 ** private routines
 **/

/* ------------------------------------------------------------------------- */
/*
 * Return MA_TRUE if ad can satisfy ar, else return MA_FALSE.
 * If ad can satisfy ar, set its client_space and nbytes fields
 * after performing any splitting.
 */
/* ------------------------------------------------------------------------- */

private Boolean ad_big_enough(ad, ar)
    AD        *ad;        /* the AD to test */
    Pointer    ar;        /* allocation request */
{
    Pointer    client_space;    /* location of client_space */
    ulongi    nbytes;        /* length of block for ar */

    /* perform trial allocation to determine size */
    balloc_after((AR *)ar, (Pointer)ad, &client_space, &nbytes);

    if (nbytes <= ad->nbytes)
    {
        /* ad is big enough; split block if necessary */
        (void)block_split(ad, nbytes, MA_TRUE);

        /* set fields appropriately */
        ad->client_space = client_space;

        /* success */
        return MA_TRUE;
    }
    else
        /* ad is not big enough */
        return MA_FALSE;
}

/* ------------------------------------------------------------------------- */
/*
 * Return MA_TRUE if ad == ad_target, else return MA_FALSE.
 */
/* ------------------------------------------------------------------------- */

private Boolean ad_eq(ad, ad_target)
    AD        *ad;        /* the AD to test */
    Pointer    ad_target;    /* the AD to match */
{
    return (ad == (AD *)ad_target) ? MA_TRUE : MA_FALSE;
}

/* ------------------------------------------------------------------------- */
/*
 * Return MA_TRUE if ad > ad_target, else return MA_FALSE.
 */
/* ------------------------------------------------------------------------- */

private Boolean ad_gt(ad, ad_target)
    AD        *ad;        /* the AD to test */
    Pointer    ad_target;    /* the AD to match */
{
    return (ad > (AD *)ad_target) ? MA_TRUE : MA_FALSE;
}

/* ------------------------------------------------------------------------- */
/*
 * Return MA_TRUE if ad <= ad_target, else return MA_FALSE.
 */
/* ------------------------------------------------------------------------- */

private Boolean ad_le(ad, ad_target)
    AD        *ad;        /* the AD to test */
    Pointer    ad_target;    /* the AD to match */
{
    return (ad <= (AD *)ad_target) ? MA_TRUE : MA_FALSE;
}

/* ------------------------------------------------------------------------- */
/*
 * Return MA_TRUE if ad < ad_target, else return MA_FALSE.
 */
/* ------------------------------------------------------------------------- */

private Boolean ad_lt(ad, ad_target)
    AD        *ad;        /* the AD to test */
    Pointer    ad_target;    /* the AD to match */
{
    return (ad < (AD *)ad_target) ? MA_TRUE : MA_FALSE;
}

/* ------------------------------------------------------------------------- */
/*
 * Print identifying information about the given AD to stdout.
 */
/* ------------------------------------------------------------------------- */

private void ad_print(ad, block_type)
    AD        *ad;        /* to print */
    char    *block_type;    /* for output */
{
    Integer    memhandle;    /* memhandle for AD */

    /* convert AD to memhandle */
    memhandle = ma_table_lookup_assoc((TableData)ad);

    /* print to stdout */
    (void)printf("%s block '%s', handle ",
        block_type,
        ad->name);
    if (memhandle == TABLE_HANDLE_NONE)
        (void)printf("unknown");
    else
        (void)printf("%ld",
            (long)memhandle);
    (void)printf(", address 0x%lx",
        (long)ad);
}

/* ------------------------------------------------------------------------- */
/*
 * Allocate a block suitable for ar starting at address.  No fields of
 * the new block are modified.
 */
/* ------------------------------------------------------------------------- */

private void balloc_after(ar, address, client_space, nbytes)
    AR        *ar;        /* allocation request */
    Pointer    address;    /* to allocate after */
    Pointer    *client_space;    /* RETURN: location of client_space */
    ulongi    *nbytes;    /* RETURN: length of block */
{
    Integer    datatype;    /* of elements in this block */
    ulongi    L_client_space;    /* length of client_space */
    Pointer    A_client_space;    /* address of client_space */
    int        L_gap1;        /* length of gap1 */
    int        L_gap2;        /* length of gap2 */

    ulongi    B_address;    /* byte equivalent of address */
    ulongi    B_base;        /* byte equivalent of ma_base[datatype] */
    ulongi    B_client_space;    /* byte equivalent of A_client_space */

    datatype = ar->datatype;

    B_address = p2b(address);
    B_base = p2b(ma_base[datatype]);

    /*
     * To ensure that client_space is properly aligned:
     *
     *    (A(client_space) - ma_base[datatype]) % ma_sizeof[datatype] == 0
     *
     * where
     *
     *    A(client_space) == address + L(AD) + L(gap1) + L(guard1)
     */

    L_client_space = ar->nelem * ma_sizeof[datatype];

    L_gap1 = ((long)B_base
        - (long)B_address
        - (long)sizeof(AD)
        - (long)sizeof(Guard))
        % (long)ma_sizeof[datatype];

    if (L_gap1 < 0)
        L_gap1 += ma_sizeof[datatype];

    B_client_space = B_address + sizeof(AD) + L_gap1 + sizeof(Guard);
    A_client_space = b2p(B_client_space);
    B_client_space = p2b(A_client_space);

    /*
     * To align client space according to overall alignment of absolute
     * address on user requested 2^ma_numalign boundary.
     * Note that if the base arrays are not aligned accordingly then
     * this alignement request is not satisfiable and will be quietly
     * ignored.
     */

    if (ma_numalign > 0) {
      unsigned long mask = (1<<ma_numalign)-1;
      int diff = ((unsigned long) B_client_space) & mask;
      
      /* Check that the difference is a multiple of the type size.
       * If so, then we can shift the client space which is already
       * aligned to satisfy this requirement.
       */

      if (diff) {
    diff = (1<<ma_numalign) - diff;
    if ((diff % ma_sizeof[datatype]) == 0 ) {
      /*printf("bafter realigned diff=%d\n",diff);*/
      A_client_space = b2p(B_client_space + diff);
      B_client_space = p2b(A_client_space);
    }    
    /*    else {
      printf("did not realign diff=%d typelen=%d mod=%d\n",
         diff, ma_sizeof[datatype], (diff % ma_sizeof[datatype]));
         }*/
      }
    }

    /*
     * To ensure that the AD is properly aligned:
     *
     *    L(block) % ALIGNMENT == 0
     *
     * where
     *
     *    L(block) == A(client_space) + L(client_space) + L(guard2) + L(gap2)
     *        - address
     */

    L_gap2 = ((long)B_address
        - (long)B_client_space
        - (long)L_client_space
        - (long)sizeof(Guard))
        % (long)ALIGNMENT;

    if (L_gap2 < 0)
        L_gap2 += ALIGNMENT;

    /*
     * set the return values
     */

    *client_space = A_client_space;
    *nbytes = (ulongi)(B_client_space
        + L_client_space
        + sizeof(Guard)
        + L_gap2
        - B_address);
}

/* ------------------------------------------------------------------------- */
/*
 * Allocate a block suitable for ar ending before address.  No fields of
 * the new block are modified.
 */
/* ------------------------------------------------------------------------- */

private void balloc_before(ar, address, client_space, nbytes)
    AR        *ar;        /* allocation request */
    Pointer    address;    /* to allocate before */
    Pointer    *client_space;    /* RETURN: location of client_space */
    ulongi    *nbytes;    /* RETURN: length of block */
{
    Integer    datatype;    /* of elements in this block */
    ulongi    L_client_space;    /* length of client_space */
    Pointer    A_client_space;    /* address of client_space */
    int        L_gap1;        /* length of gap1 */
    int        L_gap2;        /* length of gap2 */

    ulongi    B_address;    /* byte equivalent of address */
    ulongi    B_base;        /* byte equivalent of ma_base[datatype] */
    ulongi    B_client_space;    /* byte equivalent of A_client_space */

    datatype = ar->datatype;

    B_address = p2b(address);
    B_base = p2b(ma_base[datatype]);

    /*
     * To ensure that client_space is properly aligned:
     *
     *    (A(client_space) - ma_base[datatype]) % ma_sizeof[datatype] == 0
     *
     * where
     *
     *    A(client_space) == address - L(gap2) - L(guard2) - L(client_space)
     */

    L_client_space = ar->nelem * ma_sizeof[datatype];

    L_gap2 = (B_address
        - sizeof(Guard)
        - L_client_space
        - B_base)
        % ma_sizeof[datatype];

    if (L_gap2 < 0)
        L_gap2 += ma_sizeof[datatype];

    B_client_space = B_address - L_gap2 - sizeof(Guard) - L_client_space;
    A_client_space = b2p(B_client_space);
    B_client_space = p2b(A_client_space);

    /*
     * To align client space according to overall alignment of absolute
     * address on user requested 2^ma_numalign boundary.
     * Note that if the base arrays are not aligned accordingly then
     * this alignement request is not satisfiable and will be quietly
     * ignored.
     */

    if (ma_numalign > 0) {
      unsigned long mask = (1<<ma_numalign)-1;
      int diff = ((unsigned long) B_client_space) & mask;
      
      /* Check that the difference is a multiple of the type size.
       * If so, then we can shift the client space which is already
       * aligned to satisfy this requirement.
       */

      if (diff) {
    if ((diff % ma_sizeof[datatype]) == 0 ) {
      /* printf("bbefore realigned diff=%d\n",diff); */
      A_client_space = b2p(B_client_space - diff);
      B_client_space = p2b(A_client_space);
    }    
    /*    else {
      printf("did not realign diff=%d typelen=%d mod=%d\n",
         diff, ma_sizeof[datatype], (diff % ma_sizeof[datatype]));
         }*/
      }
    }

    /*
     * To ensure that the AD is properly aligned:
     *
     *    A(AD) % ALIGNMENT == 0
     *
     * where
     *
     *    A(AD) == A(client_space) - L(guard1) - L(gap1) - L(AD)
     */

    L_gap1 = (B_client_space
        - sizeof(Guard)
        - sizeof(AD))
        % ALIGNMENT;

    /*
     * set the return values
     */

    *client_space = A_client_space;
    *nbytes = (ulongi)(B_address
        - B_client_space
        + sizeof(Guard)
        + L_gap1
        + sizeof(AD));
}

/* ------------------------------------------------------------------------- */
/*
 * Reclaim the given block by updating ma_hp and ma_hfree.
 */
/* ------------------------------------------------------------------------- */

private void block_free_heap(ad)
    AD        *ad;        /* AD to free */
{
    AD        *ad2;        /* traversal pointer */
    AD        *max_ad;    /* rightmost AD */

    /* find rightmost heap block */
    for (max_ad = (AD *)NULL, ad2 = ma_hused; ad2; ad2 = ad2->next)
    {
        if (ad2 > max_ad)
            max_ad = ad2;
    }

    if (max_ad)
    {
        /* at least 1 block is in use */

        /* set ma_hp to first address past end of max_ad */
        ma_hp = (Pointer)max_ad + max_ad->nbytes;

        /* delete any free list blocks that are no longer in heap region */
        (void)list_delete_many(
            &ma_hfree,
            ad_gt,
            (Pointer)max_ad,
            (void (*)())NULL);

        /* if ad is in the heap region, add it to free list */
        if (ad < max_ad)
        {
            list_insert_ordered(ad, &ma_hfree, ad_lt);
            list_coalesce(ma_hfree);
        }
    }
    else
    {
        /* no blocks are in use */

        /* set ma_hp to start of segment */
        ma_hp = ma_segment;

        /* clear the free list */
        ma_hfree = (AD *)NULL;
    }
}

/* ------------------------------------------------------------------------- */
/*
 * If ad is sufficiently bigger than bytes_needed bytes, create a new
 * block from the remainder, optionally insert it in the free list,
 * and set the lengths of both blocks.
 *
 * Return a pointer to the new block (NULL if not created).
 */
/* ------------------------------------------------------------------------- */

private AD *block_split(ad, bytes_needed, insert_free)
    AD        *ad;        /* the AD to split */
    ulongi    bytes_needed;    /* from ad */
    Boolean    insert_free;    /* insert in free list? */
{
    ulongi    bytes_extra;    /* in ad */
    AD        *ad2;        /* the new AD */

    /* caller ensures that ad->nbytes >= bytes_needed */
    bytes_extra = ad->nbytes - bytes_needed;

    if (bytes_extra >= ((ulongi)MINBLOCKSIZE))
    {
        /* create a new block */
        ad2 = (AD *)((Pointer)ad + bytes_needed);

        /* set the length of ad2 */
        ad2->nbytes = bytes_extra;

        if (insert_free)
        {
            /* insert ad2 into free list */
            list_insert_ordered(ad2, &ma_hfree, ad_lt);
        }

        /* set the length of ad */
        ad->nbytes = bytes_needed;

        return ad2;
    }
    else
    {
        /*
         * If 0 <= bytes_extra < MINBLOCKSIZE then there are too few
         * extra bytes to form a new block.  In this case, we simply
         * do nothing; ad will retain its original length (which is
         * either perfect or slightly too big), and the entire block
         * will be reclaimed upon deallocation, preventing any
         * memory leakage.
         */

        return (AD *)NULL;
    }
}

/* ------------------------------------------------------------------------- */
/*
 * Compute and return a checksum for ad.  Include all fields except name,
 * next, and checksum.
 */
/* ------------------------------------------------------------------------- */

private ulongi checksum(ad)
    AD        *ad;        /* the AD to compute checksum for */
{
    return (ulongi)(
                ad->datatype +
                ad->nelem +
        (ulongi)ad->client_space +
                ad->nbytes);
}

/* ------------------------------------------------------------------------- */
/*
 * Print to stderr the addresses of the fields of the given ad.
 */
/* ------------------------------------------------------------------------- */

#ifdef DEBUG

private void debug_ad_print(ad)
    AD        *ad;        /* the AD to print */
{
#define NUMADFIELDS 7

    char    *fn[NUMADFIELDS];    /* field names */
    long    fa[NUMADFIELDS];    /* field addresses */
    int        i;            /* loop index */
    long    address;        /* other addresses */

    /* set field names */
    fn[0] = "datatype";
    fn[1] = "nelem";
    fn[2] = "name";
    fn[3] = "client_space";
    fn[4] = "nbytes";
    fn[5] = "next";
    fn[6] = "checksum";

    /* set field addresses */
    fa[0] = (long)(&(ad->datatype));
    fa[1] = (long)(&(ad->nelem));
    fa[2] = (long)(&(ad->name));
    fa[3] = (long)(&(ad->client_space));
    fa[4] = (long)(&(ad->nbytes));
    fa[5] = (long)(&(ad->next));
    fa[6] = (long)(&(ad->checksum));

    /* print AD fields to stderr */
    (void)fprintf(stderr, "debug_ad_print:\n");
    for (i = 0; i < NUMADFIELDS; i++)
        (void)fprintf(stderr, "\t0x%lx  mod4,8,16=%d,%d,%-2d  ad->%s\n",
            fa[i],
            fa[i] % 4,
            fa[i] % 8,
            fa[i] % 16,
            fn[i]);

    /* print other addresses to stderr */
    address = (long)guard1(ad);
    (void)fprintf(stderr, "\t0x%lx  mod4,8,16=%d,%d,%-2d  guard1\n",
        address,
        address % 4,
        address % 8,
        address % 16);
    address = (long)ad->client_space;
    (void)fprintf(stderr, "\t0x%lx  mod4,8,16=%d,%d,%-2d  client_space\n",
        address,
        address % 4,
        address % 8,
        address % 16);
    address = (long)guard2(ad);
    (void)fprintf(stderr, "\t0x%lx  mod4,8,16=%d,%d,%-2d  guard2\n",
        address,
        address % 4,
        address % 8,
        address % 16);

    (void)fflush(stderr);
}

#endif /* DEBUG */

/* ------------------------------------------------------------------------- */
/*
 * Return MA_TRUE if the guards associated with ad contain valid signatures,
 * else return MA_FALSE.
 */
/* ------------------------------------------------------------------------- */

private Boolean guard_check(ad)
    AD        *ad;        /* the AD to check guards for */
{
    Guard    signature;    /* value to be read */
    Pointer    guard;        /* address to read from */

    guard = guard1(ad);
    guard_read(guard, &signature);
    if (signature != GUARD1)
        return MA_FALSE;

    guard = guard2(ad);
    guard_read(guard, &signature);
    if (signature != GUARD2)
        return MA_FALSE;

    /* success */
    return MA_TRUE;
}

/* ------------------------------------------------------------------------- */
/*
 * Write signatures into the guards associated with ad.
 */
/* ------------------------------------------------------------------------- */

private void guard_set(ad)
    AD        *ad;        /* the AD to set guards for */
{
    Guard    signature;    /* value to be written */
    Pointer    guard;        /* address to write to */

    signature = GUARD1;
    guard = guard1(ad);
    guard_write(guard, &signature);

    signature = GUARD2;
    guard = guard2(ad);
    guard_write(guard, &signature);
}

/* ------------------------------------------------------------------------- */
/*
 * Coalesce list by merging any adjacent elements that are contiguous.
 * The list is assumed to be ordered by increasing addresses, i.e.,
 * addressOf(element i) < addressOf(element i+1).
 */
/* ------------------------------------------------------------------------- */

private void list_coalesce(list)
    AD        *list;        /* the list to coalesce */
{
    AD        *ad1;        /* lead traversal pointer */
    AD        *ad2;        /* trailing traversal pointer */

    for (ad2 = list; ad2;)
    {
        /* compute first address beyond ad2 */
        ad1 = (AD *)((Pointer)ad2 + ad2->nbytes);

        /* are ad2 and ad1 contiguous? */
        if (ad1 == ad2->next)
        {
            /* yes; merge ad1 into ad2 */
            ad2->nbytes += ad1->nbytes;
            ad2->next = ad1->next;
        }
        else
        {
            /* no; advance ad2 */
            ad2 = ad2->next;
        }
    }
}

/* ------------------------------------------------------------------------- */
/*
 * Delete and return the first occurrence of ad from list.  If ad is not
 * a member of list, return NULL.
 */
/* ------------------------------------------------------------------------- */

private AD *list_delete(ad, list)
    AD        *ad;        /* the AD to delete */
    AD        **list;        /* the list to delete from */
{
    return list_delete_one(list, ad_eq, (Pointer)ad);
}

/* ------------------------------------------------------------------------- */
/*
 * Apply pred (with closure) to each element of list.  Delete each element
 * that satisfies pred, after applying action to the element (if action is
 * not NULL).  Return the number of elements deleted.
 */
/* ------------------------------------------------------------------------- */

private int list_delete_many(list, pred, closure, action)
    AD        **list;        /* the list to search */
    Boolean    (*pred)();    /* predicate */
    Pointer    closure;    /* for pred */
    void    (*action)();    /* to apply before deletion */
{
    AD        *ad1;        /* lead traversal pointer */
    AD        *ad2;        /* trailing traversal pointer */
    int        ndeleted = 0;    /* # of elements deleted from list */

    for (ad2 = (AD *)NULL, ad1 = *list; ad1; ad1 = ad1->next)
    {
        /* does ad1 match? */
        if ((*pred)(ad1, closure))
        {
            /* yes; apply action, then delete ad1 from list */
            if (action != (void (*)())NULL)
                (*action)(ad1);
            if (ad2)
            {
                /* ad1 is second or later element */
                ad2->next = ad1->next;
            }
            else
            {
                /* ad1 is first element */
                *list = ad1->next;
            }

            ndeleted++;
        }
        else
        {
            /* no; ad1 survives, so scoot ad2 along */
            ad2 = ad1;
        }
    }

    /* return the # of elements deleted from list */
    return ndeleted;
}

/* ------------------------------------------------------------------------- */
/*
 * Apply pred (with closure) to each element of list.  Delete and return
 * the first element that satisfies pred.  If no element satisfies pred,
 * return NULL.
 */
/* ------------------------------------------------------------------------- */

private AD *list_delete_one(list, pred, closure)
    AD        **list;        /* the list to search */
    Boolean    (*pred)();    /* predicate */
    Pointer    closure;    /* for pred */
{
    AD        *ad1;        /* lead traversal pointer */
    AD        *ad2;        /* trailing traversal pointer */

    for (ad2 = (AD *)NULL, ad1 = *list; ad1; ad2 = ad1, ad1 = ad1->next)
    {
        /* does ad1 match? */
        if ((*pred)(ad1, closure))
        {
            /* yes; delete ad1 from list */
            if (ad2)
            {
                /* ad1 is second or later element */
                ad2->next = ad1->next;
            }
            else
            {
                /* ad1 is first element */
                *list = ad1->next;
            }

            /* success */
            return ad1;
        }
    }

    /* failure */
    return (AD *)NULL;
}

/* ------------------------------------------------------------------------- */
/*
 * Insert ad into list.
 */
/* ------------------------------------------------------------------------- */

private void list_insert(ad, list)
    AD        *ad;        /* the AD to insert */
    AD        **list;        /* the list to insert into */
{
    /* push ad onto list */
    ad->next = *list;
    *list = ad;
}

/* ------------------------------------------------------------------------- */
/*
 * Insert ad into list, immediately before the first element e
 * for which pred(ad, e) returns true.  If there is no such element,
 * insert ad after the last element of list.
 */
/* ------------------------------------------------------------------------- */

private void list_insert_ordered(ad, list, pred)
    AD        *ad;        /* the AD to insert */
    AD        **list;        /* the list to insert into */
    Boolean    (*pred)();    /* predicate */
{
    AD        *ad1;        /* lead traversal pointer */
    AD        *ad2;        /* trailing traversal pointer */

    if (*list == (AD *)NULL)
    {
        /* empty list */
        ad->next = (AD *)NULL;
        *list = ad;
        return;
    }

    /* list has at least one element */
    for (ad2 = (AD *)NULL, ad1 = *list; ad1; ad2 = ad1, ad1 = ad1->next)
    {
        /* does ad1 match? */
        if ((*pred)(ad, ad1))
        {
            /* yes; insert ad before ad1 */
            if (ad2)
            {
                /* ad1 is second or later element */
                ad2->next = ad;
            }
            else
            {
                /* ad1 is first element */
                *list = ad;
            }
            ad->next = ad1;

            /* success */
            return;
        }
    }

    /* append ad to list */
    ad2->next = ad;
    ad->next = (AD *)NULL;
}

/* ------------------------------------------------------------------------- */
/*
 * Return MA_TRUE if ad is a member of list, else return MA_FALSE.
 */
/* ------------------------------------------------------------------------- */

private Boolean list_member(ad, list)
    AD        *ad;        /* the AD to search for */
    AD        *list;        /* the list to search */
{
    AD        *ad1;        /* traversal pointer */

    for (ad1 = list; ad1; ad1 = ad1->next)
        if (ad1 == ad)
            /* ad is a member of list */
            return MA_TRUE;

    /* ad is not a member of list */
    return MA_FALSE;
}

/* ------------------------------------------------------------------------- */
/*
 * Print information to stdout about each block on list.  Return the
 * number of blocks on list.
 */
/* ------------------------------------------------------------------------- */

private int list_print(list, block_type, index_base)
    AD        *list;        /* to print */
    char    *block_type;    /* for output */
    int        index_base;    /* 0 (C) or 1 (FORTRAN) */
{
    AD        *ad;        /* traversal pointer */
    int        nblocks;    /* # of blocks on list */

    /* print each block on list */
    for (ad = list, nblocks = 0; ad; ad = ad->next, nblocks++)
    {
        /* print to stdout */
        ad_print(ad, block_type);
        (void)printf(":\n");
        (void)printf("\ttype of elements:\t\t%s\n",
            ma_datatype[ad->datatype]);
        (void)printf("\tnumber of elements:\t\t%ld\n",
            (long)ad->nelem);
        (void)printf("\taddress of client space:\t0x%lx\n",
            (long)ad->client_space);
        (void)printf("\tindex for client space:\t\t%ld\n",
            (long)(client_space_index(ad) + index_base));
        (void)printf("\ttotal number of bytes:\t\t%lu\n",
            ad->nbytes);
    }

    /* return the number of blocks on list */
    return nblocks;
}

/* ------------------------------------------------------------------------- */
/*
 * Check each block on list for checksum and guard errors.  For each error
 * found, print a message to stdout.  Return counts of the various errors
 * in the bad_ parameters.
 */
/* ------------------------------------------------------------------------- */

private void list_verify(list, block_type, preamble, blocks,
                         bad_blocks, bad_checksums, bad_lguards, bad_rguards)
    AD        *list;        /* to verify */
    char    *block_type;    /* for error messages */
    char    *preamble;    /* printed before first error message */
    int        *blocks;    /* RETURN: # of blocks */
    int        *bad_blocks;    /* RETURN: # of blocks having errors */
    int        *bad_checksums;    /* RETURN: # of blocks having bad checksum */
    int        *bad_lguards;    /* RETURN: # of blocks having bad guard1 */
    int        *bad_rguards;    /* RETURN: # of blocks having bad guard2 */
{
    AD        *ad;        /* traversal pointer */
    Boolean    first_bad_block;/* first bad block found? */
    Boolean    bad_block;    /* problem in current block? */
    Guard    signature;    /* value to be read */
    Pointer    guard;        /* address to read from */

    /* initialize */
    *blocks = 0;
    *bad_blocks = 0;
    *bad_checksums = 0;
    *bad_lguards = 0;
    *bad_rguards = 0;
    first_bad_block = MA_TRUE;

    /* check each block on list */
    for (ad = list; ad; ad = ad->next)
    {
        (*blocks)++;
        bad_block = MA_FALSE;

        /* check for checksum error */
        if (checksum(ad) != ad->checksum)
        {
            /* print preamble if necessary */
            if (first_bad_block && (preamble != (char *)NULL))
            {
                (void)printf(preamble);
                first_bad_block = MA_FALSE;
            }

            /* print error message to stdout */
            ad_print(ad, block_type);
            (void)printf(":\n\t");
            (void)printf("current checksum %lu != stored checksum %lu\n",
                checksum(ad),
                ad->checksum);

            /* do bookkeeping */
            (*bad_checksums)++;
            bad_block = MA_TRUE;
        }

        /* check for bad guard1 */
        guard = guard1(ad);
        guard_read(guard, &signature);
        if (signature != GUARD1)
        {
            /* print preamble if necessary */
            if (first_bad_block && (preamble != (char *)NULL))
            {
                (void)printf(preamble);
                first_bad_block = MA_FALSE;
            }

            /* print error message to stdout */
            ad_print(ad, block_type);
            (void)printf(":\n\t");
            (void)printf(
                "current left signature %u != proper left signature %u\n",
                signature,
                GUARD1);

            /* do bookkeeping */
            (*bad_lguards)++;
            bad_block = MA_TRUE;
        }

        /* check for bad guard2 */
        guard = guard2(ad);
        guard_read(guard, &signature);
        if (signature != GUARD2)
        {
            /* print preamble if necessary */
            if (first_bad_block && (preamble != (char *)NULL))
            {
                (void)printf(preamble);
                first_bad_block = MA_FALSE;
            }

            /* print error message to stdout */
            ad_print(ad, block_type);
            (void)printf(":\n\t");
            (void)printf(
                "current right signature %u != proper right signature %u\n",
                signature,
                GUARD2);

            /* do bookkeeping */
            (*bad_rguards)++;
            bad_block = MA_TRUE;
        }

        /* if any errors, bump bad block count */
        if (bad_block)
            (*bad_blocks)++;
    }
}

/* ------------------------------------------------------------------------- */
/*
 * Return the maximum number of datatype elements that can currently be
 * accomodated in a heap fragment (a block on the heap free list) entirely
 * within the heap region, or 0 if this number is less than min_nelem.
 */
/* ------------------------------------------------------------------------- */

private Integer ma_max_heap_frag_nelem(datatype, min_nelem)
    Integer    datatype;    /* of elements */
    Integer    min_nelem;    /* for fragment to be considered */
{
    ulongi    min_bytes;    /* for fragment to be considered */
    AD        *ad;        /* traversal pointer */
    ulongi    nbytes;        /* in current fragment */
    Integer    nelem;        /* in current fragment */
    Integer    max_nelem;    /* result */

    /* set the threshold */
    min_bytes = (min_nelem * ma_sizeof[datatype]) + BLOCK_OVERHEAD_FIXED;

    /* search the heap free list */
    max_nelem = 0;
    for (ad = ma_hfree; ad; ad = ad->next)
    {
        /*
         * There are 3 cases to consider:
         *
         * (a) fragment is outside heap region
         * (b) fragment straddles partition between heap and stack regions
         * (c) fragment is inside heap region
         */

        if ((Pointer)ad >= ma_partition)
        {
            /* case (a): reject */
            continue;
        }
        else if (((Pointer)ad + ad->nbytes) >= ma_partition)
        {
            /* case (b): truncate fragment at partition */
            nbytes = (ulongi)(ma_partition - (Pointer)ad);
        }
        else
        {
            /* case (c): accept */
            nbytes = ad->nbytes;
        }

        if (nbytes >= min_bytes)
        {
            nelem = ma_nelem((Pointer)ad, nbytes, datatype);
            max_nelem = max(max_nelem, nelem);
        }
    }

    /* return the result */
    return max_nelem;
}

/* ------------------------------------------------------------------------- */
/*
 * Return the maximum number of datatype elements that can currently
 * be accomodated in length bytes starting at address.
 */
/* ------------------------------------------------------------------------- */

private Integer ma_nelem(address, length, datatype)
    Pointer    address;    /* location of hypothetical block */
    ulongi    length;        /* length of hypothetical block */
    Integer    datatype;    /* of elements in hypothetical block */
{
    AR        ar;        /* allocation request */
    Pointer    client_space;    /* location of client_space */
    ulongi    nbytes;        /* length of block for ar */

    if (length <= BLOCK_OVERHEAD_FIXED)
        /* no point in computing anything */
        return (Integer)0;

    /* compute initial request */
    ar.datatype = datatype;
    ar.nelem = (length - BLOCK_OVERHEAD_FIXED) / ma_sizeof[datatype];

    /* make requests until one succeeds or we give up */
    while (ar.nelem > 0)
    {
        /* perform trial allocation to determine size */
        balloc_after(&ar, address, &client_space, &nbytes);

        if (nbytes > length)
            /* not enough space for ar.nelem elements */
            ar.nelem--;
        else
            /* enough space for ar.nelem elements */
            break;
    }

    /* return the result */
    return ar.nelem;
}

/* ------------------------------------------------------------------------- */
/*
 * Perform operations necessary to allow certain functions to be invoked
 * before MA_init is called.
 */
/* ------------------------------------------------------------------------- */

private void ma_preinitialize(caller)
    char    *caller;    /* name of calling routine */
{
    if (ma_preinitialized)
        return;

    /* call a FORTRAN routine to set bases and sizes of FORTRAN datatypes */
    if (ma_set_sizes_() == 0)
    {
        (void)sprintf(ma_ebuf,
            "unable to set sizes of FORTRAN datatypes");
        ma_error(EL_Fatal, ET_Internal, caller, ma_ebuf);
        return;
    }

    /* success */
    ma_preinitialized = MA_TRUE;
}

/* ------------------------------------------------------------------------- */
/*
 * If memhandle is valid according to location, return the corresponding AD
 * in adout and return MA_TRUE, else print an error message and return
 * MA_FALSE.
 */
/* ------------------------------------------------------------------------- */

private Boolean mh2ad(memhandle, adout, location, caller)
    Integer    memhandle;    /* the handle to verify and convert */
    AD        **adout;    /* RETURN: AD corresponding to memhandle */
    BlockLocation location;    /* where AD must reside */
    char    *caller;    /* name of calling routine */
{
    AD        *ad;
    Boolean    check_checksum = MA_TRUE;
    Boolean    check_guards = MA_TRUE;
    Boolean    check_heap = MA_FALSE;
    Boolean    check_stack = MA_FALSE;
    Boolean    check_stacktop = MA_FALSE;
    Boolean    check_heapandstack = MA_FALSE;

    switch (location)
    {
        case BL_HeapOrStack:
            check_heapandstack = MA_TRUE;
            break;
        case BL_Heap:
            check_heap = MA_TRUE;
            break;
        case BL_Stack:
            check_stack = MA_TRUE;
            break;
        case BL_StackTop:
            check_stacktop = MA_TRUE;
            break;
        default:
            (void)sprintf(ma_ebuf,
                "invalid location: %d",
                (int)location);
            ma_error(EL_Nonfatal, ET_Internal, "mh2ad", ma_ebuf);
            return MA_FALSE;
    }

    /* convert memhandle to AD */
    if (!ma_table_verify(memhandle, caller))
        return MA_FALSE;
    else
        ad = (AD *)ma_table_lookup(memhandle);

    /* attempt to avoid crashes due to corrupt addresses */
    if (!reasonable_address((Pointer)ad))
    {
        (void)sprintf(ma_ebuf,
            "invalid block address (0x%lx) for memhandle %ld",
            (long)ad, (long)memhandle);
        ma_error(EL_Nonfatal, ET_External, caller, ma_ebuf);
        return MA_FALSE;
    }

    if (check_checksum)
    {
        if (checksum(ad) != ad->checksum)
        {
            (void)sprintf(ma_ebuf,
                "invalid checksum for memhandle %ld (name: '%s')",
                (long)memhandle, ad->name);
            ma_error(EL_Nonfatal, ET_External, caller, ma_ebuf);
            return MA_FALSE;
        }
    }

    if (check_guards)
    {
        if (!guard_check(ad))
        {
            (void)sprintf(ma_ebuf,
                "invalid guard(s) for memhandle %ld (name: '%s')",
                (long)memhandle, ad->name);
            ma_error(EL_Nonfatal, ET_External, caller, ma_ebuf);
            return MA_FALSE;
        }
    }

    if (check_heap)
    {
        if (!list_member(ad, ma_hused))
        {
            (void)sprintf(ma_ebuf,
                "memhandle %ld (name: '%s') not in heap",
                (long)memhandle, ad->name);
            ma_error(EL_Nonfatal, ET_External, caller, ma_ebuf);
            return MA_FALSE;
        }
    }
    else if (check_stack)
    {
        if (!list_member(ad, ma_sused))
        {
            (void)sprintf(ma_ebuf,
                "memhandle %ld (name: '%s') not in stack",
                (long)memhandle, ad->name);
            ma_error(EL_Nonfatal, ET_External, caller, ma_ebuf);
            return MA_FALSE;
        }
    }
    else if (check_stacktop)
    {
        /* is it in the stack? */
        if (!list_member(ad, ma_sused))
        {
            (void)sprintf(ma_ebuf,
                "memhandle %ld (name: '%s') not in stack",
                (long)memhandle, ad->name);
            ma_error(EL_Nonfatal, ET_External, caller, ma_ebuf);
            return MA_FALSE;
        }

        /* is it on top of the stack? */
        if ((Pointer)ad != ma_sp)
        {
            (void)sprintf(ma_ebuf,
                "memhandle %ld (name: '%s') not top of stack",
                (long)memhandle, ad->name);
            ma_error(EL_Nonfatal, ET_External, caller, ma_ebuf);
            return MA_FALSE;
        }
    }
    else if (check_heapandstack)
    {
        if ((!list_member(ad, ma_hused)) && (!list_member(ad, ma_sused)))
        {
            (void)sprintf(ma_ebuf,
                "memhandle %ld (name: '%s') not in heap or stack",
                (long)memhandle, ad->name);
            ma_error(EL_Nonfatal, ET_External, caller, ma_ebuf);
            return MA_FALSE;
        }
    }

    /* ad is valid */
    *adout = ad;
    return MA_TRUE;
}

/* ------------------------------------------------------------------------- */
/*
 * Free the memhandle corresponding to the given AD.
 */
/* ------------------------------------------------------------------------- */

private void mh_free(ad)
    AD        *ad;        /* the AD whose memhandle to free */
{
    Integer    memhandle;    /* memhandle for AD */

    /* convert AD to memhandle */
    if ((memhandle = ma_table_lookup_assoc((TableData)ad)) == TABLE_HANDLE_NONE)
    {
        (void)sprintf(ma_ebuf,
            "cannot find memhandle for block address 0x%lx",
            (long)ad);
        ma_error(EL_Nonfatal, ET_Internal, "mh_free", ma_ebuf);
    }
    else
        /* free memhandle */
        ma_table_deallocate(memhandle);
}

/* ------------------------------------------------------------------------- */
/*
 * Return the first multiple of unit which is >= value.
 */
/* ------------------------------------------------------------------------- */

private long mai_round(value, unit)
    long    value;        /* to round */
    ulongi    unit;        /* to round to */
{
    /* voodoo ... */
    unit--;
    value += unit;
    value &= ~(long)unit;
    return value;
}

/* ------------------------------------------------------------------------- */
/*
 * Copy at most maxchars-1 non-NUL chars from from to to; NUL-terminate to.
 */
/* ------------------------------------------------------------------------- */

private void str_ncopy(to, from, maxchars)
    char    *to;        /* space to copy to */
    char    *from;        /* space to copy from */
    int        maxchars;    /* max # of chars (including NUL) copied */
{
    if (from == (char *)NULL)
    {
        to[0] = '\0';
        return;
    }

    /* copy up to maxchars chars */
    (void)strncpy(to, from, maxchars);

    /* ensure to is NUL-terminated */
    to[maxchars-1] = '\0';
}

/**
 ** public routines for internal use only
 **/

/* ------------------------------------------------------------------------- */
/*
 * Set the base address and size of the given datatype.
 */
/* ------------------------------------------------------------------------- */

public Boolean MAi_inform_base(datatype, address1, address2)
    Integer    datatype;    /* to set size of */
    Pointer    address1;    /* of datatype element base */
    Pointer    address2;    /* of an adjacent datatype element */
{
    /* verify uninitialization */
    if (ma_initialized)
    {
        (void)sprintf(ma_ebuf,
            "MA already initialized");
        ma_error(EL_Nonfatal, ET_Internal, "MAi_inform_base", ma_ebuf);
        return MA_FALSE;
    }

    /* verify datatype */
    if (!mt_valid(datatype))
    {
        (void)sprintf(ma_ebuf,
            "invalid datatype: %ld",
            (long)datatype);
        ma_error(EL_Nonfatal, ET_Internal, "MAi_inform_base", ma_ebuf);
        return MA_FALSE;
    }

    /* convert datatype to internal (index-suitable) value */
    datatype = mt_import(datatype);

    /* set the base address of datatype */
    ma_base[datatype] = address1;

    /* set the size of datatype */
    ma_sizeof[datatype] = (int)(address2 - address1);

    /* success */
    return MA_TRUE;
}

#if NOFORT
Integer ma_set_sizes_()
{
    MAi_inform_base(MT_F_BYTE, &ma_cb_char[0],  &ma_cb_char[1]);
    MAi_inform_base(MT_F_INT,  &ma_cb_int[0],   &ma_cb_int[1]);
    MAi_inform_base(MT_F_LOG,  &ma_cb_int[0],   &ma_cb_int[1]);
    MAi_inform_base(MT_F_REAL, &ma_cb_float[0], &ma_cb_float[1]);
    MAi_inform_base(MT_F_DBL,  &ma_cb_dbl[0],   &ma_cb_dbl[1]);
    MAi_inform_base(MT_F_SCPL, &ma_cb_scpl[0],  &ma_cb_scpl[1]);
    MAi_inform_base(MT_F_DCPL, &ma_cb_dcpl[0],  &ma_cb_dcpl[1]);
    return 1;
}
#endif

/* ------------------------------------------------------------------------- */
/*
 * Print debugging information about all blocks currently in use
 * on the heap or the stack.
 */
/* ------------------------------------------------------------------------- */

public void MAi_summarize_allocated_blocks(index_base)
    int        index_base;    /* 0 (C) or 1 (FORTRAN) */
{
    int        heap_blocks;    /* # of blocks on heap used list */
    int        stack_blocks;    /* # of blocks on stack used list */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_summarize_allocated_blocks]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return;
#endif /* VERIFY */

    /* verify index_base */
    if ((index_base != 0) && (index_base != 1))
    {
        (void)sprintf(ma_ebuf,
            "invalid index_base: %d",
            index_base);
        ma_error(EL_Nonfatal, ET_Internal, "MAi_summarize_allocated_blocks", ma_ebuf);
        return;
    }

    (void)printf("MA_summarize_allocated_blocks: starting scan ...\n");

    /* print blocks on the heap used list */
    heap_blocks = list_print(ma_hused, "heap", index_base);

    /* print blocks on the stack used list */
    stack_blocks = list_print(ma_sused, "stack", index_base);

    (void)printf("MA_summarize_allocated_blocks: scan completed: ");
    (void)printf("%d heap block%s, %d stack block%s\n",
        heap_blocks,
        plural(heap_blocks),
        stack_blocks,
        plural(stack_blocks));
}

/**
 ** public routines
 **/

/* ------------------------------------------------------------------------- */
/*
 * Convenience function that combines MA_allocate_heap and MA_get_index.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_alloc_get(
    Integer    datatype,    /* of elements in this block */
    Integer    nelem,        /* # of elements in this block */
    const char    *name,        /* assigned to this block by client */
    Integer    *memhandle,    /* RETURN: handle for this block */
    MA_AccessIndex    *index  /* RETURN: index for this block */   )
{
#ifdef STATS
    ma_stats.calls[(int)FID_MA_alloc_get]++;
#endif /* STATS */

    if (MA_allocate_heap(datatype, nelem, name, memhandle))
        /* MA_allocate_heap succeeded; try MA_get_index */
        return MA_get_index(*memhandle, index);
    else
        /* MA_allocate_heap failed */
        return MA_FALSE;
}

/* ------------------------------------------------------------------------- */
/*
 * Allocate a heap block big enough to hold nelem elements
 * of the given datatype.
 *
 * Return MA_TRUE upon success, or MA_FALSE upon failure.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_allocate_heap(
    Integer    datatype,    /* of elements in this block */
    Integer    nelem,        /* # of elements in this block */
    const char    *name,        /* assigned to this block by client */
    Integer    *memhandle    /* RETURN: handle for this block */ )
{
    AR        ar;        /* allocation request */
    AD        *ad;        /* AD for newly allocated block */
    Pointer    client_space;    /* location of client_space */
    ulongi    nbytes;        /* length of block for ar */
    Pointer    new_hp;        /* new ma_hp */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_allocate_heap]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return MA_FALSE;
#endif /* VERIFY */

    if (ma_trace) 
    (void)printf("MA: allocating '%s' (%d)\n", name, (int)nelem);

    /* verify initialization */
    if (!ma_initialized)
    {
        (void)sprintf(ma_ebuf,
            "block '%s', MA not yet initialized",
            name);
        ma_error(EL_Nonfatal, ET_External, "MA_allocate_heap", ma_ebuf);
        return MA_FALSE;
    }

    /* verify datatype */
    if (!mt_valid(datatype))
    {
        (void)sprintf(ma_ebuf,
            "block '%s', invalid datatype: %ld",
            name, (long)datatype);
        ma_error(EL_Nonfatal, ET_External, "MA_allocate_heap", ma_ebuf);
        return MA_FALSE;
    }

    /* verify nelem */
    if (nelem < 0)
    {
        (void)sprintf(ma_ebuf,
            "block '%s', invalid nelem: %ld",
            name, (long)nelem);
        ma_error(EL_Nonfatal, ET_External, "MA_allocate_heap", ma_ebuf);
        return MA_FALSE;
    }

    /* convert datatype to internal (index-suitable) value */
    datatype = mt_import(datatype);

    /*
     * attempt to allocate space
     */

    ar.datatype = datatype;
    ar.nelem = nelem;

    /* search the free list */
    ad = list_delete_one(&ma_hfree, ad_big_enough, (Pointer)&ar);

    /* if search of free list failed, try expanding heap region */
    if (ad == (AD *)NULL)
    {
        /* perform trial allocation to determine size */
        balloc_after(&ar, ma_hp, &client_space, &nbytes);

        new_hp = ma_hp + nbytes;
        if (new_hp > ma_sp)
        {
            (void)sprintf(ma_ebuf,
                "block '%s', not enough space to allocate %lu bytes",
                name, nbytes);
            ma_error(EL_Nonfatal, ET_External, "MA_allocate_heap", ma_ebuf);
            return MA_FALSE;
        }
        else
        {
            /* heap region expanded successfully */
            ad = (AD *)ma_hp;

            /* set fields appropriately */
            ad->client_space = client_space;
            ad->nbytes = nbytes;
        }
    }

    /*
     * space has been allocated
     */

    /* initialize the AD */
    ad->datatype = datatype;
    ad->nelem = nelem;
    str_ncopy(ad->name, (char*)name, MA_NAMESIZE);
    /* ad->client_space is already set */
    /* ad->nbytes is already set */
    list_insert(ad, &ma_hused);
    ad->checksum = checksum(ad);

    /* set the guards */
    guard_set(ad);

#ifdef DEBUG
    debug_ad_print(ad);
#endif /* DEBUG */

    /* update ma_hp if necessary */
    new_hp = (Pointer)ad + ad->nbytes;
    if (new_hp > ma_hp)
    {
        ma_hp = new_hp;
    }

#ifdef STATS
    ma_stats.hblocks++;
    ma_stats.hblocks_max = max(ma_stats.hblocks, ma_stats.hblocks_max);
    ma_stats.hbytes += ad->nbytes;
    ma_stats.hbytes_max = max(ma_stats.hbytes, ma_stats.hbytes_max);
#endif /* STATS */

    /* convert AD to memhandle */
    if ((*memhandle = ma_table_allocate((TableData)ad)) == TABLE_HANDLE_NONE)
        /* failure */
        return MA_FALSE;
    else
        /* success */
        return MA_TRUE;
}

/* ------------------------------------------------------------------------- */
/*
 * Deallocate the given stack block and all stack blocks allocated
 * after it.
 *
 * Return MA_TRUE upon success, or MA_FALSE upon failure.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_chop_stack(Integer memhandle)/*the block to deallocate up to*/
{
    AD        *ad;        /* AD for memhandle */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_chop_stack]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return MA_FALSE;
#endif /* VERIFY */

    /* verify memhandle and convert to AD */
    if (!mh2ad(memhandle, &ad, BL_Stack, "MA_chop_stack"))
        return MA_FALSE;

    /* delete block and all blocks above it from used list */
#ifdef STATS
    ma_stats.sblocks -=
        list_delete_many(&ma_sused, ad_le, (Pointer)ad, mh_free);
#else
    (void)list_delete_many(&ma_sused, ad_le, (Pointer)ad, mh_free);
#endif /* STATS */

    /* pop block and all blocks above it from stack */
#ifdef STATS
    ma_stats.sbytes -= (((Pointer)ad + ad->nbytes) - ma_sp);
#endif /* STATS */
    ma_sp = (Pointer)ad + ad->nbytes;

    /* success */
    return MA_TRUE;
}

/* ------------------------------------------------------------------------- */
/*
 * Deallocate the given heap block.
 *
 * Return MA_TRUE upon success, or MA_FALSE upon failure.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_free_heap(Integer memhandle) /* the block to deallocate */
{
    AD        *ad;        /* AD for memhandle */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_free_heap]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return MA_FALSE;
#endif /* VERIFY */

    /* verify memhandle and convert to AD */
    if (!mh2ad(memhandle, &ad, BL_Heap, "MA_free_heap"))
        return MA_FALSE;

    if (ma_trace) 
    (void)printf("MA: freeing '%s'\n", ad->name);

    /* delete block from used list */
    if (list_delete(ad, &ma_hused) != ad)
    {
        (void)sprintf(ma_ebuf,
            "memhandle %ld (name: '%s') not on heap used list",
            (long)memhandle, ad->name);
        ma_error(EL_Nonfatal, ET_Internal, "MA_free_heap", ma_ebuf);
        return MA_FALSE;
    }

#ifdef STATS
    ma_stats.hblocks--;
    ma_stats.hbytes -= ad->nbytes;
#endif /* STATS */

    /* reclaim the deallocated block */
    block_free_heap(ad);

    /* free memhandle */
    ma_table_deallocate(memhandle);

    /* success */
    return MA_TRUE;
}

/* ------------------------------------------------------------------------- */
/*
 * Deallocate nelem elements from the given heap block.
 *
 * The nelem elements (of the datatype specified when the heap block
 * was allocated) to be deallocated are assumed to be at the rightmost
 * (i.e., higher addresses) edge of the heap block.
 *
 * Return MA_TRUE upon success, or MA_FALSE upon failure.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_free_heap_piece(
    Integer    memhandle,    /* the block to deallocate a piece of */
    Integer    nelem         /* # of elements to deallocate */)
{
    AD        *ad;        /* AD for memhandle */
    AD        *ad_reclaim;    /* AD for data returned */
    AR        ar;        /* AR for data kept */
    Pointer    client_space;    /* location of client_space */
    ulongi    nbytes;        /* length of block for data kept */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_free_heap_piece]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return MA_FALSE;
#endif /* VERIFY */

    /* verify memhandle and convert to AD */
    if (!mh2ad(memhandle, &ad, BL_Heap, "MA_free_heap_piece"))
        return MA_FALSE;

    /* verify nelem */
    if (nelem < 0)
    {
        (void)sprintf(ma_ebuf,
            "block '%s', invalid nelem: %ld",
            ad->name, (long)nelem);
        ma_error(EL_Nonfatal, ET_External, "MA_free_heap_piece", ma_ebuf);
        return MA_FALSE;
    }
    else if (nelem >= ad->nelem)
    {
        /* deallocate the whole block */
        return MA_free_heap(memhandle);
    }

    if (ma_trace) 
    (void)printf("MA: freeing %ld elements of '%s'\n",
            (long)nelem, ad->name);

    /* set AR for data to keep */
    ar.datatype = ad->datatype;
    ar.nelem = ad->nelem - nelem;

    /* perform trial allocation to determine size */
    balloc_after(&ar, (Pointer)ad, &client_space, &nbytes);

    if (nbytes < ad->nbytes)
    {
        /* ad has extra space; split block if possible */
        ad_reclaim = block_split(ad, nbytes, (Boolean)MA_FALSE);

        if (ad_reclaim)
        {
#ifdef STATS
            ma_stats.hbytes -= ad_reclaim->nbytes;
#endif /* STATS */

            /* reclaim the deallocated (new) block */
            block_free_heap(ad_reclaim);
        }
    }

    /* update surviving block */
    ad->nelem = ar.nelem;
    ad->checksum = checksum(ad);

    /* set the guards */
    guard_set(ad);

#ifdef DEBUG
    debug_ad_print(ad);
#endif /* DEBUG */

    /* success */
    return MA_TRUE;
}

/* ------------------------------------------------------------------------- */
/*
 * Get the base index for the given block.
 *
 * Return MA_TRUE upon success, or MA_FALSE upon failure.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_get_index(
    Integer    memhandle,    /* block to get index for */
    MA_AccessIndex    *index         /* RETURN: base index */)
{
    AD        *ad;        /* AD for memhandle */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_get_index]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return MA_FALSE;
#endif /* VERIFY */

    /* verify memhandle and convert to AD */
    if (mh2ad(memhandle, &ad, BL_HeapOrStack, "MA_get_index"))
    {
        /* compute index */
        *index = client_space_index(ad);

        /* success */
        return MA_TRUE;
    }
    else
    {
        /* failure */
        return MA_FALSE;
    }
}

/* ------------------------------------------------------------------------- */
/*
 * Return the base address of the given datatype.
 */
/* ------------------------------------------------------------------------- */

public Pointer MA_get_mbase(Integer datatype)    /* to get base address of */
{
#ifdef STATS
    ma_stats.calls[(int)FID_MA_get_mbase]++;
#endif /* STATS */

    /* preinitialize if necessary */
    ma_preinitialize("MA_get_mbase");

    /* verify datatype */
    if (!mt_valid(datatype))
    {
        (void)sprintf(ma_ebuf,
            "invalid datatype: %ld",
            (long)datatype);
        ma_error(EL_Fatal, ET_External, "MA_get_mbase", ma_ebuf);
        return NULL;
    }

    /* convert datatype to internal (index-suitable) value */
    datatype = mt_import(datatype);

    return ma_base[datatype];
}

/* ------------------------------------------------------------------------- */
/*
 * Get the handle for the next block in the scan of currently allocated
 * blocks.
 *
 * Return MA_TRUE upon success, or MA_FALSE upon failure.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_get_next_memhandle(
    Integer    *ithandle,    /* handle for this iterator */
    Integer    *memhandle     /* RETURN: handle for the next block */)
{
#ifdef STATS
    ma_stats.calls[(int)FID_MA_get_next_memhandle]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return MA_FALSE;
#endif /* VERIFY */

    /* not yet implemented */
    (void)sprintf(ma_ebuf,
        "not yet implemented");
    ma_error(EL_Nonfatal, ET_External, "MA_get_next_memhandle", ma_ebuf);
    return MA_FALSE;
}

/* ------------------------------------------------------------------------- */
/*
 * Get the requested alignment.
 *
 * Return MA_TRUE upon success, or MA_FALSE upon failure.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_get_numalign(Integer *value)
    /* RETURN: requested alignment */
{
#ifdef STATS
    ma_stats.calls[(int)FID_MA_get_numalign]++;
#endif /* STATS */

    *value = ma_numalign;
    return MA_TRUE;
}

/* ------------------------------------------------------------------------- */
/*
 * Get the base pointer for the given block.
 *
 * Return MA_TRUE upon success, or MA_FALSE upon failure.
 */
/* ------------------------------------------------------------------------- */

/* JN  converted to void* to avoid calling hassles */
public Boolean MA_get_pointer(
    Integer    memhandle,    /* block to get pointer for */
    void    *pointer     /* RETURN: base pointer */)
{
    AD        *ad;        /* AD for memhandle */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_get_pointer]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return MA_FALSE;
#endif /* VERIFY */

    /* verify memhandle and convert to AD */
    if (mh2ad(memhandle, &ad, BL_HeapOrStack, "MA_get_pointer"))
    {
        /* compute pointer */
#if 0
        *pointer = ad->client_space;
#endif
        *(char**)pointer = ad->client_space;

        /* success */
        return MA_TRUE;
    }
    else
    {
        /* failure */
        return MA_FALSE;
    }
}

/* ------------------------------------------------------------------------- */
/*
 * Initialize the memory allocator.
 *
 * Return MA_TRUE upon success, or MA_FALSE upon failure.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_init(
    Integer    datatype,    /* for computing storage requirement */
    Integer    nominal_stack,    /* # of datatype elements desired for stack */
    Integer    nominal_heap     /* # of datatype elements desired for heap */)
{
    ulongi    heap_bytes;    /* # of bytes for heap */
    ulongi    stack_bytes;    /* # of bytes for stack */
    ulongi    total_bytes;    /* total # of bytes */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_init]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return MA_FALSE;
#endif /* VERIFY */

    /* preinitialize if necessary */
    ma_preinitialize("MA_init");

    /* verify uninitialization */
    if (ma_initialized)
    {
        (void)sprintf(ma_ebuf,
            "MA already initialized");
        ma_error(EL_Nonfatal, ET_External, "MA_init", ma_ebuf);
        return MA_FALSE;
    }

    /* verify datatype */
    if (!mt_valid(datatype))
    {
        (void)sprintf(ma_ebuf,
            "invalid datatype: %ld",
            (long)datatype);
        ma_error(EL_Nonfatal, ET_External, "MA_init", ma_ebuf);
        return MA_FALSE;
    }

    /* convert datatype to internal (index-suitable) value */
    datatype = mt_import(datatype);

    /* compute # of bytes in heap */
    if (nominal_heap < 0)
    {
        heap_bytes = DEFAULT_TOTAL_HEAP;
    }
    else
    {
        heap_bytes = (nominal_heap * ma_sizeof[datatype]) +
            (DEFAULT_REQUESTS_HEAP * max_block_overhead(datatype));
    }
    heap_bytes = (unsigned long)mai_round((long)heap_bytes, (ulongi)ALIGNMENT);

    /* compute # of bytes in stack */
    if (nominal_stack < 0)
    {
        stack_bytes = DEFAULT_TOTAL_STACK;
    }
    else
    {
        stack_bytes = (nominal_stack * ma_sizeof[datatype]) +
            (DEFAULT_REQUESTS_STACK * max_block_overhead(datatype));
    }
    stack_bytes = (unsigned long)mai_round((long)stack_bytes, (ulongi)ALIGNMENT);

    /* segment consists of heap and stack */
    total_bytes = heap_bytes + stack_bytes;
#ifdef NOUSE_MMAP
    /* disable memory mapped malloc */
    mallopt(M_MMAP_MAX, 0);
    mallopt(M_TRIM_THRESHOLD, -1);
#endif
    /* allocate the segment of memory */
#ifdef ENABLE_ARMCI_MEM_OPTION
    if(getenv("MA_USE_ARMCI_MEM"))
    {
        ma_segment = (Pointer)ARMCI_Malloc_local(total_bytes);
    }
    else
#endif
        ma_segment = (Pointer)bytealloc(total_bytes);
    if (ma_segment == (Pointer)NULL)
    {
        (void)sprintf(ma_ebuf,
            "could not allocate %lu bytes",
            total_bytes);
        ma_error(EL_Nonfatal, ET_External, "MA_init", ma_ebuf);
        return MA_FALSE;
    }

    /*
     * initialize management stuff
     */

    /* partition is at (first address past) end of heap */
    ma_partition = ma_segment + heap_bytes;
    /* compute (first address past) end of segment */
    ma_eos = ma_segment + total_bytes;
    /* ma_hp initially points at start of segment */
    ma_hp = ma_segment;
    /* ma_sp initially points at end of segment */
    ma_sp = ma_eos;

    /* lists are all initially empty */
    ma_hfree = (AD *)NULL;
    ma_hused = (AD *)NULL;
    ma_sused = (AD *)NULL;

    /* we are now initialized */
    ma_initialized = MA_TRUE;

    /* success */
    return MA_TRUE;
}

/* ------------------------------------------------------------------------- */
/*
 * Return MA_TRUE if MA_init has been called successfully,
 * else return MA_FALSE.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_initialized()
{
#ifdef STATS
    ma_stats.calls[(int)FID_MA_initialized]++;
#endif /* STATS */

    return ma_initialized;
}

/* ------------------------------------------------------------------------- */
/*
 * Initialize a scan of currently allocated blocks.
 *
 * Return MA_TRUE upon success, or MA_FALSE upon failure.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_init_memhandle_iterator( Integer *ithandle)
{
#ifdef STATS
    ma_stats.calls[(int)FID_MA_init_memhandle_iterator]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return MA_FALSE;
#endif /* VERIFY */

    /* not yet implemented */
    (void)sprintf(ma_ebuf,
        "not yet implemented");
    ma_error(EL_Nonfatal, ET_External, "MA_init_memhandle_iterator", ma_ebuf);
    return MA_FALSE;
}

/* ------------------------------------------------------------------------- */
/*
 * Return the maximum number of datatype elements that can currently
 * be allocated in the space between the heap and the stack, in a single
 * allocation request, ignoring the partition defined at initialization.
 *
 * Note that this might not be the largest piece of memory available;
 * the heap may contain deallocated blocks that are larger.
 */
/* ------------------------------------------------------------------------- */

public Integer MA_inquire_avail(Integer datatype)
{
    long    gap_length;    /* # of bytes between heap and stack */
    Integer    nelem_gap;    /* max elements containable in gap */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_inquire_avail]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return DONTCARE;
#endif /* VERIFY */

    /* verify initialization */
    if (!ma_initialized)
    {
        (void)sprintf(ma_ebuf,
            "MA not yet initialized");
        ma_error(EL_Nonfatal, ET_External, "MA_inquire_avail", ma_ebuf);
        return (Integer)0;
    }

    /* verify datatype */
    if (!mt_valid(datatype))
    {
        (void)sprintf(ma_ebuf,
            "invalid datatype: %ld",
            (long)datatype);
        ma_error(EL_Fatal, ET_External, "MA_inquire_avail", ma_ebuf);
        return DONTCARE;
    }

    /* convert datatype to internal (index-suitable) value */
    datatype = mt_import(datatype);

    /*
     * compute the # of elements for which space is available
     */

    /* try space between heap and stack */
    gap_length = (long)(ma_sp - ma_hp);
    if (gap_length > 0)
        nelem_gap = ma_nelem(ma_hp, (ulongi)gap_length, datatype);
    else
        nelem_gap = 0;

    /* success */
    return nelem_gap;
}

/* ------------------------------------------------------------------------- */
/*
 * Return the maximum number of datatype elements that can currently
 * be allocated in the heap, in a single allocation request,
 * honoring the partition defined at initialization.
 *
 * This routine does not check the stack.  Therefore, if the stack
 * has overgrown the partition, the answer returned by this routine
 * might be incorrect, i.e., there might be less memory available
 * than this routine indicates.
 */
/* ------------------------------------------------------------------------- */

public Integer MA_inquire_heap(Integer datatype)
{
    long    gap_length;    /* # of bytes between heap and partition */
    Integer    nelem_gap;    /* max elements containable in gap */
    Integer    nelem_frag;    /* max elements containable in any frag */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_inquire_heap]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return DONTCARE;
#endif /* VERIFY */

    /* verify initialization */
    if (!ma_initialized)
    {
        (void)sprintf(ma_ebuf,
            "MA not yet initialized");
        ma_error(EL_Nonfatal, ET_External, "MA_inquire_heap", ma_ebuf);
        return (Integer)0;
    }

    /* verify datatype */
    if (!mt_valid(datatype))
    {
        (void)sprintf(ma_ebuf,
            "invalid datatype: %ld",
            (long)datatype);
        ma_error(EL_Fatal, ET_External, "MA_inquire_heap", ma_ebuf);
        return DONTCARE;
    }

    /* convert datatype to internal (index-suitable) value */
    datatype = mt_import(datatype);

    /*
     * compute the # of elements for which space is available
     */

    /* try space between heap and partition */
    gap_length = (long)(ma_partition - ma_hp);
    if (gap_length > 0)
        nelem_gap = ma_nelem(ma_hp, (ulongi)gap_length, datatype);
    else
        nelem_gap = 0;

    /* try heap fragments */
    nelem_frag = ma_max_heap_frag_nelem(datatype, nelem_gap);

    /* success */
    return max(nelem_gap, nelem_frag);
}

/* ------------------------------------------------------------------------- */
/*
 * Return the maximum number of datatype elements that can currently
 * be allocated in the heap, in a single allocation request,
 * honoring the partition defined at initialization.
 *
 * This routine does check the stack.  Therefore, whether or not the stack
 * has overgrown the partition, the answer returned by this routine
 * will be correct, i.e., there will be at least the memory available
 * that this routine indicates.
 *
 * Note that this might not be the largest piece of memory available;
 * the space between the heap and the stack may be larger.
 */
/* ------------------------------------------------------------------------- */

public Integer MA_inquire_heap_check_stack(Integer datatype)
{
    long    gap_length;    /* # of bytes between heap and partition */
    Integer    nelem_gap;    /* max elements containable in gap */
    Integer    nelem_frag;    /* max elements containable in any frag */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_inquire_heap_check_stack]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return DONTCARE;
#endif /* VERIFY */

    /* verify initialization */
    if (!ma_initialized)
    {
        (void)sprintf(ma_ebuf,
            "MA not yet initialized");
        ma_error(EL_Nonfatal, ET_External, "MA_inquire_heap_check_stack", ma_ebuf);
        return (Integer)0;
    }

    /* verify datatype */
    if (!mt_valid(datatype))
    {
        (void)sprintf(ma_ebuf,
            "invalid datatype: %ld",
            (long)datatype);
        ma_error(EL_Fatal, ET_External, "MA_inquire_heap_check_stack", ma_ebuf);
        return DONTCARE;
    }

    /* convert datatype to internal (index-suitable) value */
    datatype = mt_import(datatype);

    /*
     * compute the # of elements for which space is available
     */

    /* try space between heap and partition or heap and stack */
    gap_length = min((long)(ma_partition - ma_hp), (long)(ma_sp - ma_hp));
    if (gap_length > 0)
        nelem_gap = ma_nelem(ma_hp, (ulongi)gap_length, datatype);
    else
        nelem_gap = 0;

    /* try heap fragments */
    nelem_frag = ma_max_heap_frag_nelem(datatype, nelem_gap);

    /* success */
    return max(nelem_gap, nelem_frag);
}

/* ------------------------------------------------------------------------- */
/*
 * Return the maximum number of datatype elements that can currently
 * be allocated in the heap, in a single allocation request,
 * ignoring the partition defined at initialization.
 *
 * This routine does check the stack.  Therefore, whether or not the stack
 * has overgrown the partition, the answer returned by this routine
 * will be correct, i.e., there will be at least the memory available
 * that this routine indicates.
 *
 * Note that this will be the largest piece of memory available.
 */
/* ------------------------------------------------------------------------- */

public Integer MA_inquire_heap_no_partition(Integer datatype)
{
    long    gap_length;    /* # of bytes between heap and partition */
    Integer    nelem_gap;    /* max elements containable in gap */
    Integer    nelem_frag;    /* max elements containable in any frag */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_inquire_heap_no_partition]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return DONTCARE;
#endif /* VERIFY */

    /* verify initialization */
    if (!ma_initialized)
    {
        (void)sprintf(ma_ebuf,
            "MA not yet initialized");
        ma_error(EL_Nonfatal, ET_External, "MA_inquire_heap_no_partition", ma_ebuf);
        return (Integer)0;
    }

    /* verify datatype */
    if (!mt_valid(datatype))
    {
        (void)sprintf(ma_ebuf,
            "invalid datatype: %ld",
            (long)datatype);
        ma_error(EL_Fatal, ET_External, "MA_inquire_heap_no_partition", ma_ebuf);
        return DONTCARE;
    }

    /* convert datatype to internal (index-suitable) value */
    datatype = mt_import(datatype);

    /*
     * compute the # of elements for which space is available
     */

    /* try space between heap and stack */
    gap_length = (long)(ma_sp - ma_hp);
    if (gap_length > 0)
        nelem_gap = ma_nelem(ma_hp, (ulongi)gap_length, datatype);
    else
        nelem_gap = 0;

    /* try heap fragments */
    nelem_frag = ma_max_heap_frag_nelem(datatype, nelem_gap);

    /* success */
    return max(nelem_gap, nelem_frag);
}

/* ------------------------------------------------------------------------- */
/*
 * Return the maximum number of datatype elements that can currently
 * be allocated in the stack, in a single allocation request,
 * honoring the partition defined at initialization.
 *
 * This routine does not check the heap.  Therefore, if the heap
 * has overgrown the partition, the answer returned by this routine
 * might be incorrect, i.e., there might be less memory available
 * than this routine indicates.
 */
/* ------------------------------------------------------------------------- */

public Integer MA_inquire_stack(Integer datatype)
{
    long    gap_length;    /* # of bytes between partition and stack */
    Integer    nelem_gap;    /* max elements containable in gap */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_inquire_stack]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return DONTCARE;
#endif /* VERIFY */

    /* verify initialization */
    if (!ma_initialized)
    {
        (void)sprintf(ma_ebuf,
            "MA not yet initialized");
        ma_error(EL_Nonfatal, ET_External, "MA_inquire_stack", ma_ebuf);
        return (Integer)0;
    }

    /* verify datatype */
    if (!mt_valid(datatype))
    {
        (void)sprintf(ma_ebuf,
            "invalid datatype: %ld",
            (long)datatype);
        ma_error(EL_Fatal, ET_External, "MA_inquire_stack", ma_ebuf);
        return DONTCARE;
    }

    /* convert datatype to internal (index-suitable) value */
    datatype = mt_import(datatype);

    /*
     * compute the # of elements for which space is available
     */

    /* try space between partition and stack */
    gap_length = (long)(ma_sp - ma_partition);
    if (gap_length > 0)
        nelem_gap = ma_nelem(ma_partition, (ulongi)gap_length, datatype);
    else
        nelem_gap = 0;

    /* success */
    return nelem_gap;
}

/* ------------------------------------------------------------------------- */
/*
 * Return the maximum number of datatype elements that can currently
 * be allocated in the stack, in a single allocation request,
 * honoring the partition defined at initialization.
 *
 * This routine does check the heap.  Therefore, whether or not the heap
 * has overgrown the partition, the answer returned by this routine
 * will be correct, i.e., there will be at least the memory available
 * that this routine indicates.
 *
 * Note that this might not be the largest piece of memory available;
 * the space between the heap and the stack may be larger.
 */
/* ------------------------------------------------------------------------- */

public Integer MA_inquire_stack_check_heap(Integer datatype)
{
    long    gap_length;    /* # of bytes between partition and stack */
    Integer    nelem_gap;    /* max elements containable in gap */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_inquire_stack_check_heap]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return DONTCARE;
#endif /* VERIFY */

    /* verify initialization */
    if (!ma_initialized)
    {
        (void)sprintf(ma_ebuf,
            "MA not yet initialized");
        ma_error(EL_Nonfatal, ET_External, "MA_inquire_stack_check_heap", ma_ebuf);
        return (Integer)0;
    }

    /* verify datatype */
    if (!mt_valid(datatype))
    {
        (void)sprintf(ma_ebuf,
            "invalid datatype: %ld",
            (long)datatype);
        ma_error(EL_Fatal, ET_External, "MA_inquire_stack_check_heap", ma_ebuf);
        return DONTCARE;
    }

    /* convert datatype to internal (index-suitable) value */
    datatype = mt_import(datatype);

    /*
     * compute the # of elements for which space is available
     */

    /* try space between partition and stack or heap and stack */
    gap_length = min((long)(ma_sp - ma_partition), (long)(ma_sp - ma_hp));
    if (gap_length > 0)
        nelem_gap = ma_nelem(ma_partition, (ulongi)gap_length, datatype);
    else
        nelem_gap = 0;

    /* success */
    return nelem_gap;
}

/* ------------------------------------------------------------------------- */
/*
 * Return the maximum number of datatype elements that can currently
 * be allocated in the stack, in a single allocation request,
 * ignoring the partition defined at initialization.
 *
 * This routine does check the heap.  Therefore, whether or not the heap
 * has overgrown the partition, the answer returned by this routine
 * will be correct, i.e., there will be at least the memory available
 * that this routine indicates.
 *
 * Note that this might not be the largest piece of memory available;
 * the heap may contain deallocated blocks that are larger.
 *
 * This routine is equivalent to MA_inquire_avail.
 */
/* ------------------------------------------------------------------------- */

public Integer MA_inquire_stack_no_partition(Integer datatype)
{
    long    gap_length;    /* # of bytes between heap and partition */
    Integer    nelem_gap;    /* max elements containable in gap */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_inquire_stack_no_partition]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return DONTCARE;
#endif /* VERIFY */

    /* verify initialization */
    if (!ma_initialized)
    {
        (void)sprintf(ma_ebuf,
            "MA not yet initialized");
        ma_error(EL_Nonfatal, ET_External, "MA_inquire_stack_no_partition", ma_ebuf);
        return (Integer)0;
    }

    /* verify datatype */
    if (!mt_valid(datatype))
    {
        (void)sprintf(ma_ebuf,
            "invalid datatype: %ld",
            (long)datatype);
        ma_error(EL_Fatal, ET_External, "MA_inquire_stack_no_partition", ma_ebuf);
        return DONTCARE;
    }

    /* convert datatype to internal (index-suitable) value */
    datatype = mt_import(datatype);

    /*
     * compute the # of elements for which space is available
     */

    /* try space between heap and stack */
    gap_length = (long)(ma_sp - ma_hp);
    if (gap_length > 0)
        nelem_gap = ma_nelem(ma_hp, (ulongi)gap_length, datatype);
    else
        nelem_gap = 0;

    /* success */
    return nelem_gap;
}

/* ------------------------------------------------------------------------- */
/*
 * Deallocate the given stack block, which must be the one most recently
 * allocated.
 *
 * Return MA_TRUE upon success, or MA_FALSE upon failure.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_pop_stack(Integer memhandle) /* the block to deallocate */
{
    AD        *ad;        /* AD for memhandle */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_pop_stack]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return MA_FALSE;
#endif /* VERIFY */

    /* verify memhandle and convert to AD */
    if (!mh2ad(memhandle, &ad, BL_StackTop, "MA_pop_stack"))
        return MA_FALSE;

    if (ma_trace) 
    (void)printf("MA: popping '%s'\n", ad->name);

    /* delete block from used list */
    if (list_delete(ad, &ma_sused) != ad)
    {
        (void)sprintf(ma_ebuf,
            "memhandle %ld (name: '%s') not on stack used list",
            (long)memhandle, ad->name);
        ma_error(EL_Nonfatal, ET_Internal, "MA_pop_stack", ma_ebuf);
        return MA_FALSE;
    }

    /* pop block from stack */
    ma_sp += ad->nbytes;

#ifdef STATS
    ma_stats.sblocks--;
    ma_stats.sbytes -= ad->nbytes;
#endif /* STATS */

    /* free memhandle */
    ma_table_deallocate(memhandle);

    /* success */
    return MA_TRUE;
}

/* ------------------------------------------------------------------------- */
/*
 * Print usage statistics.
 */
/* ------------------------------------------------------------------------- */

public void MA_print_stats(Boolean printroutines)
{
#ifdef STATS

    int        i;        /* loop index */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_print_stats]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return;
#endif /* VERIFY */

    (void)printf("MA usage statistics:\n");
    (void)printf("\n\tallocation statistics:\n");
    (void)printf("\t\t\t\t\t      heap\t     stack\n");
    (void)printf("\t\t\t\t\t      ----\t     -----\n");
    (void)printf("\tcurrent number of blocks\t%10lu\t%10lu\n",
        ma_stats.hblocks,
        ma_stats.sblocks);
    (void)printf("\tmaximum number of blocks\t%10lu\t%10lu\n",
        ma_stats.hblocks_max,
        ma_stats.sblocks_max);
    (void)printf("\tcurrent total bytes\t\t%10lu\t%10lu\n",
        ma_stats.hbytes,
        ma_stats.sbytes);
    (void)printf("\tmaximum total bytes\t\t%10lu\t%10lu\n",
        ma_stats.hbytes_max,
        ma_stats.sbytes_max);
    (void)printf("\tmaximum total K-bytes\t\t%10lu\t%10lu\n",
        ((ma_stats.hbytes_max+999)/1000),
        ((ma_stats.sbytes_max+999)/1000));
    (void)printf("\tmaximum total M-bytes\t\t%10lu\t%10lu\n",
        ((ma_stats.hbytes_max+999999)/1000000),
        ((ma_stats.sbytes_max+999999)/1000000));
    if (printroutines)
    {
        (void)printf("\n\tcalls per routine:\n");
        for (i = 0; i < NUMROUTINES; i++)
            (void)printf("\t\t%10lu  %s\n",
                ma_stats.calls[i],
                ma_routines[i]);
    }

#else

    (void)sprintf(ma_ebuf,
        "unavailable; recompile MA with -DSTATS");
    ma_error(EL_Nonfatal, ET_External, "MA_print_stats", ma_ebuf);

#endif /* STATS */
}

/* ------------------------------------------------------------------------- */
/*
 * Convenience function that combines MA_push_stack and MA_get_index.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_push_get(
    Integer    datatype,    /* of elements in this block */
    Integer    nelem,        /* # of elements in this block */
    const char    *name,        /* assigned to this block by client */
    Integer    *memhandle,    /* RETURN: handle for this block */
    MA_AccessIndex  *index    /* RETURN: index for this block */)
{
#ifdef STATS
    ma_stats.calls[(int)FID_MA_push_get]++;
#endif /* STATS */

    if (MA_push_stack(datatype, nelem, name, memhandle))
        /* MA_push_stack succeeded; try MA_get_index */
        return MA_get_index(*memhandle, index);
    else
        /* MA_push_stack failed */
        return MA_FALSE;
}

/* ------------------------------------------------------------------------- */
/*
 * Allocate a stack block big enough to hold nelem elements
 * of the given datatype.
 *
 * Return MA_TRUE upon success, or MA_FALSE upon failure.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_push_stack(
    Integer    datatype,    /* of elements in this block */
    Integer    nelem,        /* # of elements in this block */
    const char    *name,        /* assigned to this block by client */
    Integer    *memhandle     /* RETURN: handle for this block */)
{
    AR        ar;        /* allocation request */
    AD        *ad;        /* AD for newly allocated block */
    Pointer    client_space;    /* location of client_space */
    ulongi    nbytes;        /* length of block for ar */
    Pointer    new_sp;        /* new ma_sp */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_push_stack]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return MA_FALSE;
#endif /* VERIFY */

    if (ma_trace) 
    (void)printf("MA: pushing '%s' (%d)\n", name, (int)nelem);

    /* verify initialization */
    if (!ma_initialized)
    {
        (void)sprintf(ma_ebuf,
            "block '%s', MA not yet initialized",
            name);
        ma_error(EL_Nonfatal, ET_External, "MA_push_stack", ma_ebuf);
        return MA_FALSE;
    }

    /* verify datatype */
    if (!mt_valid(datatype))
    {
        (void)sprintf(ma_ebuf,
            "block '%s', invalid datatype: %ld",
            name, (long)datatype);
        ma_error(EL_Nonfatal, ET_External, "MA_push_stack", ma_ebuf);
        return MA_FALSE;
    }

    /* verify nelem */
    if (nelem < 0)
    {
        (void)sprintf(ma_ebuf,
            "block '%s', invalid nelem: %ld",
            name, (long)nelem);
        ma_error(EL_Nonfatal, ET_External, "MA_push_stack", ma_ebuf);
        return MA_FALSE;
    }

    /* convert datatype to internal (index-suitable) value */
    datatype = mt_import(datatype);

    /*
     * attempt to allocate space
     */

    ar.datatype = datatype;
    ar.nelem = nelem;

    balloc_before(&ar, ma_sp, &client_space, &nbytes);

    new_sp = ma_sp - nbytes;
    /* if (new_sp < ma_hp) */
    if (((ulongi)(ma_sp - ma_hp)) < nbytes)
    {
        (void)sprintf(ma_ebuf,
            "block '%s', not enough space to allocate %lu bytes",
            name, nbytes);
        ma_error(EL_Nonfatal, ET_External, "MA_push_stack", ma_ebuf);
        return MA_FALSE;
    }
    else
    {
        ad = (AD *)new_sp;
    }

    /*
     * space has been allocated
     */

    /* initialize the AD */
    ad->datatype = datatype;
    ad->nelem = nelem;
    str_ncopy(ad->name, (char*)name, MA_NAMESIZE);
    ad->client_space = client_space;
    ad->nbytes = nbytes;
    list_insert(ad, &ma_sused);
    ad->checksum = checksum(ad);

    /* set the guards */
    guard_set(ad);

#ifdef DEBUG
    debug_ad_print(ad);
#endif /* DEBUG */

    /* update ma_sp */
    ma_sp = new_sp;

#ifdef STATS
    ma_stats.sblocks++;
    ma_stats.sblocks_max = max(ma_stats.sblocks, ma_stats.sblocks_max);
    ma_stats.sbytes += ad->nbytes;
    ma_stats.sbytes_max = max(ma_stats.sbytes, ma_stats.sbytes_max);
#endif /* STATS */

    /* convert AD to memhandle */
    if ((*memhandle = ma_table_allocate((TableData)ad)) == TABLE_HANDLE_NONE)
        /* failure */
        return MA_FALSE;
    else
        /* success */
        return MA_TRUE;
}

/* ------------------------------------------------------------------------- */
/*
 * Set the ma_auto_verify flag to value and return its previous value.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_set_auto_verify(Boolean  value /* to set flag to */)
{
    Boolean    old_value;    /* of flag */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_set_auto_verify]++;
#endif /* STATS */

    old_value = ma_auto_verify;
    ma_auto_verify = value;
    return old_value;
}

/* ------------------------------------------------------------------------- */
/*
 * Set the ma_error_print flag to value and return its previous value.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_set_error_print(Boolean value /* to set flag to */)
{
    Boolean    old_value;    /* of flag */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_set_error_print]++;
#endif /* STATS */

    old_value = ma_error_print;
    ma_error_print = value;
    return old_value;
}

/* ------------------------------------------------------------------------- */
/*
 * Set the ma_hard_fail flag to value and return its previous value.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_set_hard_fail( Boolean value /* to set flag to */)
{
    Boolean    old_value;    /* of flag */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_set_hard_fail]++;
#endif /* STATS */

    old_value = ma_hard_fail;
    ma_hard_fail = value;
    return old_value;
}

/* ------------------------------------------------------------------------- */
/*
 * Set the requested alignment.
 *
 * Return MA_TRUE upon success, or MA_FALSE upon failure.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_set_numalign(Integer  value)
{
#ifdef STATS
    ma_stats.calls[(int)FID_MA_set_numalign]++;
#endif /* STATS */

    if ((value < 0) || (value > 30))
    {
        (void)sprintf(ma_ebuf,
            "invalid alignment: %ld",
            (long)value);
        ma_error(EL_Nonfatal, ET_External, "MA_set_numalign", ma_ebuf);
        return MA_FALSE;
    }
    ma_numalign = value;
    return MA_TRUE;
}

/* ------------------------------------------------------------------------- */
/*
 * Return the number of elements of datatype2 required to contain
 * nelem1 elements of datatype1.
 */
/* ------------------------------------------------------------------------- */

public Integer MA_sizeof(
    Integer    datatype1,    /* of source elements */
    Integer    nelem1,        /* # of source elements */
    Integer    datatype2     /* of target elements */)
{
    ulongi    source_bytes;    /* nelem1 * ma_sizeof[datatype1] */
    int        ceiling;    /* 1 iff ceiling alters result */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_sizeof]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return DONTCARE;
#endif /* VERIFY */

    /* preinitialize if necessary */
    ma_preinitialize("MA_sizeof");

    /* verify datatype1 */
    if (!mt_valid(datatype1))
    {
        (void)sprintf(ma_ebuf,
            "invalid datatype: %ld",
            (long)datatype1);
        ma_error(EL_Fatal, ET_External, "MA_sizeof", ma_ebuf);
        return DONTCARE;
    }

    /* verify nelem1 */
    if (nelem1 < 0)
    {
        (void)sprintf(ma_ebuf,
            "invalid nelem: %ld",
            (long)nelem1);
        ma_error(EL_Fatal, ET_External, "MA_sizeof", ma_ebuf);
        return DONTCARE;
    }

    /* verify datatype2 */
    if (!mt_valid(datatype2))
    {
        (void)sprintf(ma_ebuf,
            "invalid datatype: %ld",
            (long)datatype2);
        ma_error(EL_Fatal, ET_External, "MA_sizeof", ma_ebuf);
        return DONTCARE;
    }

    /* convert datatype1 to internal (index-suitable) value */
    datatype1 = mt_import(datatype1);

    /* convert datatype2 to internal (index-suitable) value */
    datatype2 = mt_import(datatype2);

    /* compute and return the result */
    source_bytes = nelem1 * ma_sizeof[datatype1];
    ceiling = (source_bytes % ma_sizeof[datatype2]) ? 1 : 0;
    return (Integer)((source_bytes / ma_sizeof[datatype2]) + ceiling);
}

/* ------------------------------------------------------------------------- */
/*
 * Return the number of elements of datatype required to contain
 * the worst case number of bytes of overhead for any block.
 */
/* ------------------------------------------------------------------------- */

public Integer MA_sizeof_overhead(Integer datatype) 
{
    int        overhead_bytes;    /* max bytes of overhead for any block */
    int        ceiling;    /* 1 iff ceiling alters result */
    int        max_sizeof;    /* max over i of ma_sizeof[i] */
    int        biggest_datatype=0; /* corresponds to max_sizeof */
    int        i;        /* loop index */

#ifdef STATS
    ma_stats.calls[(int)FID_MA_sizeof_overhead]++;
#endif /* STATS */

#ifdef VERIFY
    if (ma_auto_verify && !MA_verify_allocator_stuff())
        return DONTCARE;
#endif /* VERIFY */

    /* preinitialize if necessary */
    ma_preinitialize("MA_sizeof_overhead");

    /* verify datatype */
    if (!mt_valid(datatype))
    {
        (void)sprintf(ma_ebuf,
            "invalid datatype: %ld",
            (long)datatype);
        ma_error(EL_Fatal, ET_External, "MA_sizeof_overhead", ma_ebuf);
        return DONTCARE;
    }

    /* convert datatype to internal (index-suitable) value */
    datatype = mt_import(datatype);

    /* compute and return the result */
    for (max_sizeof = 0, i = 0; i < MT_NUMTYPES; i++)
        if (ma_sizeof[i] > max_sizeof)
        {
            max_sizeof = ma_sizeof[i];
            biggest_datatype = i;
        }
    overhead_bytes = max_block_overhead(biggest_datatype);
    ceiling = (overhead_bytes % ma_sizeof[datatype]) ? 1 : 0;
    return (Integer)((overhead_bytes / ma_sizeof[datatype]) + ceiling);
}

/* ------------------------------------------------------------------------- */
/*
 * Print debugging information about all blocks currently in use
 * on the heap or the stack.
 */
/* ------------------------------------------------------------------------- */

public void MA_summarize_allocated_blocks()
{
    /* C indices are 0-based */
    MAi_summarize_allocated_blocks(0);
}

/* ------------------------------------------------------------------------- */
/*
 * Control tracing of allocation and deallocation operations.
 */
/* ------------------------------------------------------------------------- */

public void MA_trace(Boolean value)
{
    ma_trace = value;
}

/* ------------------------------------------------------------------------- */
/*
 * Sanity check the internal state of MA and print the results.
 *
 * Return MA_TRUE upon success, or MA_FALSE upon failure.
 */
/* ------------------------------------------------------------------------- */

public Boolean MA_verify_allocator_stuff()
{
#ifdef VERIFY

    char    *preamble;    /* printed before block error messages */

    int        heap_blocks;
    int        bad_heap_blocks;
    int        bad_heap_checksums;
    int        bad_heap_lguards;
    int        bad_heap_rguards;
    int        stack_blocks;
    int        bad_stack_blocks;
    int        bad_stack_checksums;
    int        bad_stack_lguards;
    int        bad_stack_rguards;

#ifdef STATS
    ma_stats.calls[(int)FID_MA_verify_allocator_stuff]++;
#endif /* STATS */

    preamble = "MA_verify_allocator_stuff: starting scan ...\n";

    /* check each block on the heap used list */
    list_verify(ma_hused,
        "heap",
        preamble,
        &heap_blocks,
        &bad_heap_blocks,
        &bad_heap_checksums,
        &bad_heap_lguards,
        &bad_heap_rguards);

    if (bad_heap_blocks > 0)
        /* only print preamble once */
        preamble = (char *)NULL;

    /* check each block on the stack used list */
    list_verify(ma_sused,
        "stack",
        preamble,
        &stack_blocks,
        &bad_stack_blocks,
        &bad_stack_checksums,
        &bad_stack_lguards,
        &bad_stack_rguards);

    if ((bad_heap_blocks > 0) || (bad_stack_blocks > 0))
    {
        Boolean    old_ma_error_print;

        /* print postamble */
        (void)printf("MA_verify_allocator_stuff: scan completed\n");

        /* construct a summary of the results */
        (void)sprintf(ma_ebuf, "\n\t\t\t\theap\tstack\n\t\t\t\t----\t-----\n\tchecksum errors\t\t%4d\t%5d\n\tleft signature errors\t%4d\t%5d\n\tright signature errors\t%4d\t%5d\n\ttotal bad blocks\t%4d\t%5d\n\ttotal blocks\t\t%4d\t%5d",
            bad_heap_checksums,
            bad_stack_checksums,
            bad_heap_lguards,
            bad_stack_lguards,
            bad_heap_rguards,
            bad_stack_rguards,
            bad_heap_blocks,
            bad_stack_blocks,
            heap_blocks,
            stack_blocks);

        /* print the summary on stderr */
        old_ma_error_print = ma_error_print;
        ma_error_print = MA_TRUE;
        ma_error(EL_Nonfatal, ET_External, "MA_verify_allocator_stuff", ma_ebuf);
        ma_error_print = old_ma_error_print;

        /* problems were found */
        return MA_FALSE;
    }
    else
        /* no problems found */
        return MA_TRUE;

#else

#ifdef STATS
    ma_stats.calls[(int)FID_MA_verify_allocator_stuff]++;
#endif /* STATS */

    (void)sprintf(ma_ebuf,
        "unavailable; recompile MA with -DVERIFY");
    ma_error(EL_Nonfatal, ET_External, "MA_verify_allocator_stuff", ma_ebuf);
    return MA_FALSE;

#endif /* VERIFY */
}
