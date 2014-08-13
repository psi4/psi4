/** @file
 * $ID:$
 */
#ifndef RA_COMMON_H_
#define RA_COMMON_H_

#if SIZEOF_LONG == 8
#   define POLY 0x0000000000000007UL
#   define PERIOD 1317624576693539401L
    typedef unsigned long u64Int;
    typedef long s64Int;
#   define FSTR64 "%ld"
#   define FSTRU64 "%lu"
#   define ZERO64B 0L
#elif SIZEOF_LONG_LONG == 8
#   define POLY 0x0000000000000007ULL
#   define PERIOD 1317624576693539401LL
    typedef unsigned long long u64Int;
    typedef long long s64Int;
#   define FSTR64 "%lld"
#   define FSTRU64 "%llu"
#   define ZERO64B 0LL
#else
#   error long nor long long are 8 bytes in size
#endif

/* Macros for timing */
#define CPUSEC() (HPL_timer_cputime())
#define RTSEC() (MPI_Wtime())

#define MAX_TOTAL_PENDING_UPDATES 1024
#define LOCAL_BUFFER_SIZE MAX_TOTAL_PENDING_UPDATES
#define MAX_OUTSTANDING_HANDLES 64
extern u64Int **HPCC_Table;

#endif /* RA_COMMON_H_ */
