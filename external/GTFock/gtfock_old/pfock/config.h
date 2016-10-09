#ifndef __CONFIG_H__
#define __CONFIG_H__


#define GA_NB
#define MIN(a, b)    ((a) < (b) ? (a) : (b))
#define MAX(a, b)    ((a) > (b) ? (a) : (b))

#define _DEBUG_LEVEL_    2   //  0 to 10,
                             //  0 is no debug print info at all,
                             //  10 is full info

#define alignsize  64
#ifdef __INTEL_COMPILER
   #define PFOCK_MALLOC(size)    _mm_malloc(size, alignsize)
   #define PFOCK_FREE(addr)      _mm_free(addr)
#else
   #define PFOCK_MALLOC(size)    malloc(size)
   #define PFOCK_FREE(addr)      free(addr)
#endif

#if ( _DEBUG_LEVEL_ == -1 )
#define PFOCK_PRINTF( level, fmt, args... )        {}
#else
#define PFOCK_PRINTF( level, fmt, args... )                                \
        do                                                                 \
        {                                                                  \
            if ( (unsigned)(level) <= _DEBUG_LEVEL_ )                      \
            {                                                              \
                fprintf( stdout, "%s() line %d: "fmt,                      \
                    __FUNCTION__, __LINE__, ##args);                       \
                fflush( stdout );                                          \
            }                                                              \
        } while ( 0 )
#endif


#if ( _DEBUG_LEVEL_ > 1 )
#define PFOCK_INFO( fmt, args... )                                         \
        do                                                                 \
        {                                                                  \
            fprintf( stdout, "**** PFock: "fmt, ##args );                  \
            fflush( stdout );                                              \
        } while ( 0 )
#else
#define PFOCK_INFO( fmt, args... )        {}
#endif


#endif /* __CONFIG_H__ */
