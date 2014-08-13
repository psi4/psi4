#ifndef __CONFIG_H__
#define __CONFIG_H__


#include <stdio.h>


#define _DEBUG_LEVEL_    0   //  0 to 10,
                             //  0 is no debug print info at all,
                             //  10 is full info


#if ( _DEBUG_LEVEL_ == -1 )
#define CINT_PRINTF( level, fmt, args... )        {}
#else
#define CINT_PRINTF( level, fmt, args... )                                          \
        do                                                                          \
        {                                                                           \
            if ( (unsigned)(level) <= _DEBUG_LEVEL_ )                               \
            {                                                                       \
                sprintf( basis->str_buf, "%s() line %d ", __FUNCTION__, __LINE__ ); \
                sprintf( basis->str_buf + strlen(basis->str_buf), fmt, ##args );      \
                fprintf( stdout, "%s", basis->str_buf );                            \
                fflush( stdout );                                                   \
            }                                                                       \
        } while ( 0 )
#endif


#if ( _DEBUG_LEVEL_ > 1 )
#define CINT_INFO( fmt, args... )                                            \
        do                                                                   \
        {                                                                    \
            sprintf( basis->str_buf, "**** CInt: ");                         \
            sprintf( basis->str_buf + strlen("**** CInt: "), fmt, ##args );  \
            fprintf( stdout, "%s", basis->str_buf );                         \
            fflush( stdout );                                                \
        } while ( 0 )
#else
#define CINT_INFO( fmt, args... )        {}
#endif


#endif /* __CONFIG_H__ */
