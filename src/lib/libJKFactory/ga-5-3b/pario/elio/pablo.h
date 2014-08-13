/* Header file for PABLO instrumentation */

#if defined(PABLO)
#  define IOTRACE
#  include "IOTrace.h"
#endif


/*  Pablo profiler definitions */

#  define PABLO_elio_write    710000
#  define PABLO_elio_awrite    710001
#  define PABLO_elio_read    710002
#  define PABLO_elio_aread    710003
#  define PABLO_elio_wait    710004
#  define PABLO_elio_probe    710005
#  define PABLO_elio_stat    710006
#  define PABLO_elio_open    710007
#  define PABLO_elio_close    710009
#  define PABLO_elio_set_cb    710010
#  define PABLO_elio_delete    710011
#  define PABLO_elio_truncate    710012
#  define PABLO_elio_length    710014


#if defined(PABLO)
#  define PABLO_init        initIOTrace
#  define PABLO_start(_id)    startTimeEvent( _id )
#  define PABLO_end(_id)    endTimeEvent( _id )
#  define PABLO_terminate    {endIOTrace(); endTracing();}
#else
#  define PABLO_init 
#  define PABLO_start(_id) 
#  define PABLO_end( _id ) 
#  define PABLO_terminate  
#endif
