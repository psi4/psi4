#ifndef __massert_h__
#define __massert_h__

#include <cassert>

/* #define massertl(cond,label) do { \ */
/* if(!(cond)) { goto label;}        \ */
/* } while(0) */

#define massertl(cond,label) do { \
    if(!(cond)) { assert(0);}	  \
} while(0)


#define massert(cond) massertl(cond,error)

enum TslError {
  TSL_OK = 0,
  TSL_ERR
};

#endif /*__massert_h__*/

