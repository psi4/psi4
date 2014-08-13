#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: datatypes.c,v 1.9.10.1 2006-12-14 13:24:47 manoj Exp $
 * conversion of MA identifiers between C to Fortran data types 
 * Note that ga_type_c2f(MT_F_INT) == MT_F_INT 
 */

#include "gacommon.h"
#include "macommon.h"
#include "ga-papi.h"

Integer pnga_type_f2c(Integer type)
{
Integer ctype;
   switch(type){
   case MT_F_INT: 
#       if   SIZEOF_F77_INTEGER == SIZEOF_INT
                ctype = C_INT;
#       elif SIZEOF_F77_INTEGER == SIZEOF_LONG
                ctype = C_LONG;
#       elif SIZEOF_F77_INTEGER == SIZEOF_LONG_LONG
                ctype = C_LONGLONG;
#       else
#           error SIZEOF_F77_INTEGER == SIZEOF_???
#       endif
                break;
   case MT_F_REAL: 
#       if   SIZEOF_F77_REAL == SIZEOF_FLOAT
                ctype = C_FLOAT;
#       elif SIZEOF_F77_REAL == SIZEOF_DOUBLE
                ctype = C_DBL;
#       elif SIZEOF_F77_REAL == SIZEOF_LONG_DOUBLE
#           error SIZEOF_F77_REAL == SIZEOF_LONG_DOUBLE not supported?
#       else
#           error SIZEOF_F77_REAL == SIZEOF_???
#       endif
                break;
   case MT_F_DBL: 
#       if   SIZEOF_F77_DOUBLE_PRECISION == SIZEOF_DOUBLE
                ctype = C_DBL;
#       elif SIZEOF_F77_DOUBLE_PRECISION == SIZEOF_LONG_DOUBLE
#           error SIZEOF_F77_DOUBLE_PRECISION == SIZEOF_LONG_DOUBLE not supported?
#       else
#           error SIZEOF_F77_DOUBLE_PRECISION == SIZEOF_???
#       endif
                break;
   case MT_F_DCPL: 
		ctype = C_DCPL;
                break;
   case MT_F_SCPL: 
#       if   SIZEOF_F77_REAL == SIZEOF_FLOAT
                ctype = C_SCPL;
#       elif SIZEOF_F77_REAL == SIZEOF_DOUBLE
                ctype = C_DCPL;
#       elif SIZEOF_F77_REAL == SIZEOF_LONG_DOUBLE
#           error SIZEOF_F77_REAL == SIZEOF_LONG_DOUBLE not supported?
#       else
#           error SIZEOF_F77_REAL == SIZEOF_???
#       endif
                break;
   default:     ctype = type;
                break;
   }
   
   return(ctype);
}


Integer pnga_type_c2f(Integer type)
{
Integer ftype;
   switch(type){
   case C_INT: 
                ftype = (sizeof(int) != sizeof(Integer))? -1: MT_F_INT;
                break;
   case C_LONG: 
                ftype = (sizeof(long) != sizeof(Integer))? -1: MT_F_INT;
                break;
   case C_LONGLONG: 
                ftype = (sizeof(long long) != sizeof(Integer))? -1: MT_F_INT;
                break;
   case C_FLOAT:
#       if   SIZEOF_FLOAT == SIZEOF_F77_REAL
                ftype = MT_F_REAL; 
#       elif SIZEOF_FLOAT == SIZEOF_F77_DOUBLE_PRECISION
                ftype = MT_F_DBL; 
#       else
                ftype = -1;
#       endif
                break;
   case C_DBL: 
                ftype = MT_F_DBL;
                break;
   case C_DCPL:
                ftype = MT_F_DCPL;
                break;
   case C_SCPL:
#       if   SIZEOF_FLOAT == SIZEOF_F77_REAL
                ftype = MT_F_SCPL;
#       elif SIZEOF_FLOAT == SIZEOF_F77_DOUBLE_PRECISION
                ftype = MT_F_DCPL;
#       else
                ftype = -1;
#       endif
                break;
   default:     ftype = type;
                break;
   }
   
   return(ftype);
}
