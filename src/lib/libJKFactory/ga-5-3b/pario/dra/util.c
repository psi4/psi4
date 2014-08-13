#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "ga.h"
#include "typesf2c.h"

/* Return the no. of bytes that n doubles occupy */
#define util_mdtob_ F77_FUNC_(util_mdtob,UTIL_MDTOB)
Integer FATR util_mdtob_(Integer *n)
{
    if (*n < 0)
        GA_Error("util_MDTOB_: negative argument",*n);

    return (Integer) (*n * sizeof(double));
}


/* Return the no. of bytes that n ints=Integers occupy */
#define util_mitob_ F77_FUNC_(util_mitob,UTIL_MITOB)
Integer FATR util_mitob_(Integer *n)
{
    if (*n < 0)
        GA_Error("util_MITOB_: negative argument",*n);

    return (Integer) (*n * sizeof(Integer));
}
