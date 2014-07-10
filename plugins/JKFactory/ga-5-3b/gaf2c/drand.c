#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#include "typesf2c.h"

#ifdef CRAY_YMP
/** on YMP/J90 need to use thread safe version of rand */
static double ran(unsigned int flag)
{
    static unsigned long seed = 76521;

    if(flag != 0) {
        seed = flag;
    }

    seed = seed *1812433253 + 12345;

    return ((double) (seed & 0x7fffffff)) * 4.6566128752458e-10;
}
#endif

static DoublePrecision gai_drand_(Integer *flag)
{
#ifdef CRAY_YMP
    return ran((unsigned int)*flag);
#else
    if (*flag)
        srandom((unsigned) *flag);

    return ((DoublePrecision) random()) * 4.6566128752458e-10;
#endif
}

#define drand_ F77_FUNC(drand,DRAND)
DoublePrecision drand_(Integer *flag)
{
    return (gai_drand_(flag));
}
