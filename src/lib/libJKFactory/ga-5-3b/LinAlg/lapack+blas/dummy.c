#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* dummy routine to make linalg.a legitimate when we do not compile lapack or blas */
void __dummy_linalg(int *a, int *b)
{
     if(a) *a=*b;
}
