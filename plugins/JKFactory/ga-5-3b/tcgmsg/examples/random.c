#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*$Id: random.c,v 1.2 1995-02-02 23:24:25 d3g681 Exp $*/
double random_(seed)
     unsigned long *seed;
{
  *seed = *seed *1812433253 + 12345;
  return (double) ((*seed & 0x7fffffff) * 4.6566128752458e-10);
}
