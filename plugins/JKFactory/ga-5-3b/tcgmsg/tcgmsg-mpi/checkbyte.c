#if HAVE_CONFIG_H
#   include "config.h"
#endif

unsigned char CheckByte(unsigned char *c, long n)
{
    /*
       unsigned char sum = (char) 0;
       while (n-- > 0) {
       sum  = sum ^ *c++;
       }

       return sum;
       */

    unsigned int sum = 0;
    unsigned int mask = 0xff;

    while (n-- > 0) {
        sum += (int) *c++;
    }

    sum = (sum + (sum>>8) + (sum>>16) + (sum>>24)) & mask;
    return (unsigned char) sum;
}
