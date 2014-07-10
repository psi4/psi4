#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file
*/

/**
 * Cut from SUN man-page.
 *
 * memalign() allocates size bytes  on  a  specified  alignment
 * boundary, and returns a pointer to the allocated block.  The
 * value of the returned address is guaranteed to  be  an  even
 * multiple of alignment.  Note: the value of alignment must be
 * a power of two, and must be greater than  or  equal  to  the
 * size of a word.
 *
 * No checking is done on the value of alignment ... should really.
 */
char *memalign(unsigned alignment, unsigned size)
{
    union screwup {
        unsigned long integer;
        char *address;
    } fiddle;

    unsigned long offset;

    alignment += alignment;   /* Actually align on twice requested boundary */

    fiddle.address = malloc((unsigned) (alignment+size));

    if (fiddle.address != (char *) 0) {
        offset = fiddle.integer & (alignment-1);
        if (offset != 0) {
            fiddle.address += alignment - offset;
        }
    }

    return fiddle.address;
}
