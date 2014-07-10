#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#    include <stdio.h>
#endif
#if HAVE_SHMALLOC
#    define SHMALLOC shmalloc
#    define SHFREE   shfree
#endif
#if HAVE_SHARE_MALLOC
#    define SHMALLOC share_malloc
#    define SHFREE   share_free
#endif

extern char *SHMALLOC();
extern int SHFREE();

#define MAX_ADDR 20
static int next_id = 0;              /* Keep track of id */
static char *shaddr[MAX_ADDR];       /* Keep track of addresses */


char *CreateSharedRegion(long *id, long *size)
{
    char *temp;

    if (next_id >= MAX_ADDR) {
        Error("CreateSharedRegion: too many shared regions", (long) next_id);
    }

    if ( (temp = SHMALLOC((unsigned) *size)) == (char *) NULL) {
        Error("CreateSharedRegion: failed in SHMALLOC", (long) *size);
    }

    *id = next_id++;
    shaddr[*id] = temp;

    return temp;
}


long DetachSharedRegion(long id, long size, char *addr)
{
    /* This needs improving to make more robust */
    return SHFREE(addr);
}


long DeleteSharedRegion(long id)
{
    /* This needs improving to make more robust */
    return SHFREE(shaddr[id]);
}


char *AttachSharedRegion(long id, long size)
{
    Error("AttachSharedRegion: cannot do this on SEQUENT or BALANCE",
            (long) -1);
}
