#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/file.h>
#include <sys/mman.h>

extern char *strdup();
extern char *mktemp();

#define MAX_ID 20
static struct id_list_struct {
    char *addr;                      /* pointer to shmem region */
    unsigned size;                   /* size of region */
    char *filename;                  /* associated file name */
    int fd;                          /*            file descriptor */
    int status;                      /* = 1 if in use */
} id_list[MAX_ID];

static int next_id = 0;
static char template[] = "/tmp/SHMEM.XXXXXX";

char *CreateSharedRegion(long *id, long *size)
{
    char *temp;

    if (next_id == MAX_ID) {
        Error("CreateSharedRegion: MAX_ID exceeded ", MAX_ID);
    }
    *id = next_id;

#ifdef APOLLO
    id_list[*id].fd = -1;
#else
    if ( (temp = strdup(template)) == (char *) NULL) {
        Error("CreateSharedRegion: failed to get space for filename", 0);
    }

    /* Generate scratch file to identify region ... need to know this
       name to attach to the region so need to establish some policy
       before AttachtoSharedRegion can work */

    id_list[*id].filename = mktemp(temp);
    if ( (id_list[*id].fd = open(id_list[*id].filename, 
                    O_RDWR|O_CREAT, 0666)) < 0) {
        Error("CreateSharedRegion: failed to open temporary file",0);
    }
#endif

    id_list[*id].addr = mmap((caddr_t) 0, (unsigned *) size, 
            PROT_READ|PROT_WRITE, 
            MAP_ANON|MAP_SHARED, id_list[*id].fd, 0);
#ifdef APOLLO
    if (id_list[*id].addr == (char *) 0) {
        Error("CreateSharedRegion: mmap failed",-1);
    }
#else
    if (id_list[*id].addr == (char *) -1) {
        Error("CreateSharedRegion: mmap failed",-1);
    }
#endif

    id_list[*id].size = *size;
    id_list[*id].status = 1;

    next_id++;
    return id_list[*id].addr;
}

long DetachSharedRegion(long id, long size, char *addr)
{
    if ( (id < 0) || (id > next_id)) {
        return (long) -1;
    }

    if (id_list[id].status != 1) {
        return (long) -1;
    }

    id_list[id].status = 0;

    return (long) munmap(id_list[id].addr, 0);
}

long DeleteSharedRegion(long id)
{
    if ( (id < 0) || (id > next_id) ) {
        return (long) -1;
    }

    if (id_list[id].status != 1) {
        return (long) -1;
    }

    (void) DetachSharedRegion(id, 0, (char *) 0);

    if (id_list[id].fd >= 0) {
        (void) close(id_list[id].fd);
        (void) unlink(id_list[id].filename);
    }

    return (long) 0;
}

char *AttachSharedRegion(long id, long size)
{
    Error("AttachSharedRegion: need mods for this to work on CONVEX",
            (long) -1);
}

long DeleteSharedAll()
{
    long id;
    long status = 0;

    for (id=0; id<next_id; id++) {
        if (id_list[id].status == 1) {
            status += DeleteSharedRegion(id);
        }
    }

    if (status) {
        return (long) -1;
    } else {
        return (long) 0;
    }
}
