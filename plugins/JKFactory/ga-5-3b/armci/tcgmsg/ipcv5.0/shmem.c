#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file
 *
 * $Header: /tmp/hpctools/ga/tcgmsg/ipcv5.0/shmem.c,v 1.4 2002-01-24 22:07:27 d3h325 Exp $
 *
 * This stuff attempts to provide a simple interface to temporary shared
 * memory regions, loosely modelled after that of Alliant Concentrix 5.0 
 *
 *
 * Note that the input arguments switch between integers and pointers
 * to integers depending on if they are modified on return.
 *
 *
 * Create a shared region of at least size bytes, returning the actual size,
 * the id associated with the region. The return value is a pointer to the
 * the region. Any error is a hard fail.
 *
 * (char *) CreateSharedRegion((long *) id, (long *) size)
 *
 *
 * Detach a process from a shared memory region. 0 is returned on success,
 * -1 for failure. id, size, and addr must match exactly those items returned
 * from CreateSharedRegion
 *
 * long DetachSharedRegion((long) id, (long) size, (char *) addr)
 *
 *
 * Delete a shared region from the system. This has to be done on the SUN
 * to remove it from the system. On the Alliant the shared region disappears
 * when the last process dies or detaches. Returns 0 on success, -1 on error.
 *
 * long DeleteSharedRegion((long) id)
 *
 *
 * Delete all the shared regions associated with this process.
 *
 * long DeleteSharedAll()
 *
 *
 * Attach to a shared memory region of known id and size. Returns the
 * address of the mapped memory. Size must exactly match the size returned
 * from CreateSharedRegion (which in turn is the requested size rounded
 * up to a multiple of 4096). Any error is a hard fail. 
 *
 * (char *) AttachSharedRegion((long) id, (long) size))
 */

extern void Error(const char *, long);

#if !defined(MMAP) || defined(MACX)

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif
#if HAVE_SYS_IPC_H
#   include <sys/ipc.h>
#endif
#if HAVE_SYS_SHM_H
#   include <sys/shm.h>
#endif

#ifdef SUN
extern int shmget(key_t, int, int);
extern int shmdt(void *);
extern int shmctl(int, int, struct shmid_ds *);
extern void *shmat(int, const void *, int);
#endif


char *CreateSharedRegion(long *id, long *size)
{
    char *temp;

    /* Create the region */

    if ( (*id = shmget(IPC_PRIVATE, (int) *size, 
                    (int) (IPC_CREAT | 00600))) < 0 )
        Error("CreateSharedRegion: failed to create shared region", (long) *id);

    /* Attach to the region */

    if ( (long) (temp = shmat((int) *id, (char *) NULL, 0)) == -1L)
        Error("CreateSharedRegion: failed to attach to shared region", 0L);

    return temp;
}


long DetachSharedRegion(long id, long size, char *addr)
{
    return shmdt(addr);
}


long DeleteSharedRegion(long id)
{
    return shmctl((int) id, IPC_RMID, (struct shmid_ds *) NULL);
}


char *AttachSharedRegion(long id, long size)
{
    char *temp;

    if ( (long) (temp = shmat((int) id, (char *) NULL, 0)) == -1L)
        Error("AttachSharedRegion: failed to attach to shared region", 0L);

    return temp;
}

#else /* MMAP */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_SYS_TIME_H
#   include <sys/time.h>
#endif
#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif
#if HAVE_SYS_FILE_H
#   include <sys/file.h>
#endif
#if HAVE_SYS_MMAN_H
#   include <sys/mman.h>
#endif

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

char *CreateSharedRegion(id, size)
    long *size, *id;
{
    char *temp;

    if (next_id == MAX_ID)
        Error("CreateSharedRegion: MAX_ID exceeded ", MAX_ID);
    *id = next_id;

    if ( (temp = strdup(template)) == (char *) NULL)
        Error("CreateSharedRegion: failed to get space for filename", 0);

    /* Generate scratch file to identify region ... need to know this
       name to attach to the region so need to establish some policy
       before AttachtoSharedRegion can work */

    id_list[*id].filename = mktemp(temp);
    if ( (id_list[*id].fd = open(id_list[*id].filename, 
                    O_RDWR|O_CREAT, 0666)) < 0)
        Error("CreateSharedRegion: failed to open temporary file",0);

    id_list[*id].addr = mmap((caddr_t) 0, (size_t)*size, 
            PROT_READ|PROT_WRITE, 
            MAP_ANON|MAP_SHARED, id_list[*id].fd, 0);
    if (id_list[*id].addr == (char *) -1)
        Error("CreateSharedRegion: mmap failed",-1);

    id_list[*id].size = *size;
    id_list[*id].status = 1;

    next_id++;
    return id_list[*id].addr;
}


long DetachSharedRegion(long id, long size, char *addr)
{
    if ( (id < 0) || (id > next_id))
        return (long) -1;

    if (id_list[id].status != 1)
        return (long) -1;

    id_list[id].status = 0;

    return (long) munmap(id_list[id].addr, 0);
}


long DeleteSharedRegion(long id)
{
    if ( (id < 0) || (id > next_id) )
        return (long) -1;

    if (id_list[id].status != 1)
        return (long) -1;

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

    for (id=0; id<next_id; id++)
        if (id_list[id].status == 1)
            status += DeleteSharedRegion(id);

    if (status)
        return (long) -1;
    else
        return (long) 0;
}

#endif /* MMAP */
