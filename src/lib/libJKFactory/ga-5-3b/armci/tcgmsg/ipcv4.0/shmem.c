#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/shmem.c,v 1.13 2000-10-13 20:55:40 d3h325 Exp $ */

/*
  This stuff attempts to provide a simple interface to temporary shared
  memory regions, loosely modelled after that of Alliant Concentrix 5.0 


  Note that the input arguments switch between integers and pointers
  to integers depending on if they are modified on return.


  Create a shared region of at least size bytes, returning the actual size,
  the id associated with the region. The return value is a pointer to the
  the region. Any error is a hard fail.

  (char *) CreateSharedRegion((long *) id, (long *) size)


  Detach a process from a shared memory region. 0 is returned on success,
  -1 for failure. id, size, and addr must match exactly those items returned
  from CreateSharedRegion

  long DetachSharedRegion((long) id, (long) size, (char *) addr)


  Delete a shared region from the system. This has to be done on the SUN
  to remove it from the system. On the Alliant the shared region disappears
  when the last process dies or detaches. Returns 0 on success, -1 on error.

  long DeleteSharedRegion((long) id)


  Delete all the shared regions associated with this process.

  long DeleteSharedAll()


  Attach to a shared memory region of known id and size. Returns the
  address of the mapped memory. Size must exactly match the size returned
  from CreateSharedRegion (which in turn is the requested size rounded
  up to a multiple of 4096). Any error is a hard fail. 

  (char *) AttachSharedRegion((long) id, (long) size))

*/

extern void Error();

#ifdef ALLIANT

#include <stdio.h>
#include <sys/time.h>

extern char *valloc();

char *CreateSharedRegion(id, size)
     long *size, *id;
{
  struct timeval tp;
  struct timezone tzp;
  char *temp;
  int status;

  /* Have to round up to a multiple of page size before allocating
     on a page boundary */

  *size = ( (*size + 4095) / 4096 ) * 4096;

  if ( (temp = valloc((unsigned) *size)) == (char *) NULL)
    Error("CreateSharedRegion: failed in valloc", (long) 0);

  /* Now have to get a unique id ... try using time of day in centi-sec */

  if ( (status = gettimeofday(&tp, &tzp)) != 0)
    Error("CreateSharedRegion: error from gettimeofday", (long) status);

  *id = (tp.tv_sec + 10000*tp.tv_usec) & 0xffffff;

  /* Now make the region */

  if ( (status = create_shared_region(*id, temp, *size, 0)) != 0)
    Error("CreateSharedRegion: error from create_shared_region", (long) status);

  return temp;
}

long DetachSharedRegion( id, size, addr)
     long id, size;
     char *addr;
{
  return detach_shared_region( id, addr, size);
}

long DeleteSharedRegion(id)
     long id;
{
  return delete_shared_region(id);
}

char *AttachSharedRegion(id, size)
     long id, size;
{
  char *temp;
  int status;

  if (size !=  (((size + 4095) / 4096) * 4096))
    Error("AttachSharedRegion: input size is not multiple of 4096",
          (long) size);

  if ( (temp = valloc((unsigned) size)) == (char *) NULL)
    Error("AttachSharedRegion: failed in valloc", (long) 0);

  /* Now try to attach */

  if ( (status = attach_shared_region(id, temp, size)) != 0)
    Error("AttachSharedRegion: error from attach_shared_region",
          (long) status);

  return temp;
}

#endif
#if defined(SEQUENT) || defined(ENCORE) /* @**!ing SEQUENT and CRAY no elif */

#include <stdio.h>

#ifdef SEQUENT
#define SHMALLOC shmalloc
#define SHFREE   shfree
#endif
#ifdef ENCORE
#define SHMALLOC share_malloc
#define SHFREE   share_free
#endif

extern char *SHMALLOC();
extern int SHFREE();

#define MAX_ADDR 20
static int next_id = 0;              /* Keep track of id */
static char *shaddr[MAX_ADDR];       /* Keep track of addresses */

char *CreateSharedRegion(id, size)
     long *size, *id;
{
  char *temp;

  if (next_id >= MAX_ADDR)
	Error("CreateSharedRegion: too many shared regions", (long) next_id);

  if ( (temp = SHMALLOC((unsigned) *size)) == (char *) NULL)
    Error("CreateSharedRegion: failed in SHMALLOC", (long) *size);

  *id = next_id++;
  shaddr[*id] = temp;

  return temp;
}

/*ARGSUSED*/
long DetachSharedRegion( id, size, addr)
     long id, size;
     char *addr;
{
  /* This needs improving to make more robust */
  return SHFREE(addr);
}

long DeleteSharedRegion(id)
     long id;
{
  /* This needs improving to make more robust */
  return SHFREE(shaddr[id]);
}

/*ARGSUSED*/
char *AttachSharedRegion(id, size)
     long id, size;
{
  Error("AttachSharedRegion: cannot do this on SEQUENT or BALANCE", (long) -1);
}


#endif
   /* Bizarre sequent has sysv semaphores but proprietary shmem */
   /* Encore has sysv shmem but is limited to total of 16384bytes! */
#if defined(SYSV) && !defined(SEQUENT) && !defined(ENCORE)

#include <stdio.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

char *CreateSharedRegion(id, size)
     long *size, *id;
{
  char *temp;

  /* Create the region */

  if ( (*id = shmget(IPC_PRIVATE, (int) *size, 
                     (int) (IPC_CREAT | 00600))) < 0 )
    Error("CreateSharedRegion: failed to create shared region", (long) *id);

  /* Attach to the region */

  if ( (long) (temp = shmat((int) *id, (char *) NULL, 0)) == -1L)
    Error("CreateSharedRegion: failed to attach to shared region", (long) 0);

  return temp;
}

/*ARGSUSED*/
long DetachSharedRegion( id, size, addr)
     long id, size;
     char *addr;
{
  return shmdt(addr);
}

long DeleteSharedRegion(id)
     long id;
{
  return shmctl((int) id, IPC_RMID, (struct shmid_ds *) NULL);
}

/*ARGSUSED*/
char *AttachSharedRegion(id, size)
     long id, size;
{
  char *temp;

  if ( (long) (temp = shmat((int) id, (char *) NULL, 0)) == -1L)
    Error("AttachSharedRegion: failed to attach to shared region", (long) 0);

  return temp;
}

#endif
#if (defined(CONVEX) || defined(APOLLO)) && !defined(HPUX) 

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

char *CreateSharedRegion(id, size)
     long *size, *id;
{
  char *temp;

  if (next_id == MAX_ID)
    Error("CreateSharedRegion: MAX_ID exceeded ", MAX_ID);
  *id = next_id;

#ifdef APOLLO
  id_list[*id].fd = -1;
#else
  if ( (temp = strdup(template)) == (char *) NULL)
    Error("CreateSharedRegion: failed to get space for filename", 0);

/* Generate scratch file to identify region ... need to know this
   name to attach to the region so need to establish some policy
   before AttachtoSharedRegion can work */

  id_list[*id].filename = mktemp(temp);
  if ( (id_list[*id].fd = open(id_list[*id].filename, 
                                   O_RDWR|O_CREAT, 0666)) < 0)
    Error("CreateSharedRegion: failed to open temporary file",0);
#endif

  id_list[*id].addr = mmap((caddr_t) 0, (unsigned *) size, 
                           PROT_READ|PROT_WRITE, 
                           MAP_ANON|MAP_SHARED, id_list[*id].fd, 0);
#ifdef APOLLO
  if (id_list[*id].addr == (char *) 0)
    Error("CreateSharedRegion: mmap failed",-1);
#else
  if (id_list[*id].addr == (char *) -1)
    Error("CreateSharedRegion: mmap failed",-1);
#endif

  id_list[*id].size = *size;
  id_list[*id].status = 1;

  next_id++;
  return id_list[*id].addr;
}

/*ARGSUSED*/
long DetachSharedRegion( id, size, addr)
     long id, size;
     char *addr;
{
 if ( (id < 0) || (id > next_id))
   return (long) -1;

 if (id_list[id].status != 1)
   return (long) -1;

 id_list[id].status = 0;

 return (long) munmap(id_list[id].addr, 0);
}

long DeleteSharedRegion(id)
     long id;
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

/*ARGSUSED*/
char *AttachSharedRegion(id, size)
     long id, size;
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

#endif
