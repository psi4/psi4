/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/shmem.h,v 1.4 1995-02-24 02:17:44 d3h325 Exp $ */

/*
  Header file which declares stubs for the shared memory interface.

  Note that the input arguments switch between integers and pointers
  to integers depending on if they are modified on return.
*/


/*
  Create a shared region of at least size bytes, returning the actual size,
  the id associated with the region. The return vaue is a pointer to the
  the region. Any error is a hard fail.

  (char *) CreateSharedRegion((long *) id, (long *) size)
*/
extern char *CreateSharedRegion();

/*
  Detach a process from a shared memory region. 0 is returned on success,
  -1 for failure. id, size, and addr much match exactly those items returned
  from CreateSharedRegion

  long DetachSharedRegion((long) id, (long) size, (char *) addr)
*/
extern long DetachSharedRegion();


/*
  Delete a shared region from the system. This has to be done on the SUN
  to remove it from the system. On the Alliant the shared region disappears
  when the last process dies or detaches. Returns 0 on success, -1 on error.

  long DeleteSharedRegion( (long) id)
*/
extern long DeleteSharedRegion();


/*
  Attach to a shared memory region of known id and size. Returns the
  address of the mapped memory. Size must exactly match the size returned
  from CreateSharedRegion (which in turn is the requested size rounded
  up to a multiple of 4096). Any error is a hard fail. 

  (char *) AttachSharedRegion( (long) id, (long) size))
*/
extern char *AttachSharedRegion();
