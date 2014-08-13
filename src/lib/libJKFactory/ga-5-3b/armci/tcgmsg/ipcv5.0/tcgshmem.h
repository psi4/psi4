/** @file
 * Header file which declares stubs for the shared memory interface.
 * Note that the input arguments switch between integers and pointers
 * to integers depending on if they are modified on return.
 */
#ifndef SHMEM_H_
#define SHMEM_H_

/**
 * Create a shared region of at least size bytes, returning the actual size,
 * the id associated with the region. The return vaue is a pointer to the
 * the region. Any error is a hard fail.
 */
extern char *CreateSharedRegion(long *id, long *size);

/**
 * Detach a process from a shared memory region. 0 is returned on success,
 * -1 for failure. id, size, and addr much match exactly those items returned
 * from CreateSharedRegion
 *
 */
extern long DetachSharedRegion(long id, long size, char *addr);

/**
 * Delete a shared region from the system. This has to be done on the SUN
 * to remove it from the system. On the Alliant the shared region disappears
 * when the last process dies or detaches. Returns 0 on success, -1 on error.
 */
extern long DeleteSharedRegion(long id);

/**
 * Attach to a shared memory region of known id and size. Returns the
 * address of the mapped memory. Size must exactly match the size returned
 * from CreateSharedRegion (which in turn is the requested size rounded
 * up to a multiple of 4096). Any error is a hard fail. 
 */
extern char *AttachSharedRegion(long id, long size);

#endif /* SHMEM_H_ */
