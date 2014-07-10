#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* The program was written to address missing "ipcrm" command on Mac X
 * It probes a range of System V shared memory id from 0 to MAXID
 * and if exist, it attempts to delete them.
 */
#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif
#if HAVE_SYS_IPC_H
#   include <sys/ipc.h>
#endif
#if HAVE_SYS_SHM_H
#   include <sys/shm.h>
#endif
#if HAVE_SYS_PARAM_H
#   include <sys/param.h>
#endif
#if HAVE_ERRNO_H
#   include <errno.h>
#endif
#if HAVE_STDIO_H
#   include <stdio.h>
#endif

#ifdef SUN
char *shmat();
#endif

#define MAXID 1000000
int main(int argc, char **argv)
{
  int i;
  for (i = 0; i < MAXID; i++) {
#if SIZEOF_VOIDP == SIZEOF_INT
    int rc = (int) shmat(i, (char *) NULL, 0);
#elif SIZEOF_VOIDP == SIZEOF_LONG
    long rc = (long) shmat(i, (char *) NULL, 0);
#elif SIZEOF_VOIDP == SIZEOF_LONGLONG
    long long rc = (long long) shmat(i, (char *) NULL, 0);
#endif
    if (rc < 0) {
      continue;
    }
    printf("found %d\n", i);
    shmdt((void *)rc);
    /* delete segment id */
    if (shmctl(i, IPC_RMID, (struct shmid_ds *)NULL)) {
      printf("failed to remove shm id=%d\n", i);
    }
  }

  return 0;
}
