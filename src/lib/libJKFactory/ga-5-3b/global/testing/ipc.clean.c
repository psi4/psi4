#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif
#if HAVE_SYS_IPC_H
#   include <sys/ipc.h>
#endif
#if HAVE_SYS_SEM_H
#   include <sys/sem.h>
#endif
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_ERRNO_H
#   include <errno.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#if !HAVE_UNION_SEMUN
union semun {
        int val;                    /* value for SETVAL */
        struct semid_ds *buf;       /* buffer for IPC_STAT, IPC_SET */
        unsigned short int *array;  /* array for GETALL, SETALL */
        struct seminfo *__buf;      /* buffer for IPC_INFO */
};
#endif

#define MAX_SEM  10 
 
struct sembuf sops;
int semaphoreID;
int sem_init=0;

#define P_      -1
#define V_       1
#define P(s)  \
{\
  sops.sem_num = (s);\
  sops.sem_op  =  P_;\
  sops.sem_flg =  0; \
  semop(semaphoreID,&sops,1);\
}
#define V(s) \
{\
  sops.sem_num = (s);\
  sops.sem_op  =  V_;\
  sops.sem_flg =  0; \
  semop(semaphoreID,&sops,1);\
}


int SemGet(num_sem)
    int num_sem;
{
    if(num_sem<1)return(0);
    if(num_sem>MAX_SEM)return(0);
 
    semaphoreID = semget(IPC_PRIVATE,num_sem,0600);
    if(semaphoreID<0){
       fprintf(stderr,"SemGet failed \n");
       perror((char*)0);
    }
       
    sem_init = num_sem;
    return(semaphoreID);
}

void SemInit(id,value)
    int id,value;
{
    union semun semctl_arg;
    fprintf(stderr,"SemInit %d %d\n",id,value);
   
    semctl_arg.val = value;
    if(id >= sem_init || id<0 ) 
      fprintf(stderr,"attempt to intialize invalid semaphore %d %d\n",
                                                         id,sem_init);
    else if( semctl(semaphoreID, id,SETVAL,semctl_arg )<0){ 
         fprintf(stderr,"SemInit error\n");
         perror((char*)0);
    }
    fprintf(stderr,"exiting SemInit \n");
}


/*  release semaphore(s) */
void SemDel()
{
     semctl(semaphoreID,NULL,IPC_RMID,NULL);
}



/*\
 * (char *) CreateSharedRegion((long *) id, (long *) size)
 * long DetachSharedRegion((long) id, (long) size, (char *) addr)
 * long DeleteSharedRegion((long) id)
 * long DeleteSharedAll()
 * (char *) AttachSharedRegion((long) id, (long) size))
\*/

void Error( str, code)
     char *str;
     int code;
{
fprintf(stderr,"%s %d\n",str, code);
exit(0);
}

#ifdef ALLIANT

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
  *size = ( (*size + (PAGE_SIZE -1)) / PAGE_SIZE ) * PAGE_SIZE;

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

  if (size !=  (((size + (PAGE_SIZE -1)) / PAGE_SIZE) * PAGE_SIZE))
    Error("AttachSharedRegion: input size is not multiple of PAGE_SIZE",
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
#if defined(SEQUENT) || defined(ENCORE)

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

char *AttachSharedRegion(id, size)
     long id, size;
{
  Error("AttachSharedRegion: cannot do this on SEQUENT or BALANCE", (long) -1);
}


#endif
   /* Bizarre sequent has sysv semaphores but proprietary shmem */
   /* Encore has sysv shmem but is limited to total of 16384bytes! */
#if defined(SYSV) && !defined(SEQUENT) && !defined(ENCORE)

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <errno.h>
#include <stdio.h>


#ifdef SUN
extern char *shmat();
#endif

char *CreateSharedRegion(id, size)
     long *size, *id;
{
  char *temp;

  /* Create the region */
  if ( (*id = shmget(IPC_PRIVATE, (int) *size, 
                     (int) (IPC_CREAT | 00600))) < 0 ){
    fprintf(stderr,"id=%d size=%d\n",*id, (int) *size);
    perror((char*)0);
    Error("CreateSharedRegion: failed to create shared region", (long) *id);
  }

  /* Attach to the region */
  if ( (temp = shmat((int) *id, (char *) NULL, 0)) == (char *) NULL){
    perror((char*)0);
    Error("CreateSharedRegion: failed to attach to shared region", (long) 0);
  }

  return temp;
}

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

char *AttachSharedRegion(id, size)
     long id, size;
{
  char *temp;

  if ( (temp = shmat((int) id, (char *) NULL, 0)) == (char *) NULL)
    Error("AttachSharedRegion: failed to attach to shared region", (long) 0);

  return temp;
}

#endif
#if defined(CONVEX) || defined(APOLLO)

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

char *AttachSharedRegion(id, size)
     long id, size;
{
  Error("AttachSharedRegion: need mods for this to work on CONVEX", (long) -1);
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
int main(int argc, char **argv)
{
int from=0, to, i;
    if(argc<2){
      printf("Usage:\n ipc.clean [<from>] <to> \n single argument is interpreted as <to> with <from> = 0 assumed\n");
      return 1;
    }
    if(argc=2) sscanf(argv[1],"%d",&to);
    else {
         sscanf(argv[1],"%d",&from);
         sscanf(argv[2],"%d",&to);
    }
    if(from>to && to <0){
       printf("wrong arguments\n");
       return 1;
    }
    for(i=from;i<=to;i++){ 
      semaphoreID =i;
      SemDel();
      DeleteSharedRegion((long)i);
    }

    return 0;
}
