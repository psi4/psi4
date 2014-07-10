#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* OS specific server process/thread creation and destruction  
 * JN/03.25.2000
 */
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_ERRNO_H
#   include <errno.h>
#endif
#include "armcip.h"

#ifdef WIN32
/************************** Windows threads **************************/
#if HAVE_WINDOWS_H
#   include <windows.h>
#endif
#if HAVE_PROCESS_H
#   include <process.h>
#endif

thread_id_t armci_serv_tid;
unsigned long armci_serv_handle;
#ifndef NO_STDLIBC
#define NEWTHREAD CreateThread
#else
#define NEWTHREAD _beginthreadex 
#endif

unsigned __stdcall armci_wrap_func(void *arg)
{
void (*func)(void*);
     func = arg;

     /* boost the server thread priority be better responsiveness */
     (void)SetThreadPriority(GetCurrentThread(),THREAD_PRIORITY_HIGHEST);

     func(NULL);

     return (unsigned)0;
}


void armci_create_server_thread ( void* (* func)(void*) )
{
     /* as we need to use std C rt library we cannot use CreateThread */
     armci_serv_handle = NEWTHREAD(NULL, 0, armci_wrap_func, (void*)func, 
                                         0, &armci_serv_tid);
     if(!armci_serv_handle) 
         armci_die("armci_create_server_thread: create failed",0);



}

void armci_terminate_server_thread()
{
/*int rc;*/
/*    TerminateThread(armci_serv_handle,&rc);*/
}

/****************************** PTHREADS *****************************/
#elif defined(PTHREADS)
#include <pthread.h>

thread_id_t armci_serv_tid;

void armci_create_server_thread ( void* (* func)(void*) )
{
pthread_attr_t attr;
int rc;

    if(pthread_attr_init(&attr))
       armci_die("armci_create_server_thread: attr init failed",0);

#if defined(AIX) || defined(SOLARIS)
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
#endif

    rc = pthread_create(&armci_serv_tid, &attr, func, NULL);
    if(rc) armci_die("armci_create_server_thread: create failed",errno);

    pthread_attr_destroy(&attr);
}


void armci_terminate_server_thread()
{
  if(pthread_join(armci_serv_tid,NULL))
     armci_die("armci_terminate_server_thread: failed",0);
}


#else
/**************************** Unix processes ******************************/

#if HAVE_UNISTD_H
#   include <unistd.h>
#endif
#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif
#if HAVE_SYS_WAIT_H
#   include <sys/wait.h>
#endif

pid_t server_pid= (pid_t)0;


void armci_create_server_process ( void* (* func)(void*) )
{
pid_t pid;     

     if ( (pid = fork() ) < 0)
        armci_die("fork failed", (int)pid);

     else if(pid == 0){

        armci_me = SOFFSET - armci_me; /* server id derived from parent id */

        func(NULL);

     } else server_pid = pid;
}


void armci_wait_server_process()
{
     int stat;
     pid_t rc;

     if(!server_pid) return;
     rc = wait (&stat);
     if (rc != server_pid){
         perror("ARMCI master: wait for child process (server) failed:");
     }
     server_pid = (pid_t)0;
}

#endif
