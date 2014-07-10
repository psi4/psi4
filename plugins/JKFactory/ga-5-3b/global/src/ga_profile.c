#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: ga_profile.c,v 1.5 2005-07-21 08:13:26 manoj Exp $ */
/**
 * Note #1: Right now, only process 0's profile is printed.
 * Each and every process saves its profile in the correspoding data struture.
 * However profiler prints process 0's profile when ga_profile_terminate()
 * is called. Do the corresponding changes in ga_profile_terminate() to
 * print the profile of other processes.
 *
 * Note #2: By default profiles prints message ranges #s 21. Example: range 10
 * corresponds to message ranges from 1024 bytes to 2047 bytes.
 * Message ranges are in the power of 2. for ex:
 * ------------------------------------
 *  MSG_RANGE (r)        BYTES (2^r to 2^(r+1)-1)
 * ------------------------------------
 *      0                    0-1
 *      1                    2-3
 *      2                    4-7
 *     ...                   ...
 *      10                1024-2047 bytes
 *     ...                   ...
 *      20                1MB - (2MB-1)
 *      21                  >= 2MB
 * -------------------------------------
 *
 * Note#3: If Stride information needs to be printed, set GA_PRINT_STRIDE.
 * Stride information is printed in ga_profile_terminate() for a various
 * selective message ranges and event types. Modify according to your needs.
 */


#ifdef ENABLE_PROFILE

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif
#include "globalp.h"
#include "base.h" 
#include "ga_profile.h" 
#include "ga-papi.h"
#include "ga-wapi.h"

#ifndef MPI
#  include "tcgmsg.h"
#   define MP_TIMER tcg_time
#else
#  include "mpi.h"
#   define MP_TIMER MPI_Wtime
#endif


#define GA_PRINT_STRIDE 1
#define GA_MAX_MSG_RANGE 21

#if GA_PRINT_STRIDE
#define STRIDE_COUNT 1000
  typedef struct ga_stride {
    int ndim;
    int lo[GA_MAX_DIM];
    int hi[GA_MAX_DIM];
    char name[FNAM+1];
    double time;
  }ga_stride_t;
#endif

#define GA_EVENTS 6 /*  get, put, acc, Non-Contiguous get, put, acc*/
enum events {GET,    /* Contiguous get */
	     PUT, 
	     ACC, 
	     NC_GET, /* Non contiguous Get */
	     NC_PUT,
	     NC_ACC
};

char *event_name[GA_EVENTS] = {"GET", "PUT", "ACC", "NON CONTIGUOUS GET",
                               "NON CONTIGUOUS PUT", "NON CONTIGUOUS ACC"};

typedef struct ga_profile {
  int count;          /* number of times called */
  double exectime;  /* total execution time for "count" calls */
#if GA_PRINT_STRIDE
  ga_stride_t stride[STRIDE_COUNT];
#endif
}ga_profile_t;

/* profile get/put/acc for various message ranges (i.e GA_MAX_MSG_RANGE) */
static ga_profile_t GA_PROF[GA_EVENTS][GA_MAX_MSG_RANGE]; 

/* Current event */
struct event_info {
  int event_type;
  int range;
  int is_set;
  double start_time;
} gCURRENT_EVNT; 

void ga_profile_init() {
    int i,j;
    if(pnga_nodeid()==0) {printf("\nProfiling GA - ON\n");fflush(stdout);}
    for(i=0; i<GA_EVENTS; i++)
       for(j=0; j<GA_MAX_MSG_RANGE; j++) {
	  GA_PROF[i][j].count = 0;  GA_PROF[i][j].exectime = 0.0; 
       }
}

static void ga_profile_set_event(int event_type, int range) {
    gCURRENT_EVNT.event_type = event_type;
    gCURRENT_EVNT.range      = range;
    gCURRENT_EVNT.is_set     = 1;
    gCURRENT_EVNT.start_time = MP_TIMER();
}

void ga_profile_start(int g_a, long bytes, int ndim, Integer *lo, Integer *hi,
                      int comm_type) {
    int i, count=0, non_contig=0, event_type, range;

    /* find the message range */
    if(bytes<=0) range=0;
    else range = (int) (log((double)bytes)/log(2.0));
    if(range>=GA_MAX_MSG_RANGE) range = GA_MAX_MSG_RANGE;
 
    /* check contiguous or non-contiguous */
    for(i=0; i<ndim; i++) if(hi[0]-lo[0]) count++;
    if(count>1) non_contig=1; /* i.e. non-contiguous */
 
    switch(comm_type) {
       case ENABLE_PROFILE_PUT:
	  if(non_contig) event_type = NC_PUT;
	  else event_type = PUT;
	  break;
       case ENABLE_PROFILE_GET: 
	  if(non_contig) event_type = NC_GET;
	  else event_type = GET;
	  break;
       case ENABLE_PROFILE_ACC: 
	  if(non_contig) event_type = NC_ACC;
	  else event_type = ACC;
	  break;
       default: pnga_error("ENABLE_PROFILE: Invalid communication type", 0L);
    }

    /* set the curent event for timer */
    ga_profile_set_event(event_type, range);
    
    /* profile update: i.e. update event count */
    GA_PROF[event_type][range].count++;

#if GA_PRINT_STRIDE 
    {
       int idx = GA_PROF[event_type][range].count-1;
       if(idx<STRIDE_COUNT) {
	  GA_PROF[event_type][range].stride[idx].ndim = ndim;
	  strcpy(GA_PROF[event_type][range].stride[idx].name, GA[g_a].name);
	  for(i=0;i<ndim;i++) {
	     GA_PROF[event_type][range].stride[idx].lo[i] = (int)lo[i];
	     GA_PROF[event_type][range].stride[idx].hi[i] = (int)hi[i];
	  }
       }
    }
#endif
}

void ga_profile_stop() {
    int event_type = gCURRENT_EVNT.event_type;
    int idx, range = gCURRENT_EVNT.range;
    double time = MP_TIMER() - gCURRENT_EVNT.start_time;

    if(gCURRENT_EVNT.is_set) { /* Yep, there is an event set */
       GA_PROF[event_type][range].exectime += time;
       gCURRENT_EVNT.is_set = 0; /* clear the event */
    }
    else
       pnga_error("ENABLE_PROFILE: No event set. Probably ga_profile_stop() is called before ga_profile_start()", 0L);

#if GA_PRINT_STRIDE
    {  /* measure the time of each strided data transfer */
       idx = GA_PROF[event_type][range].count-1;
       if(idx<STRIDE_COUNT)  GA_PROF[event_type][range].stride[idx].time = time;
    }
#endif
}

#define GA_HDR1() printf("\n\n************ CONTIGUOUS DATA TRANSFER ************\n\n");
#define GA_HDR2() printf("\n\n********** NON-CONTIGUOUS DATA TRANSFER **********\n\n"); 
#define GA_HDR3() printf("RANK\t #Gets\t #puts\t #accs\t RANGE\n\n");
#define GA_HDR4() printf("RANK\t get_time\t put_time\t acc_time\t RANGE\n\n");
#define GA_HDR5() printf("SL#\tndim  time      stride_info (array name)\n\n");

/* This prints the number of contiguous get/put/acc/ calls for every 
   message range */
void ga_print_numcalls1() {
    int i; 
    GA_HDR1(); GA_HDR3();
    for(i=0; i< GA_MAX_MSG_RANGE-1; i++)
       printf("%d\t %d\t %d\t %d\t (%d-%d)\n", pnga_nodeid(),
	      GA_PROF[GET][i].count, GA_PROF[PUT][i].count,
	      GA_PROF[ACC][i].count, 1<<i, 1<<(i+1));
    printf("%d\t %d\t %d\t %d\t (>%d)\n", pnga_nodeid(),
	   GA_PROF[GET][i].count, GA_PROF[PUT][i].count,
	   GA_PROF[ACC][i].count, 1<<GA_MAX_MSG_RANGE);
}

/* This prints the number of non-contiguous get/put/acc/ calls for every 
   message range */
void ga_print_numcalls2() {
    int i; 
    GA_HDR2(); GA_HDR3();
    for(i=0; i< GA_MAX_MSG_RANGE-1; i++)
       printf("%d\t %d\t %d\t %d\t (%d-%d)\n", pnga_nodeid(),
	      GA_PROF[NC_GET][i].count, GA_PROF[NC_PUT][i].count,
	      GA_PROF[NC_ACC][i].count, 1<<i, 1<<(i+1));
    printf("%d\t %d\t %d\t %d\t (>%d)\n",pnga_nodeid(),
	   GA_PROF[NC_GET][i].count, GA_PROF[NC_PUT][i].count,
	   GA_PROF[NC_ACC][i].count, 1<<GA_MAX_MSG_RANGE);
}

/* This prints timings of all contiguous get/put/acc/ calls for every 
   message range */
void ga_print_timings1() {
    int i; 
    GA_HDR1(); GA_HDR4();
    for(i=0; i< GA_MAX_MSG_RANGE-1; i++)
       printf("%d\t %.2e\t %.2e\t %.2e\t (%d-%d)\n", pnga_nodeid(), 
	      GA_PROF[GET][i].exectime, GA_PROF[PUT][i].exectime, 
	      GA_PROF[ACC][i].exectime, 1<<i, 1<<(i+1));
    printf("%d\t %.2e\t %.2e\t %.2e\t (>%d)\n", pnga_nodeid(), 
	   GA_PROF[GET][i].exectime, GA_PROF[PUT][i].exectime, 
	   GA_PROF[ACC][i].exectime, 1<<GA_MAX_MSG_RANGE);    
}

/* This prints timings of all non-contiguous get/put/acc/ calls for every 
   message range */
void ga_print_timings2() {
    int i; 
    GA_HDR2(); GA_HDR4();
    for(i=0; i< GA_MAX_MSG_RANGE-1; i++)
       printf("%d\t %.2e\t %.2e\t %.2e\t (%d-%d)\n", pnga_nodeid(), 
	      GA_PROF[NC_GET][i].exectime, GA_PROF[NC_PUT][i].exectime, 
	      GA_PROF[NC_ACC][i].exectime, 1<<i, 1<<(i+1));
    printf("%d\t %.2e\t %.2e\t %.2e\t (>%d)\n", pnga_nodeid(), 
	   GA_PROF[NC_GET][i].exectime, GA_PROF[NC_PUT][i].exectime, 
	   GA_PROF[NC_ACC][i].exectime, 1<<GA_MAX_MSG_RANGE);    
}

void ga_print_stridedinfo(int event, int range) {
    int i, j, ndim;
    double time=0.0;
    printf("\n\nSTRIDE INFORMATION FOR MSG_RANGE %d-%d for EVENT: %s\n", 
	   1<<range, 1<<(range+1), event_name[event]);
    GA_HDR5();
    for(i=0; i< GA_PROF[event][range].count; i++) {
       if(i>=STRIDE_COUNT) break;
       time += GA_PROF[event][range].stride[i].time;
       ndim  = GA_PROF[event][range].stride[i].ndim;
       printf("%d\t%d     %.2e  (",i, ndim,
	      GA_PROF[event][range].stride[i].time);
       for(j=0;j<ndim;j++) {
	  printf("%d", GA_PROF[event][range].stride[i].hi[j] -
		       GA_PROF[event][range].stride[i].lo[j] +1);
	  if(j!=ndim-1) printf("x");
       }
       printf(") ==> ");
       for(j=0;j<ndim;j++)
	  printf("[%d-%d]", GA_PROF[event][range].stride[i].lo[j],
  		            GA_PROF[event][range].stride[i].hi[j]);
       printf(" \"%s\"\n",  GA_PROF[event][range].stride[i].name);
    }
    /*This o/p is just for verification*/
    printf("**** STRIDE_COUNT = %d ; TOTAL TIME = %.2e\n",GA_PROF[event][range].count,
	   time);
}

void ga_profile_terminate() {
    int i; 
    if(pnga_nodeid() == 0) { /* process 0's profile only */

       /* contiguous calls */
       ga_print_numcalls1();
       ga_print_timings1();

       /* non-contiguous calls */
       ga_print_numcalls2();
       ga_print_timings2();
       
#if GA_PRINT_STRIDE
       {
	  int msg_range, event_type;
	  /**
	   * printing stride info for non-contiguous get (NC_GET) for message 
	   * range #6.  2^6 - 2^(6+1) bytes. (i.e. 64-128 bytes)
	   */
	  msg_range = 6; /* message range 2^6-2^(6+1) */
	  event_type = NC_GET;
	  ga_print_stridedinfo(NC_GET, msg_range);
	  /*ga_print_stridedinfo(GET,19);*/ /* 2^19-2^20 range (524288-1MB)*/
       }
#endif
    }
}

#endif /* end of ENABLE_PROFILE */

