/*!
** \file
** \brief Obtain user and system timings for blocks of code
** \ingroup QT
**
** TIMER.CC: These functions allow one to obtain user and system
** timings for arbitrary blocks of code.  If a code block is called
** repeatedly during the course of program execution, the timer
** functions will report the block's cumulative execution time and
** the number of calls. In addition, one may time multiple code blocks
** simultaneously, and even ``overlap'' timers.  Timing data is
** written to the file "timer.dat" at the end of timer execution,
** i.e., when timer_done() is called.
**
** To use the timer functions defined here:
**
** (1) Initialize the linked list of timers at the beginning of your
** program: timer_init();
**
** (2) Start a timer at the start of the block of code:
** timer_on("My Timer");
**
** (3) Stop the timer at the end of the block: timer_off("My Timer");
**
** (4) When all timer calls are complete, dump the linked list of
** timing data to the output file, "timer.dat": timer_done();
**
** NB this code uses system functions ctime(), time(), and times(),
** which may not quite be standard on all machines.
**
** T. Daniel Crawford, August 1999.
*/

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <cstring>
#include <ctime>
#include <sys/param.h>
#include <sys/times.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <psi4-dec.h>

/* guess for HZ, if missing */
#ifndef HZ
#define HZ 60
#endif

#define TIMER_KEYLEN 128
#define TIMER_OFF 0
#define TIMER_ON 1

namespace psi {

struct timer {
    char key[TIMER_KEYLEN];
    unsigned int status;
    unsigned int calls;
    double utime;
    double stime;
    double wtime;
    struct tms ontime;
    time_t wall_start;
    struct timer *next;
    struct timer *last;
};

struct timer *global_timer;
time_t timer_start, timer_end;  /* Global wall-clock on and off times */

/*!
** timer_init(): Initialize the linked list of timers
**
** \ingroup QT
*/
void timer_init(void)
{
  extern struct timer *global_timer;
  extern time_t timer_start;

  timer_start = time(NULL);

  global_timer = NULL;
}

/*!
** timer_done(): Close down all timers and write results to timer.dat
**
** \ingroup QT
*/
void timer_done(void)
{
  FILE *timer_out;
  char *host;
  extern time_t timer_start, timer_end;
  struct timer *this_timer, *next_timer;
  extern struct timer *global_timer;

  timer_end = time(NULL);

  host = (char *) malloc(40 * sizeof(char));
  gethostname(host, 40);

  /* Dump the timing data to timer.dat and free the timers */
  ffile(&timer_out, "timer.dat", 1);
  fprintf(timer_out, "\n");
  fprintf(timer_out, "Host: %s\n", host);
  fprintf(timer_out, "\n");
  fprintf(timer_out, "Timers On : %s", ctime(&timer_start));
  fprintf(timer_out, "Timers Off: %s", ctime(&timer_end));
  fprintf(timer_out, "\nWall Time:  %10.2f seconds\n\n",
          (double) timer_end - timer_start);

  this_timer = global_timer;
  while(this_timer != NULL) {
      if(this_timer->calls > 1)
          fprintf(timer_out, "%-12s: %10.2fu %10.2fs %10.2fw %6d calls\n",
                  this_timer->key, this_timer->utime, this_timer->stime,
                  this_timer->wtime, this_timer->calls);
      else if(this_timer->calls == 1)
          fprintf(timer_out, "%-12s: %10.2fu %10.2fs %10.2fw %6d call\n",
                  this_timer->key, this_timer->utime, this_timer->stime,
                  this_timer->wtime, this_timer->calls);
      next_timer = this_timer->next;
      free(this_timer);
      this_timer = next_timer;
    }

  fprintf(timer_out,
          "\n***********************************************************\n");
  fclose(timer_out);

  free(host);

  global_timer = NULL;
}

/*!
** timer_scan(): Return a timer structure whose name matches that given
**   by supplied string
**
** \param key = name of timer to search for
**
** Returns: the timer structure with the given name, else NULL
** \ingroup QT
*/
struct timer *timer_scan(const char *key)
{
  extern struct timer *global_timer;
  struct timer *this_timer;

  this_timer = global_timer;

  while(this_timer != NULL) {
      if(!strcmp(this_timer->key,key)) return(this_timer);
      this_timer = this_timer->next;
    }

  return(this_timer);
}

/*
** timer_last(): Find the last timer in the list and return a pointer to it
**
** Returns: pointer to last timer in list, else NULL
**
** \ingroup QT
*/
struct timer *timer_last(void)
{
  extern struct timer *global_timer;
  struct timer *this_timer;

  this_timer = global_timer;

  while(this_timer != NULL) {
      if(this_timer->next == NULL) return(this_timer);
      this_timer = this_timer->next;
    }
  return(NULL);
}

/*!
** timer_on(): Turn on the timer with the name given as an argument.  Can
** be turned on and off, time will accumulate while on.
**
** \param key = Name of timer
**
** \ingroup QT
*/
void timer_on(const char *key)
{
  struct timer *this_timer;

  this_timer = timer_scan(key);

  if(this_timer == NULL) { /* New timer */
      this_timer = (struct timer *) malloc(sizeof(struct timer));
      strcpy(this_timer->key,key);
      this_timer->calls = 0;
      this_timer->utime = 0;
      this_timer->stime = 0;
      this_timer->wtime = 0;
      this_timer->next = NULL;
      this_timer->last = timer_last();
      if(this_timer->last != NULL) this_timer->last->next = this_timer;
      else global_timer = this_timer;
    }
  else {
    if((this_timer->status == TIMER_ON) && (this_timer->calls)) {
        std::string str = "Timer ";
        str += key;
        str += " is already on.";
        throw PsiException(str,__FILE__,__LINE__);
    }
  }

  this_timer->status = TIMER_ON;
  this_timer->calls++;

  times(&(this_timer->ontime));
  this_timer->wall_start = time(NULL);
}

/*!
** timer_off(): Turn off the timer with the name given as an argument.  Can
** be turned on and off, time will accumulate while on.
**
** \param key = Name of timer
**
** \ingroup QT
*/
void timer_off(const char *key)
{
  struct tms ontime, offtime;
  struct timer *this_timer;
  time_t wall_stop;

  this_timer = timer_scan(key);

  if(this_timer == NULL) {
      std::string str = "Bad timer key:";
      str += key;
      throw PsiException(str,__FILE__,__LINE__);
    }

  if(this_timer->status == TIMER_OFF) {
     std::string str = "Timer ";
     str += key;
     str += " is already off.";
     throw PsiException(str,__FILE__,__LINE__);
    }

  ontime = this_timer->ontime;

  times(&offtime);

  this_timer->utime += ((double) (offtime.tms_utime-ontime.tms_utime))/HZ;
  this_timer->stime += ((double) (offtime.tms_stime-ontime.tms_stime))/HZ;

  wall_stop = time(NULL);
  this_timer->wtime += ((double) (wall_stop - this_timer->wall_start));

  this_timer->status = TIMER_OFF;
}

}

