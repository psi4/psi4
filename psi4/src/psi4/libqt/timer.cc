/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

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
**
** Modified to use timeval structures for the module wall times, getting
** nanosecond precision on cumulated wall time instead of second.
** Useful to time integral computations where there can be millions
** to billion calls to functions.
**
** J. F. Gonthier, February 2016
*/

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <cstring>
#include <ctime>
#include <sys/param.h>
#include <sys/times.h>
#include "psi4/libciomr/libciomr.h"
#include "psi4/psifiles.h"
#include "psi4/psi4-dec.h"
#include "psi4/libparallel/ParallelPrinter.h"
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
    struct timeval wall_start;
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

  char *host;
  extern time_t timer_start, timer_end;
  struct timer *this_timer, *next_timer;
  extern struct timer *global_timer;

  timer_end = time(NULL);

  host = (char *) malloc(40 * sizeof(char));
  gethostname(host, 40);

  /* Dump the timing data to timer.dat and free the timers */
  std::shared_ptr<OutFile> printer(new OutFile("timer.dat",APPEND));
  printer->Printf( "\n");
  printer->Printf( "Host: %s\n", host);
  printer->Printf( "\n");
  printer->Printf( "Timers On : %s", ctime(&timer_start));
  printer->Printf( "Timers Off: %s", ctime(&timer_end));
  printer->Printf( "\nWall Time:  %10.2f seconds\n\n",
          (double) timer_end - timer_start);

  this_timer = global_timer;
  while(this_timer != NULL) {
      if(this_timer->calls > 1) {
          if(this_timer->wtime < 10.0) {
              printer->Printf( "%-12s: %10.2fu %10.2fs %10.6fw %6d calls\n",
                      this_timer->key, this_timer->utime, this_timer->stime,
                      this_timer->wtime, this_timer->calls);
          } else {
              printer->Printf( "%-12s: %10.2fu %10.2fs %10.2fw %6d calls\n",
                      this_timer->key, this_timer->utime, this_timer->stime,
                      this_timer->wtime, this_timer->calls);
          }
      }
      else if(this_timer->calls == 1) {
          if(this_timer->wtime < 10.0) {
              printer->Printf( "%-12s: %10.2fu %10.2fs %10.8fw %6d call\n",
                      this_timer->key, this_timer->utime, this_timer->stime,
                      this_timer->wtime, this_timer->calls);
          } else {
              printer->Printf( "%-12s: %10.2fu %10.2fs %10.2fw %6d call\n",
                      this_timer->key, this_timer->utime, this_timer->stime,
                      this_timer->wtime, this_timer->calls);
          }
      }
      next_timer = this_timer->next;
      free(this_timer);
      this_timer = next_timer;
    }

  printer->Printf(
          "\n***********************************************************\n");

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
  gettimeofday(&(this_timer->wall_start), NULL);
}

/*!
** timer_nsdiff(struct timeval endt, struct timeval begint): Returns
** the time difference between two timeval values as a double, in seconds.
** Takes into account the difference in nanoseconds accurately.
** Taken from gnu.org, supposedly the best way to do that works even on
** peculiar OS where tv_sec is unsigned.
**
** \param endt = Final timeval value as obtained through gettimeofday()
**
** \param begint = Initial timeval value as obtained through gettimeofday()
**
** Returns: the time difference in seconds as a double.
**
** \ingroup QT
*/

double timer_nsdiff(struct timeval& endt, struct timeval& begint) {
  /*Carry the difference in nanoseconds to the second if needed */
  long int nsec;
  long int million = 1000 * 1000;
  if(endt.tv_usec < begint.tv_usec) {
    nsec = (begint.tv_usec - endt.tv_usec) / million + 1;
    begint.tv_usec -= million * nsec;
    begint.tv_sec += nsec;
  }
  if(endt.tv_usec - begint.tv_usec > million) {
    nsec = (endt.tv_usec - begint.tv_usec) / million;
    begint.tv_usec += million * nsec;
    begint.tv_sec -= nsec;
  }

  return ( (double)(endt.tv_sec - begint.tv_sec) + ((endt.tv_usec - begint.tv_usec) / ((double)million)) );
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
  struct timeval wall_stop;

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

  gettimeofday(&wall_stop, NULL);
  this_timer->wtime += timer_nsdiff(wall_stop, this_timer->wall_start);

  this_timer->status = TIMER_OFF;
}

}
