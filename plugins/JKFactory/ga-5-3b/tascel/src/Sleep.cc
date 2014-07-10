#if HAVE_CONFIG_H
# include "config.h"
#endif

#ifdef _POSIX_C_SOURCE
# undef _POSIX_C_SOURCE
# define _POSIX_C_SOURCE 199309
#endif

#ifdef _XOPEN_SOURCE
# undef _XOPEN_SOURCE
# define _XOPEN_SOURCE 500
#endif

#include <time.h>
#include <unistd.h>

#include "Sleep.h"

void tascel::sleep(long seconds) {
  unsigned int s = static_cast<unsigned int>(seconds);
  (void)sleep(s);
}

void tascel::microsleep(long microseconds) {
  (void)usleep(microseconds);
}

void tascel::nanosleep(long nanoseconds) {
  int ret;
  struct timespec ts;
  ts.tv_sec = 0;
  ts.tv_nsec = nanoseconds;
  ret = nanosleep(&ts, NULL);
}
