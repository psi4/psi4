#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      double precision function timer_dontuse()
c
c     return the time since the last call to timer.
c
c     must be initialized by calling once and throwing away the
c     value 
c     ... use cpu time on multi-user machines
c     ... use elapsed time on dedicated or single user machines.
c
*mdc*if unix
*      real*4 dtime, tt(2)
*      timer = dble(dtime(tt))
*mdc*elseif tcgmsg
C$Id: timer.F,v 1.1 2005-03-08 23:58:03 d3g293 Exp $
#if 0
#include "mpif.h"
      save mlast
      data mlast/0/
      m = MPI_WTIME()
      timer = dble(m - mlast) * 0.01d0
      mlast = m
*mdc*endif
c
#endif
      timer_dontuse = 0.0d0
      end
