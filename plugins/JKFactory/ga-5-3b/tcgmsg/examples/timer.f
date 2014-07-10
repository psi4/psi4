      double precision function timer()
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
C$Id: timer.f,v 1.2 1995-02-02 23:24:33 d3g681 Exp $
      save mlast
      data mlast/0/
      m = mtime()
      timer = dble(m - mlast) * 0.01d0
      mlast = m
*mdc*endif
c
      end
