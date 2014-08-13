#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c $Id: ffflush.F,v 1.2 2002-10-07 16:50:05 d3h325 Exp $
      subroutine ffflush(unit)
      integer unit
c
#ifdef CRAY
*     if(unit.eq.6)then
*       call flush(101)
*     else
*       call flush(unit)
*     endif
#elif HAVE_F77_FLUSH
      call F77_FLUSH(unit)
#endif
c
      end
