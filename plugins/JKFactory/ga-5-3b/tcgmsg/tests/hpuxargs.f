c
c     $Header: /tmp/hpctools/ga/tcgmsg-mpi/hpuxargs.f,v 1.2 1999-06-08 21:08:29 d3h325 Exp $
c
      integer function hpargc()
      hpargc = iargc() + 1
      end
      integer function hpargv(index, arg, maxlen)
      character*256 arg
      hpargv = igetarg(index,arg,maxlen)
      end
