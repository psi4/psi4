*
* $Id: fci_tinit.F,v 1.3 1999/07/28 00:23:42 d3e129 Exp $
*
      subroutine fci_tinit
      implicit none
#include "fcitime.fh"      
c
      integer i
c
      do i = 1, maxroutine
         names(i) = ' '
         times(i) = 0.0d0
      enddo
c
      names(1)  = 'couple1'
      names(2)  = 'couple2'
      names(3)  = 'couple3'
      names(4)  = 'davidson'
      names(5)  = 'diags'
      names(6)  = 'make_eaj '
      names(7)  = 'setup'
      names(8)  = 'sigma'
      names(9)  = 'solve'
      names(10) = 'spnfrc'
      names(11) = 'cntone'
      names(12) = 'hami1'
      names(13) = 'hami2'
      names(14) = 'symfrc'
c
      end
      subroutine fci_tstart(i,oprint)
      implicit none
#include "fcitime.fh"      
      integer i
      logical oprint
      integer fci_strlen
      external fci_strlen
c
      if (i.le.0  .or.  i.gt.maxroutine) call fci_err('tstart?',i)
      ncalls(i) = ncalls(i) + 1
      if (ncalls(i) .lt. maxcalls) then   ! Avoid too many system calls
        started(i) = fci_time()
      else
        started(i) = 0.0d0
      endif
c     
      if (oprint) write(6,1) names(i)(1:fci_strlen(names(i))), 
     $     ncalls(i), started(i)
 1    format('  Started ',a,': call number',f10.0,' at', f9.1,'s')
c
      end
      subroutine fci_tfinish(i,oprint)
      implicit none
#include "fcitime.fh"      
      integer i
      logical oprint
      double precision now
      integer fci_strlen
      external fci_strlen
c
      if (i.le.0  .or.  i.gt.maxroutine) call fci_err('tfinish?',i)
c
      if (ncalls(i) .lt. maxcalls) then  ! Avoid too many system calls
        now = fci_time()
        times(i) = times(i) + now - started(i)
      else
        now = 0.0d0
      endif
c
      if (oprint) write(6,1) names(i)(1:fci_strlen(names(i))), 
     $     ncalls(i), now, now-started(i),
     $     times(i)
 1    format(' Finished ',a,': call number',f10.0,' at', f9.1,'s',
     $     ' used ',f9.1,'s', ' total ',f9.2, 's')
c
      end
      subroutine fci_tprint()
      implicit none
#include "fcitime.fh"      
      integer i
      double precision used, now, totalcalls
c
      write(6,1) 
 1    format(/'  FCI timing information '/
     $        '  ----------------------'/
     $     1x, ' name                    ncalls       time'/
     $     1x, ' --------------------  ----------  ----------')
      totalcalls = 0
      do i = 1, maxroutine
         if (ncalls(i) .gt. 0) then
            if (ncalls(i) .gt. maxcalls) times(i) = 0.0d0  ! Did not measure it
            write(6,2) names(i), ncalls(i), times(i)
 2          format(2x,a20,2x,f10.0,2x,f9.1)
            if (times(i).gt.0.0d0) totalcalls = totalcalls + 2*ncalls(i)
         endif
      enddo
c
c     Determine the overhead due to timing
c
      used = -fci_time()
      do i = 1, maxcalls
         now = fci_time()
         if (used+now .gt. 2) goto 10
      enddo
 10   used = (used+now)/dfloat(i)
      write(6,3) used, totalcalls, used*totalcalls
 3    format(/1p,'  Timing overhead=',d8.1,'/call, #calls=',
     $     d8.1, ', total ovehead=',d8.1/)
c
      end
      double precision function fci_time()
      implicit none
c
c     return the CPU time in seconds from the first call to this routine
c
c     MACHINE DEPENDENT
c
      double precision firsttime
      logical firstcall
#ifdef AIX
      integer mclock
      data firstcall /.true./, firsttime/0.0d0/
c
      if (firstcall) then
         firsttime = dfloat(mclock())*0.01d0
         firstcall = .false.
      endif
      fci_time = dfloat(mclock())*0.01d0 - firsttime
#elif defined(LINUX)
      double precision linux_cputime
      external linux_cputime
      data firstcall /.true./, firsttime/0.0d0/
      if (firstcall) then
         firsttime = linux_cputime()
         firstcall = .false.
      endif
c
      fci_time = linux_cputime() - firsttime
#else
      real*4 tarray(2)
      real*4 etime
      external etime
      data firstcall /.true./, firsttime/0.0d0/
c
      if (firstcall) then
         firsttime = etime(tarray)
         firstcall = .false.
      endif
c
      fci_time = etime(tarray) - firsttime
#endif
c
      end
      integer function fci_strlen(a)
      implicit none
      character*(*) a
      integer i
      intrinsic len
c
      do i = len(a),1,-1
         if (a(i:i) .ne. ' ') goto 10
      enddo
c
 10   fci_strlen = i
c
      end
