#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      program main
      implicit double precision (a-h,o-z)
c
c $Header: /tmp/hpctools/ga/tcgmsg/ipcv5.0/testf.f,v 1.2 2001-05-08 17:30:39 edo Exp $
c
c     FORTRAN program to test message passing routines
c
c     LOG is the FORTRAN unit number for standard output.
c
      parameter (LOG = 6)
      parameter (MAXLEN = 262144 / 8)
#include "msgtypesf.h"
      dimension buf(MAXLEN)
      integer ibuf(MAXLEN)
      character*80 fortnn
      character*80 fname
c
c     Always the first thing to do is call pbeginf
c
      call pbeginf
      call setdbg(0)
c
c     who am i and how many processes
c
      nproc = nnodes()
      me = nodeid()
c
c     now try broadcasting messages from all nodes to every other node
c     send each process my id as a message
c     
      itype = 1 + MSGINT
      do 10 iproc = 0,nproc-1
         itest = me
         call brdcst(itype, itest, mitob(1),iproc)
         if (iproc.ne.me) then
            write(LOG,1) me, itest
 1          format(' me=',i3,', itest=',i3)
         endif
 10   continue
c
c     now try using the shared counter
c
      mproc = nproc
      do 20 i = 1,40
         write(LOG,*) ' process ',me,' got ',nxtval(mproc)
 20   continue
      junk = nxtval(-mproc)
c
c     now time sending a message round a ring
c
      if (nproc.gt.1) then
        itype = 3
        left = mod(me + nproc - 1, nproc)
        iright = mod(me + 1, nproc)
c      
        lenbuf = 1
 30     if (me .eq. 0) then
           istart = mtime()
           call snd(itype, buf, lenbuf, left, 1)
           call rcv(itype, buf, lenbuf, lenmes, iright, node, 1)
           iused = mtime() - istart
           if (iused.gt.0) then
             rate = 1.0d-4 * dble(nproc * lenbuf) / dble(iused)
           else
             rate = 0.0d0
           endif
           write(LOG,31) lenbuf, iused, rate
        else
           call rcv(itype, buf, lenbuf, lenmes, iright, node, 1)
           call snd(itype, buf, lenbuf, left, 1)
        endif
        lenbuf = lenbuf * 2
        if (lenbuf .le. mdtob(MAXLEN)) goto 30
 31     format(' len=',i7,' bytes, used=',i4,' cs, rate=',f10.6,' Mb/s')
      endif
c
c     global sums
c
      do i=1,MAXLEN
         ibuf(i) = i*me
         buf(i) = dble(ibuf(i))
      enddo
      dtype=1+MSGDBL
      call igop(itype, ibuf, MAXLEN, "+")
      call dgop(dtype, buf, MAXLEN, "+")
      
      do i=1,MAXLEN
         iresult = i*nproc*(nproc-1)/2
         if (ibuf(i).ne.iresult.or.buf(i).ne.dble(iresult))
     .      call error('TestGlobals: global sum failed',  i)
      enddo
      
      if (me.eq.0) write(LOG,*) 'global sums OK'
c
c
c     Check that everyone can open, write, read and close
c     a binary FORTRAN file
c
      do ii=29,45
      write(fortnn,1234) ii
1234  format('fort',i2)
      call pfname(fortnn,fname)
      open(ii,file=fname,form='unformatted',status='unknown',
     &  err=1000)
      enddo
      do ii=11,19
      write(ii,err=1001) buf
      rewind ii
      read(ii,err=1002) buf
      close(ii,status='delete')
      enddo
c
      if (me.eq.0) call stats
c
c     Always the last thing to do is call pend
c
      call pend
c
c     check that everyone makes it thru after pend .. NODEID
c     is not actually guaranteed to work outside of pbegin/pend
c     section ... it may return junk. All you should do is exit
c     is some FORTRAN supported fashion
c
      write(LOG,32) nodeid()
 32   format(' Process ',i4,' after pend')
      stop
c
c     error returns for FORTRAN I/O
c
1000  call error('failed to open fortran binary file',-1)
1001  call error('failed to write fortran binary file',-1)
1002  call error('failed to read fortran binary file',-1)
c
      end
      subroutine pfname(name, fname)
      character*(*) name, fname
c
c     construct a unique filename by appending the process
c     number after the stub name
c     i.e. <fname> = <name>.<mynode>
c
c     find last non-blank character in name
c
      do 10 i = len(name),1,-1
      if (name(i:i).ne.' ') goto 20
10    continue
      call error('pfname: name is all blanks!',i)
c
c     check that have room to store result and then write result
c
20    if (i+4.gt.len(fname))
     &  call error('pfname: fname too short for name.id',len(fname))
      fname = name
      write(fname(i+1:i+4),1) nodeid()
1     format('.',i3.3)
c
      end
      subroutine error(s,i)
      parameter (LOG = 6)
      character*(*) s
      integer i
c
      write(LOG,1) s,i
 1    format(//
     $     ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'/
     $     1x,a,1x,i8/
     $     ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'/)
c    $     1x,a,1x,i8/
c
      call parerr(i)
c
      end
