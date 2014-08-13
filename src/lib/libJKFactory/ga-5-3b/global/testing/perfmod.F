#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c $Id: perfmod.F,v 1.4 2000-05-25 01:09:20 d3h325 Exp $
c------------------------------------------------------------------------
c Program perfmod.x is used to test performance of GA put,get,accumulate |
c It has to be executed on four processors.                              |
c remote operations access data on processes 1,2,3 in the round-robin way|
c------------------------------------------------------------------------
c
#define SLEEP sleep
      program perfmod
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer heap,stack
c
c***  Intitialize a message passing library
c
#include "mp3.fh"
c
c     Intitialize the GA package
c
      call ga_initialize()
c
      if(ga_nnodes().ne.4 .and. ga_nodeid().eq.0)
     $   call ga_error('Must be run with 4 GA processes',ga_nnodes())
c
c***  Initialize the MA package
      heap = 450000
      stack =heap
      if (.not. ma_init(MT_DBL, heap, stack))
     $    call ga_error("ma init failed",heap+stack)
c
      call testit()
c
      if(ga_nodeid().eq.0) print *, 'All tests successful'
c
      call ga_terminate()
c
      call MP_FINALIZE()
      end


      subroutine testit()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
c     
c
      integer n, nn, num_chunks
      parameter (n = 710, nn = n*n/4, num_chunks=12)
      double precision buf(nn)
c
      integer g_a
      integer ilo, ihi, jlo, jhi
      integer nproc, me, loop
      integer chunk(num_chunks)
      data    chunk /1,3,4,9,16,30,64,91,128,171,256,353/
c     
      nproc = ga_nnodes()
      me = ga_nodeid()
c
c***  Create global array
      if (.not. ga_create(MT_DBL, n, n, 'a', 0, 0, g_a))
     $     call ga_error(' ga_create failed ',n)
c     
      do loop=1,nn
         buf(loop) = .01d0
      enddo
      call ga_zero(g_a) 
c
      if (me .eq. 0) then
        write(*,*)' '
        print *,'> Performance of ga_get, ga_put & ga_acc n = ', n
        print *,' '
      endif
c
c     do loop=1,2
c
c***  local ops
c
      call ga_distribution(g_a, me, ilo, ihi, jlo, jhi)
      call TestPutGetAcc
     &     (g_a, n, chunk, num_chunks, buf, ilo, ihi, jlo, jhi, .true.)
c
c***  remote ops
c
      call TestPutGetAcc
     &     (g_a, n, chunk, num_chunks, buf, ilo, ihi, jlo, jhi,.false.)

c     enddo
      end


      subroutine TestPutGetAcc
     &      (g_a, n, chunk, num_chunks, buf, ilo, ihi, jlo,jhi, local)
      implicit none
#include "global.fh"
#include "testutil.fh"
      integer num_chunks, chunk(num_chunks)
      integer n, ilo, ihi, jlo,jhi,g_a
      double precision buf(*), tg, tp, ta
      double precision time_acc, time_get, time_put
      logical local
c
      integer me
      integer loop, jump, count, bytes
c
      me = ga_nodeid()
      if (me .eq. 0) then
        write(6,*)' '
        if(local) then
          write(6,'(21X,8hACCESS [,i3,1h:,i4,1h,,i3,1h:,i4,1h])') 
     &        ilo,ihi,jlo,jhi
        else
          write(6,'(21X,6hACCESS , 2x, 18Hremote section       )') 
        endif

        write(6,*)'bytes  loop         get                    put',
     &           '                 accumulate'
        call flush(6)
      endif
      call ga_sync()
c
      do loop = 1, num_chunks
        bytes = util_mdtob(1)*chunk(loop)*chunk(loop) !how much data is accessed
        jump  =  n/(60*loop) ! jump distance between consecutive patches
c
c       everybody touches own data 
        call ga_fill_patch(g_a, 1, n, 1, n , 1d0*me*loop)
        if (me .eq. 0) then
        tg=time_get(g_a,ilo,ihi,jlo,jhi,buf,chunk(loop),jump,count,
     $               local)
        else
          call SLEEP(2)
        endif
c
c       everybody touches own data
        call ga_fill_patch(g_a, 1, n, 1, n , 1d0*me*loop)
        if (me .eq. 0) then
        tp=time_put(g_a,ilo,ihi,jlo,jhi,buf,chunk(loop),jump,count, 
     $               local)
        else
          call SLEEP(2)
        endif
c
c       everybody touches own data
        call ga_fill_patch(g_a, 1, n, 1, n , 1d0*me*loop)
        if (me .eq. 0) then
        ta=time_acc(g_a,ilo,ihi,jlo,jhi,buf,chunk(loop),jump,count,
     $               local)
        else
          call SLEEP(2)
        endif
c
        if (me .eq. 0) then
          write(6,77)bytes, count, tg, 1d-6*bytes/tg,
     &               tp, 1d-6*bytes/tp, ta, 1d-6*bytes/ta
          call flush(6)
        endif
      enddo
c
77    format(i6, i5, 3(d10.3, d9.2,'MB/s'))
      end



      double precision function
     &   time_acc(g_a, is, ie, js, je, buf, chunk, jump, count, local)
c
      implicit none
#include "global.fh"
#include "testutil.fh"
c
      integer g_a, chunk, jump, count, is, js, ie, je
      logical local 
      integer rows, cols, indx, shifti(3), shiftj(3)
c
      integer ilo, ihi, jlo, jhi
      double precision  seconds, buf 
c
      count = 0
      rows = ie - is + 1
      cols = je - js + 1
      shifti(1) = rows
      shifti(2) = 0
      shifti(3) = rows
      shiftj(1) = 0
      shiftj(2) = cols
      shiftj(3) = cols

      seconds = util_timer()
c
c       distance between consecutive patches increased by jump
c       to destroy locality of reference
        do ilo = is, ie -chunk-jump +1, chunk+jump
           ihi = ilo + chunk -1
           do jlo = js, je -chunk-jump +1, chunk+jump
              jhi = jlo + chunk -1
              count = count + 1
              if (local) then
                 call ga_acc(g_a, ilo, ihi, jlo, jhi, buf, chunk, 1d0)
              else
                 indx = Mod(count,3) + 1 
                 call ga_acc(g_a, ilo+shifti(indx), ihi+shifti(indx),
     $                       jlo+shiftj(indx), jhi+shiftj(indx), 
     $                       buf, chunk, 1d0)
              endif
           enddo
        enddo
      seconds = util_timer() - seconds
c
      time_acc = seconds/count
      end


      double precision function
     &    time_get(g_a, is, ie, js, je, buf, chunk, jump, count, local)
c
      implicit none
#include "global.fh"
#include "testutil.fh"
c
      integer g_a, chunk, jump, count, is, js, ie, je
      integer rows, cols, indx, shifti(3), shiftj(3)
      logical local
c
      integer ilo, ihi, jlo, jhi
      double precision  seconds, buf
c
      count = 0
      rows = ie - is + 1
      cols = je - js + 1
      shifti(1) = rows
      shifti(2) = 0
      shifti(3) = rows
      shiftj(1) = 0
      shiftj(2) = cols
      shiftj(3) = cols

      seconds = util_timer()
c
c       distance between consecutive patches increased by jump
c       to destroy locality of reference
        do ilo = is, ie -chunk-jump +1, chunk+jump
           ihi = ilo + chunk -1
           do jlo = js, je -chunk-jump +1, chunk+jump
              jhi = jlo + chunk -1
              count = count + 1
              if (local) then
                 call ga_get(g_a, ilo, ihi, jlo, jhi, buf, chunk)
              else
                 indx = Mod(count,3) + 1
                 call ga_get(g_a, ilo+shifti(indx), ihi+shifti(indx),
     $                       jlo+shiftj(indx), jhi+shiftj(indx),
     $                       buf, chunk)
              endif
           enddo
        enddo
      seconds = util_timer() - seconds
c
      time_get = seconds/count
      end



      double precision function
     &   time_put(g_a, is, ie, js, je, buf, chunk, jump, count, local)
c
      implicit none
#include "global.fh"
#include "testutil.fh"
c
      integer g_a, chunk, jump, count, is, js, ie, je
      integer rows, cols, indx, shifti(3), shiftj(3)
      logical local
c
      integer ilo, ihi, jlo, jhi
      double precision  seconds, buf
c
      count = 0
      rows = ie - is + 1
      cols = je - js + 1
      shifti(1) = rows
      shifti(2) = 0
      shifti(3) = rows
      shiftj(1) = 0
      shiftj(2) = cols
      shiftj(3) = cols

      seconds = util_timer()
c
c       distance between consecutive patches increased by jump
c       to destroy locality of reference
        do ilo = is, ie -chunk-jump +1, chunk+jump
           ihi = ilo + chunk -1
           do jlo = js, je -chunk-jump +1, chunk+jump
              jhi = jlo + chunk -1
              count = count + 1
              if (local) then
                 call ga_put(g_a, ilo, ihi, jlo, jhi, buf, chunk)
              else
                 indx = Mod(count,3) + 1
                 call ga_put(g_a, ilo+shifti(indx), ihi+shifti(indx),
     $                       jlo+shiftj(indx), jhi+shiftj(indx),
     $                       buf, chunk)
              endif
           enddo
        enddo
      seconds = util_timer() - seconds
c
      time_put = seconds/count
      end


      subroutine my_sleep(seconds)
      implicit none
#include "global.fh"
#include "testutil.fh"
c
      double precision  ts, te, work
      integer seconds, loop
      common  /sleep_block/ work
c
      ts = util_timer()
      work = 0.
100   continue
         do loop= 1, 500
            work = work + 1.
         enddo
         te = util_timer()
      if(te - ts .lt. real(seconds)) goto 100
*     print *, work, ts, te 
      end
