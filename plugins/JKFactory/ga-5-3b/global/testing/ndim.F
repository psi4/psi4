#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c$Id: ndim.F,v 1.7 2005-11-23 10:25:18 manoj Exp $
#define  MAXLOOP 100

      program test
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      integer nproc
      logical status
c
c***  Intitialize a message passing library
c
#include "mp3.fh"
c     Intitialize the GA package
c
      call ga_initialize()
      nproc = ga_nnodes()
      if(ga_nodeid().eq.0)print *,nproc,' nodes'
c
c     Initialize the MA package
c
      status = ma_init(MT_DBL, 500000/nproc, 50000)
      if(.not. status) call ga_error("ma_init failed",0) 
c
c
      if(ga_nodeid().eq.0)print *,' Testing 3-dimensional Arrays'
      call testit()
      if(ga_nodeid().eq.0)print *,' Testing 4-dimensional Arrays'
      call testit4()
      if(ga_nodeid().eq.0)print *,' All tests successful'
      call ga_terminate()
c
c***  Tidy up after message-passing library
c
      call MP_FINALIZE()
      end


c-----------------



      subroutine testit()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c     
      integer n
      integer ndim
      parameter (n = 38)
      parameter (ndim = 3)
      double precision a(n,n,n),b(n,n,n)
      integer g_a
      integer  i,  lo(ndim),hi(ndim), lop(ndim),hip(ndim),elems
      integer nproc, me, proc, loop, maxloop
      integer chunk(ndim), dims(ndim), adims(ndim), ld(ndim)
      logical status, compare_patches
      integer count_elems
      double precision crap,alpha
c     
      nproc = ga_nnodes()
      me    = ga_nodeid()
c     
      
      call ifill_array(chunk,ndim,0)
      call ifill_array(adims,ndim,n-1)
      call ifill_array(dims,ndim,n)
      call ifill_array(ld,ndim,n)
      call dfill_array(a,n*n*n,dble(me))
      call dfill_array(b,n*n*n,-1d0)
c
c***  Create global arrays
      if (.not. nga_create(MT_DBL, ndim, adims, 'a', chunk, g_a))
     $     call ga_error(' ga_create failed ',1)
c     
      call ga_sync()
      if(me.eq.0)then
         write(6,'(i2,21H-dimensional Array A:  ,10i6)')
     $   ndim,(adims(i),i=1,ndim)
         print *,'distribution information for all processors'
         print *,'-------------------------------------------'
c        call ffflush(6)
      endif
      call ga_sync()
      call nga_distribution(g_a, me, lo,hi)
      elems = count_elems(lo,hi,ndim)
c
      do i = 0, nproc-1
       if (me .eq. i) then
100      format(i4,' has',i8,' elements of A, range:',10(i3,':',i3,','))
c        write(*,100)me,elems,(lo(j),hi(j),j=1,ndim)
         call print_range(me, lo, hi, ndim)
         call ffflush(6)
       endif
       call ga_sync()
      enddo
c
c------------------------------- GA_FILL ----------------------------
      call ga_fill(g_a,dble(me))
      if(me.eq.0)then 
        print *, ' '
        print *, 'Filling array A'
        call ffflush(6)
      endif
c     call ga_print(g_a)
      call ga_sync()
c
      if(elems.gt.0) then
        call nga_get(g_a,lo,hi,b(lo(1),lo(2),lo(3)),ld)

        if(compare_patches(me,0d0,ndim,a,lo,hi,dims,b,lo,hi,dims))
     $             call ga_error('bye',0)
      endif

      call ga_sync()
      if(me.eq.0)then 
        print *, 'OK'
        print *, ' '
        print *, 'Testing random PUT'
        print *,'(only process 0 prints range for its every 10-th put)'
        call ffflush(6)
      endif
      call ga_fill(g_a,-1d0)
c
c------------------------------- NGA_PUT ----------------------------
c     if(nproc.gt.0)return
      proc =  nproc-1 -me ! access other process memory 
      call nga_distribution(g_a, proc, lo,hi)
      elems = count_elems(lo,hi,ndim)
      call init_array(a,ndim,dims)
c
      call ga_sync()
      if(elems.gt.0) then
         call nga_put(g_a,lo,hi,a(lo(1),lo(2),lo(3)),ld)
         do loop = 1, MAXLOOP 
            call random_range(lo,hi,lop,hip,ndim)
            if(me.eq.0 .and. Mod(loop,10).eq.0)then
               call print_range(loop,lop,hip,ndim)
            endif
            call nga_put(g_a,lop,hip,a(lop(1),lop(2),lop(3)),ld)
         enddo

         call nga_get(g_a,lo,hi,b(lo(1),lo(2),lo(3)),ld)

         if(compare_patches(me,0d0,ndim,a,lo,hi,dims,b,lo,hi,dims))
     $             call ga_error('bye',0)

      endif
c
      call ga_sync()
      if(me.eq.0)then
        print *, 'OK'
        print *, ' '
        print *, 'Testing random GET'
        print *,'(only process 0 prints range for its every 10-th get)'
        call ffflush(6)
      endif
c------------------------------- NGA_GET ----------------------------
      call ga_sync()
      call ifill_array(lop,ndim,1)
      call ifill_array(hip,ndim,n-1)
      do loop = 1, MAXLOOP 
         call random_range(lop,hip,lo,hi,ndim)
         if(me.eq.0 .and. Mod(loop,10).eq.1)then
            call print_range(loop,lo,hi,ndim)
         endif
         call nga_get(g_a,lo,hi,b(lo(1),lo(2),lo(3)),ld)
         if(compare_patches(me,0d0,ndim,a,lo,hi,dims,b,lo,hi,dims))
     $             call ga_error('bye',0)
      enddo
c------------------------------- NGA_ACC ----------------------------
      call ga_sync()
      if(me.eq.0)then
        print *, 'OK'
        print *, ' '
        print *, 'Testing Accumulate'
        call ffflush(6)
      endif
c     
      call ga_sync()
      call ifill_array(lop,ndim,1)
      call ifill_array(hip,ndim,n-1)
      call random_range(lop,hip,lo,hi,ndim)
      crap = drand(1)
      maxloop = 10
      alpha = .1d0 ! alpha must be 1/maxloop
      call ga_sync()
c
      do loop=1, maxloop
       call nga_acc(g_a,lop,hip,a(lop(1),lop(2),lop(3)),ld,alpha)
      enddo
      call ga_sync()
      if(me.eq.0)then
        print *, 'multiple accumulate target same array section'
        call print_range(maxloop,lo,hi,ndim) 
        call nga_get(g_a,lo,hi,b(lo(1),lo(2),lo(3)),ld)
        call scale_patch(dble(nproc+1),ndim, a(lo(1),lo(2),lo(3)),
     $                    lo, hi, dims)
        if(compare_patches(me,1d-2,ndim,a,lo,hi,dims,b,lo,hi,dims))
     $             call ga_error('bye',0)
        print *, 'OK'
        call ffflush(6)
      endif
c     
      status= ga_destroy(g_a)
      end


      subroutine testit4()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c     
      integer n
      integer ndim
      parameter (n = 25)
      parameter (ndim = 4)
      double precision a(n,n,n,n),b(n,n,n,n)
      integer g_a
      integer  i, lo(ndim),hi(ndim), lop(ndim),hip(ndim),elems
      integer nproc, me, proc, loop, maxloop
      integer chunk(ndim), dims(ndim), ld(ndim)
      logical status, compare_patches
      integer count_elems
      double precision crap,alpha
c     
      nproc = ga_nnodes()
      me    = ga_nodeid()
c     
      
      call ifill_array(chunk,ndim,0)
      call ifill_array(dims,ndim,n)
      call ifill_array(ld,ndim,n)
      elems=1
      do i = 1,ndim
         elems = elems * dims(i)
      enddo
      call dfill_array(a,elems,dble(me))
      call dfill_array(b,elems,-1d0)
c
c***  Create global arrays
      if (.not. nga_create(MT_DBL, ndim, dims, 'a', chunk, g_a))
     $     call ga_error(' ga_create failed ',1)
c     
      call ga_sync()
      if(me.eq.0)then
         write(6,'(i2,21H-dimensional Array A:  ,10i6)')
     $   ndim,(dims(i),i=1,ndim)
         print *,'distribution information for all processors'
         print *,'-------------------------------------------'
c        call ffflush(6)
      endif
      call ga_sync()
      call nga_distribution(g_a, me, lo,hi)
      elems = count_elems(lo,hi,ndim)
c
      do i = 0, nproc-1
       if (me .eq. i) then
100      format(i4,' has',i8,' elements of A, range:',10(i3,':',i3,','))
c        write(*,100)me,elems,(lo(j),hi(j),j=1,ndim)
         call print_range(me, lo, hi, ndim)
         call ffflush(6)
       endif
       call ga_sync()
      enddo
c
c------------------------------- GA_FILL ----------------------------
      call ga_fill(g_a,dble(me))
      if(me.eq.0)then 
        print *, ' '
        print *, 'Filling array A'
        call ffflush(6)
      endif
c     call ga_print(g_a)
      call ga_sync()
c
      if(elems.gt.0) then
        call nga_get(g_a,lo,hi,b(lo(1),lo(2),lo(3),lo(4)),ld)

        if(compare_patches(me,0d0,ndim,a,lo,hi,dims,b,lo,hi,dims))
     $             call ga_error('bye',0)
      endif

      call ga_sync()
      if(me.eq.0)then 
        print *, 'OK'
        print *, ' '
        print *, 'Testing random PUT'
        print *,'(only process 0 prints range for its every 10-th put)'
        call ffflush(6)
      endif
      call ga_fill(g_a,-1d0)
c
c------------------------------- NGA_PUT ----------------------------
c     if(nproc.gt.0)return
      proc =  nproc-1 -me ! access other process memory 
      call nga_distribution(g_a, proc, lo,hi)
      elems = count_elems(lo,hi,ndim)
      call init_array(a,ndim,dims)
c
      call ga_sync()
      if(elems.gt.0) then
         call nga_put(g_a,lo,hi,a(lo(1),lo(2),lo(3),lo(4)),ld)
         do loop = 1, MAXLOOP 
            call random_range(lo,hi,lop,hip,ndim)
            if(me.eq.0 .and. Mod(loop,10).eq.0)then
               call print_range(loop,lop,hip,ndim)
            endif
            call nga_put(g_a,lop,hip,a(lop(1),lop(2),lop(3),lop(4)),ld)
         enddo

         call nga_get(g_a,lo,hi,b(lo(1),lo(2),lo(3),lo(4)),ld)

         if(compare_patches(me,0d0,ndim,a,lo,hi,dims,b,lo,hi,dims))
     $             call ga_error('bye',0)

      endif
c
      call ga_sync()
      if(me.eq.0)then
        print *, 'OK'
        print *, ' '
        print *, 'Testing random GET'
        print *,'(only process 0 prints range for its every 10-th get)'
        call ffflush(6)
      endif
c------------------------------- NGA_GET ----------------------------
      call ga_sync()
      call ifill_array(lop,ndim,1)
      call ifill_array(hip,ndim,n)
      do loop = 1, MAXLOOP 
         call random_range(lop,hip,lo,hi,ndim)
         if(me.eq.0 .and. Mod(loop,10).eq.0)then
            call print_range(loop,lo,hi,ndim)
         endif
         call nga_get(g_a,lo,hi,b(lo(1),lo(2),lo(3),lo(4)),ld)
         if(compare_patches(me,0d0,ndim,a,lo,hi,dims,b,lo,hi,dims))
     $             call ga_error('bye',0)
      enddo
c------------------------------- NGA_ACC ----------------------------
      call ga_sync()
      if(me.eq.0)then
        print *, 'OK'
        print *, ' '
        print *, 'Testing Accumulate'
        call ffflush(6)
      endif
c     
      call ga_sync()
      call ifill_array(lop,ndim,1)
      call ifill_array(hip,ndim,n)
      call random_range(lop,hip,lo,hi,ndim)
      crap = drand(1)
      maxloop = 10
      alpha = .1d0 ! alpha must be 1/maxloop
      call ga_sync()
c
      do loop=1, maxloop
       call nga_acc(g_a,lop,hip,a(lop(1),lop(2),lop(3),lop(4)),ld,alpha)
      enddo
      call ga_sync()
      if(me.eq.0)then
        print *, 'multiple accumulate target same array section'
        call print_range(maxloop,lo,hi,ndim) 
        call nga_get(g_a,lo,hi,b(lo(1),lo(2),lo(3),lo(4)),ld)
        call scale_patch(dble(nproc+1),ndim, a(lo(1),lo(2),lo(3),lo(4)),
     $                    lo, hi, dims)
        if(compare_patches(me,1d-2,ndim,a,lo,hi,dims,b,lo,hi,dims))
     $             call ga_error('bye',0)
        print *, 'OK'
        call ffflush(6)
      endif
c     
      status= ga_destroy(g_a)
      end







      subroutine random_range(lo,hi,lop,hip,ndim)
      implicit none
#include "testutil.fh"
      integer lo(1),hi(1),lop(1),hip(1),ndim
      integer i, range, swap, val,iran
      iran(range) = int(drand(0)*dble(range)) + 1
c     iran(range) = range/2 
      do i = 1, ndim
         range = hi(i)-lo(i)+1
         val = iran(range)
         lop(i) = lo(i) + val 
         val = iran(range)
         hip(i) = hi(i) - val
         if(hip(i) .lt. lop(i))then
            swap =hip(i)
            hip(i)=lop(i)
            lop(i)=swap
         endif
         hip(i)=MIN(hip(i),hi(i))
         lop(i)=MAX(lop(i),lo(i))
      enddo
      end


      subroutine compare(a,b,n)
      double precision a(1), b(1)
      integer n
      integer i
      do i =1, n
         if(a(i).ne.b(i))then
           print *, 'error',a(i),b(i)
           call ga_error("comparison failed",0)
         endif
      enddo
      end
         

      integer function count_elems(lo,hi,ndim)
      implicit none
      integer lo(1),hi(1),ndim,elems,i
      elems=1
      do i=1,ndim
         elems = elems*(hi(i)-lo(i)+1)
      enddo
      count_elems = elems
      end


      subroutine testit2()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c     
      integer n
      parameter (n = 5)
*     double precision a(n,n), b(n,n), c(n,n)
      integer g_a,g_b
      integer  i, ilo,ihi,jlo,jhi
      integer nproc, me
c     
      nproc = ga_nnodes()
      me    = ga_nodeid()
c     
c***  Create global arrays
      if (.not. ga_create(MT_DCPL, n, n, 'a', 0, 0, g_a))
     $     call ga_error(' ga_create failed ',2)
      if (.not. ga_create(MT_DCPL, 1, n, 'b', 1, n, g_b))
     $     call ga_error(' ga_create failed ',2)
c     
c     
      call ga_sync()
      if(me.eq.0)print *,'Array A ',n,'x',n
      do i = 0, nproc-1
         if (me .eq. i) then
            call ga_distribution(g_a, me, ilo,ihi,jlo,jhi)
            print *,  ' my portion of A  ',ilo,ihi,jlo,jhi
            call ffflush(6)
         endif
         call ga_sync()
      enddo
      call ga_sync()
      if(me.eq.0)print *,'Array B ',n/3,'x',n
      call ga_sync()
      do i = 0, nproc-1
         if (me .eq. i) then
            call ga_distribution(g_b, me, ilo,ihi,jlo,jhi)
            print *,  ' my portion of B  ',ilo,ihi,jlo,jhi
            call ffflush(6)
         endif
         call ga_sync()
      enddo

      end

      subroutine dfill_array(a,n,val)
      implicit none
      integer n
      double precision a(n),val
      integer k
      do k= 1, n
         a(k) = val
      enddo
      end                                              

      subroutine ifill_array(a,n,val)
      implicit none
      integer n
      integer a(n),val
      integer k
      do k= 1, n
         a(k) = val
      enddo                                                
      end    
