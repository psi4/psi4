#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      program main
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer dim, minutes
      integer heap, stack
      logical status
      integer proc, me 
c
c**** You can change dimension of the array and duration of the test here
      parameter (dim=500, minutes =90)
c
#include "mp3.fh"
c
c***  Initialize GA
      call ga_initialize()
c
      proc = ga_nnodes()
      heap = dim*dim/proc 
      stack= heap
c
      status = ma_init(MT_DBL, stack, heap)
      if (.not. status) call ga_error( 'ma_init failed',stack+heap) 
c
      me = ga_nodeid()
      if(me.eq.0)then
        print *, 'Testing random gets and puts'
        print *, ' array: ',dim,' x ',dim
        print *, ' using ',proc, ' process(es)'
        print *, ' test should run for ',minutes,' minutes'
        call ffflush(6)
      endif
c
      call check_dbl(dim, minutes)
c
      if(me.eq.0)then
        print *, 'Test completed succesfuly'
      endif
c
      if(ga_nodeid().eq.0)call ga_print_stats()
      call ga_terminate()
      call MP_FINALIZE()
      end


      subroutine check_dbl(dim, minutes)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c     
      integer n
      parameter (n = 10)
      integer dim, minutes
      double precision a(n,n)
      double precision t0, elapsed
      integer g_a
      integer index, ld
      integer iran, i,j, loop, maxloop, ilo, ihi, jlo, jhi, range
      integer nproc, me
      logical status
c
c**** maxloop determines number of puts/gest done before checking the clock 
c
      parameter (maxloop = 100000)
      double precision crap
      iran(i) = int(drand(0)*real(i-1)) + 1
c     
      nproc = ga_nnodes()
      me    = ga_nodeid()
      crap  = drand(real(me)) !different seed for each process
      if(n .gt. dim) call ga_error('insufficient dimension',dim)
c     
      status = ga_create(MT_DBL, dim, dim, 'a', 0, 0, g_a)
      if (.not. status) then
         write(6,*) ' ga_create failed'
         call ffflush(6)
         call ga_error('... exiting ',0)
      endif
c
c     initialize array in place
      call ga_distribution(g_a,me,ilo, ihi, jlo, jhi)
      call ga_access(g_a, ilo,ihi,jlo,jhi, index, ld)
*     print *, 'DBL_MB=', DBL_MB(1), index
      call fill_local(DBL_MB(index), ihi-ilo+1, jhi-jlo+1, ilo, jlo, ld)
c
      call ga_sync()
      t0 = util_timer()
c
      if (me .eq. 0) then 
         write(6,21) 
 21      format(/'> Start ... ')
         call ffflush(6)
      endif
c
c
      range = dim - n -1
100   continue
         do loop = 1, maxloop
c
c           always get 100x100 patches
            ilo = iran(range)
            jlo = iran(range)
            ihi = ilo+n-1
            jhi = jlo+n-1
c
            call ga_get(g_a, ilo, ihi, jlo, jhi, a, n)
c
c           check if data OK
            call check_data(a,n,n, ilo, jlo, n)    
c
c           copy the data back
            call ga_put(g_a, ilo, ihi, jlo, jhi, a, n)
#ifdef DEBUG
            print *, me, 'OK', ilo, ihi, jlo, jhi
            call ffflush(6)
#endif
         enddo
         elapsed = util_timer() -t0
     
      if (me.eq.0)then
         print *, int(100* elapsed/(minutes*60)),'% done' 
         call ffflush(6)
      endif

      if(elapsed .lt. real(minutes * 60)) goto 100
c
      call ga_sync()
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' everything looks OK'
         write(6,*)
         call ffflush(6)
      endif
      call ga_sync()
      status = ga_destroy(g_a)
      end


      subroutine fill_local(a, n,m, x, y , ld)
      implicit none
      integer ld, n,m, x,y 
      double precision a(ld,*)
      integer i,j
c
      do j=1,m
         do i=1,n
            a(i,j)= real(x+y+i+j-2)
         enddo
      enddo
      end

      subroutine check_data(a,n,m, x,y, ld)    
      implicit none
#include "global.fh"
      integer ld, n,m, x,y
      double precision a(ld,*)
      integer i,j
c
      do j=1,m
         do i=1,n
            if(a(i,j) .ne. real(x+y+i+j-2))then
                print *, 'error:',i+x-1, j+y-1, a(i,j)
                call ga_error("failed",1)
            endif
         enddo
      enddo
      end

