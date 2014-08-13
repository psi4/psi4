#if HAVE_CONFIG_H
#   include "config.fh"
#endif
C 
C     Test the minval, minloc, maxval, maxloc, and enum functions in GA.
C
      program main
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer heap, stack, fudge, ma_heap, me, nproc
      logical status
      parameter (heap=100*100*4, fudge=100, stack=100*100)
c     
c***  Intitialize a message passing library
c
#include "mp3.fh"
c
c***  Initialize GA
c
c     There are 2 choices: ga_initialize or ga_initialize_ltd.
c     In the first case, there is no explicit limit on memory usage.
c     In the second, user can set limit (per processor) in bytes.
c
      call ga_initialize()
      nproc = ga_nnodes()
      me = ga_nodeid()
c     we can also use GA_set_memory_limit BEFORE first ga_create call
c
      if(ga_nodeid().eq.0)then
         print *,' GA initialized ', nproc, stack, heap
         call ffflush(6)
      endif
c
c***  Initialize the MA package
c     MA must be initialized before any global array is allocated
c
      status = ma_init(MT_DCPL, stack, heap)
      if (.not. status) call ga_error('ma_init failed',-1) 
c
      if(me.eq.0)then
        print *, 'using ', nproc, ' process(es)'
        call ffflush(6)
      endif
c
      call test_nga_select_elem()                         ! Test global max/min
      call test_ga_patch_enum()                           ! Test enumerate

c     if(me.eq.0) call ga_print_stats()
c
c***  Tidy up the GA package
c
      call ga_terminate()
c
c***  Tidy up after message-passing library
c
      call MP_FINALIZE()
c
      end

      subroutine test_nga_select_elem()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      integer m                               ! grid size
      parameter (m = 10)
      integer g_a                             ! handles to global arrays
      integer ilo, ihi
      integer i
      integer me, nproc                       ! my processor & number of procs 
      integer ndim,dims(1),chunck(1)
      integer num                             ! number of values per proc
      integer iv(m)                           ! scatter index array 
      double precision v(m)                   ! scatter value array
      integer iran
c
      integer ilocmax, ilocmin
      integer ilocmax_ga, ilocmin_ga
      double precision xmax1, xmin1
      double precision xmax_ga, xmin_ga
c
      iran(i) = int(drand(1)*real(i)) + 1
c
c***  check parallel environment
      me = ga_nodeid()
      nproc = ga_nnodes()
c
c***  create a global 1-D array
      ndim=1
      dims(1)=100
      chunck(1)=20
      if (.not. nga_create(MT_DBL, ndim, dims, 'array', chunck, g_a))
     $     call ga_error(' ga_create failed ',0)
c
c***  compute local ilo, ihi, num for each processor
      call nga_distribution(g_a,me,ilo,ihi)
      num=ihi-ilo+1
      if(ihi.le.0)num=0
*     print *, 'me=',me, num,ilo
c
c***  scatter some values into the global array
      call ga_fill(g_a,0.0d+00)
      do i=1,Min(num,m)
         v(i)=real(ilo+i-1)
         iv(i)=ilo+i-1
*        print *,'me=',me,'val=',iv(i)
      enddo
      if(num.gt.0)call NGA_scatter(g_a,v,iv,Min(num,m))
      call ga_sync()
      if(me.eq.0)then
        ilocmax = Mod(iran(50000000),dims(1))
        ilocmin = Mod(iran(1000000),dims(1))
C        if(ilocmin.eq.ilocmax) ilocmin=Mod(ilocmin+1,dims(1)-1)
        xmax1 = drand(0)*real(ilo) + dims(1)
        xmin1= -1 
        call nga_put(g_a,ilocmin,ilocmin,xmin1,1)
        call nga_put(g_a,ilocmax,ilocmax,xmax1,1)
      endif
c     call ga_print(g_a)
c
c***  Find the maximum value and the index of the maximum value
      call nga_select_elem(g_a,'max',xmax_ga,ilocmax_ga)
c***  Find the minimum value and the index of the minimum value
      call nga_select_elem(g_a,'min',xmin_ga,ilocmin_ga)
c
c
      if(me.eq.0)then
        print *," "
        print *,"Correct Max: value=",xmax1  ," location=",ilocmax
        print *,"     GA Max: value=",xmax_ga," location=",ilocmax_ga
        print *," "
        print *,"Correct Min: value=",xmin1  ," location=",ilocmin
        print *,"     GA Min: value=",xmin_ga," location=",ilocmin_ga
      endif
c
      if(.not. ga_destroy(g_a)) call ga_error('invalid handle: g_a',0)
c
      return
      end

      subroutine test_ga_patch_enum()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer m                               ! grid size
      parameter (m = 10)
      integer g_a, g_b, g_c                   ! handles to INT global arrays
      integer ilo, ihi
      integer i,j,cmin,cmax,nelem
      integer me, nproc                       ! my processor & number of procs 
      integer ndim,dims,chunk(1)
      integer num                             ! number of values per proc
      integer iv(m)                           ! scatter index and value array 
c
c***  check parallel environment
      me = ga_nodeid()
      nproc = ga_nnodes()
c
c***  compute local ilo, ihi, num for each processor
c     ilo=1+me*m
c     ihi=(me+1)*m
c     num=ihi-ilo+1
c
c***  create the global 1-D arrays
      ndim=1
      dims=100
      chunk(1)=40
      if (.not. nga_create(MT_INT, ndim, dims, 'array a', chunk, g_a))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'array b', chunk, g_b))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'array c', chunk, g_c))
     $     call ga_error(' ga_create failed ',0)
c
c***  Enumerate a each patch to get a sequential vector.
      call ga_patch_enum(g_b,1,dims,1,1)
c     call ga_print(g_b)
c
c***  enumerate manually g_a 
      call nga_distribution(g_a,me,ilo,ihi)
      if(ilo.gt.0) then
        do i = ilo,ihi,m
           nelem = MIN(m, ihi-i+1)
           do j = 1, nelem 
              iv(j)=i + j-1
           enddo
           call nga_put(g_a,i,i+nelem -1,iv,1)
        enddo
      endif
c     call ga_print(g_a)
c
c***  g_c = -1 * g_a + 1 * g_b (hopefully all results will be zero)
      call ga_add(-1,g_a,1,g_b,g_c)
c     call ga_print(g_c)
c
c     find min and max values
      call nga_select_elem(g_c,'min',cmin,ilo)
      call nga_select_elem(g_c,'max',cmax,ihi)
      if(me.eq.0)then
         if(cmin.ne.cmax .or. cmin.ne.0)then
            print *,'Failed',cmin,cmax
         else
            print *,"GA_PATCH_ENUM successful"
         endif
      endif
c
      if(.not. ga_destroy(g_a)) call ga_error('invalid handle: g_a',0)
      if(.not. ga_destroy(g_b)) call ga_error('invalid handle: g_b',0)
      if(.not. ga_destroy(g_c)) call ga_error('invalid handle: g_c',0)
c
      return
      end
