#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c $Id: ghosts.F,v 1.1.2.1 2007-05-07 19:02:02 d3g293 Exp $
c vector boxes lack arithmetic precision 
#ifdef CRAY_YMP
# define THRESH 1d-10
# define THRESHF 1e-5
#elif defined(FUJITSU)
# define THRESH 1d-12
# define THRESHF 1e-5
#else
# define THRESH 1d-13
# define THRESHF 1e-5
#endif

#define MISMATCH(x,y) abs(x-y)/max(1d0,abs(x)).gt.THRESH
#define MISMATCHF(x,y) abs(x-y)/max(1.0,abs(x)).gt.THRESHF 

#define USE_CORNERS

#define PRINT_VAL
c#define NEW_API
c
c  Add some control over which tests are performed
c
#define TEST_1
#define TEST_2
#define TEST_3
#define TEST_4
c#define TEST_5
#define TEST_6
#define TEST_7


      program main
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      integer heap, stack, fudge, ma_heap, me, nproc
      integer inode,proclist(100),i,j,nprocs
      integer proc_group(0:100), my_proc_group, grp, num_grp
      integer midnode, splitnode, color, key
      logical status
      parameter (heap=2000*2000*4, fudge=2000, stack=2000*2000)
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
c
c***  Initialize the MA package
c     MA must be initialized before any global array is allocated
c
      ma_heap = heap + fudge 
      status = ma_init(MT_DCPL, stack, ma_heap)
      if (.not. status) call ga_error('ma_init failed',-1) 
c
c     we can also use GA_set_memory_limit BEFORE first ga_create call
c
      call GA_set_memory_limit(util_mdtob(ma_heap))
c
c***  Create process groups on SMP nodes
c
#if 1
      midnode = nproc/2
      do i = 1, midnode
        proclist(i) = i-1
      end do
      nprocs = midnode
      proc_group(0) = ga_pgroup_create(proclist, nprocs)
      do i = midnode+1, nproc
        proclist(i-midnode) = i-1
      end do
      nprocs = nproc - midnode
      proc_group(1) = ga_pgroup_create(proclist, nprocs)
c      call ga_pgroup_set_default(proc_group(inode))
      if (me.lt.midnode) then
        call ga_pgroup_set_default(proc_group(0))
        call runtest
      endif
      call ga_pgroup_set_default(ga_pgroup_get_world())
      call ga_sync()
      if (me.ge.midnode) then
        call ga_pgroup_set_default(proc_group(1))
        call runtest
      endif
#else
c     split into 2 equal groups
      if(me.eq.0) then
         print *,' ************* Testing ga_pgroup_split ************ '
         call ffflush(6)
      endif
      num_grp = 2
      grp = ga_pgroup_get_default()
      my_proc_group = ga_pgroup_split(grp, num_grp)
      call ga_pgroup_set_default(my_proc_group)
      call runtest
c     reset to world group
      call ga_pgroup_set_default(ga_pgroup_get_world())
      call ga_sync()
c     split into 2 irregular groups (33:67)
      if(me.eq.0) then
         print *,' ********* Testing ga_pgroup_split_irreg *********** '
         call ffflush(6)
      endif
      num_grp = 2
      grp = ga_pgroup_get_default()
      splitnode = nproc/3
      if(me.lt.splitnode) then
         color=0
      else
         color=1
      endif
      my_proc_group = ga_pgroup_split_irreg(grp, color)
      call ga_pgroup_set_default(my_proc_group)
      if(me.lt.splitnode) then
         call runtest
      endif
      call ffflush(6)
      call ga_sync()
      if(me.ge.splitnode) then
         call runtest
      endif
#endif

      call ga_pgroup_set_default(ga_pgroup_get_world())
c
c***  Check if memory limits are enforced 
c
c      call check_mem(util_mdtob(ma_heap*ga_nnodes()))
c
c***  Tidy up the GA package
c
      call ga_sync()
c
      if(ga_nodeid().eq.0) print *,'All tests successful '
c
      write(6,*) 'Calling ga_terminate'
      call ga_terminate()
c***  Tidy up after message-passing library
c
      call MP_FINALIZE()
c
      stop
      end

      subroutine runtest
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      integer heap, stack, fudge, ma_heap, me
      integer nmax, DIM, nwidth, MAXPROC, nloop
      parameter (nmax = 1000, DIM = 2, nwidth = 2, MAXPROC = 2000)
      parameter (nloop = 50)
      integer ndim, nproc, pdims(7), type, dcnt, g_a, maxval
      integer i, j, dims(7), width(7), map(2*nmax)
      integer lo(7), hi(7), ld(7)
      integer lo2(7), hi2(7), ld2(7)
      integer dims3(7), ld3(7), chunk(7)
      integer a(nmax, nmax), b(nmax+2*nwidth,nmax+2*nwidth)
      double precision start,t1,t2,t3,t4,t5,tmp
      double precision t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18
      double precision t19,t20,t21,t22,t23,t24
      logical status, safe_put, safe_get, has_data(0:MAXPROC-1)
      logical corner_flag
      GA_ACCESS_INDEX_TYPE index3

#ifdef USE_CORNERS
      corner_flag = .true.
#else
      corner_flag = .false.
#endif
c
c***  Initialize GA
c
c     There are 2 choices: ga_initialize or ga_initialize_ltd.
c     In the first case, there is no explicit limit on memory usage.
c     In the second, user can set limit (per processor) in bytes.
c
      nproc = ga_nnodes()
      me = ga_nodeid()
c
      if(me.eq.0)then
        print *, 'using ', nproc, ' process(es)'
        call ffflush(6)
      endif
c
c   Test ghost distributions
c
      ndim = DIM
c
c   Create irregular distribution on all nodes
c
      call factor(nproc,ndim,pdims)
      dims(1) = pdims(1) * nmax
      dims(2) = pdims(2) * nmax
      maxval = 1
      do i = 1, ndim
        maxval = dims(i)*maxval
      end do
      maxval = maxval - 1
c
      dcnt = 1
      do i = 1, pdims(1)
        map(dcnt) = (i-1)*nmax + 1
        dcnt = dcnt + 1
      end do
      do i = 1, pdims(2)
        map(dcnt) = (i-1)*nmax + 1
        dcnt = dcnt + 1
      end do
c
      do i = 1, ndim
        width(i) = nwidth
        chunk(i) = 1
        if (pdims(i).gt.dims(i)) pdims(i) = dims(i)
        if (me.eq.0) then
          write(6,*) 'Value of pdims(',i,') is ',pdims(i)
        endif
        call ffflush(6)
c        do j = 1, pdims(i)
c          if (j.eq.1) then
c            map(dcnt) = 1
c          else
c            map(dcnt) = ((j-1)*dims(i))/pdims(i) + 1
c          endif
c          dcnt = dcnt + 1
c        end do
        ld(i) = nmax
      end do
      if (me.eq.0) then
        do i = 1, dcnt - 1
          write(6,'("map(",i2,") = ",i5)') i,map(i)
          call ffflush(6)
        end do
      endif

      type = MT_INT
#ifdef NEW_API
      g_a = ga_create_handle();
      call ga_set_data(g_a,ndim,dims,type)
      call ga_set_array_name(g_a,"test_array")
      call ga_set_irreg_distr(g_a,map,pdims)
      call ga_set_ghosts(g_a,width)
      status = ga_allocate(g_a)
#else
      status = nga_create_ghosts_irreg (type, ndim, dims, width,
     +        "test_array", map, pdims, g_a)
c      status = nga_create_ghosts(type, ndim, dims, width,
c     +                           "test_array", chunk, g_a)
#endif

      if (status.and.me.eq.0) then
        write(6,*) '*'
        write(6,*) '* Global array creation was successful'
        write(6,*) '*'
      elseif (.not.status) then
        write(6,*) 'Global array creation failure on ',me
      endif
c
c   Find processors that have data
c
      call ga_sync
      do i = 0, nproc-1
        call nga_distribution(g_a, i, lo, hi)
        has_data(i) = .true.
        do j = 1, ndim
          if (lo(j).eq.0.and.hi(j).eq.-1) has_data(i) = .false.
        end do
        call ffflush(6)
        call ga_sync
      end do
c
c     initialize g_a
c
      call ga_sync
      call nga_distribution(g_a, me, lo, hi)
      do i = 1, hi(1) - lo(1) + 1
        do j = 1, hi(2) - lo(2) + 1
          a(i,j) = (i + lo(1) - 2)*dims(1) + (j + lo(2) - 2)
        end do
      end do
      safe_put = .true.
      do i = 1, ndim
        if (hi(i).lt.lo(i)) safe_put = .false.
      end do
      if (has_data(me).and.safe_put) call nga_put(g_a, lo, hi, a, ld)
c
c     get patch with ghost cells
c
      do i = 1, ndim
        lo2(i) = lo(i) - width(i)
        hi2(i) = hi(i) + width(i)
        ld2(i) = ld(i) + 2*width(i)
      end do
      call ga_sync
      call ffflush(6)
      safe_get = .true.
c      do i = 1, ndim
c        if (hi2(i)-lo2(i).ge.dims(i)) safe_get = .false.
c      end do

      t19 = 0.0d00
      t20 = 0.0d00
      t21 = 0.0d00
      do i = 1, nloop
        start = util_timer()
        call ga_ghost_barrier
        t19 = t19 + util_timer() - start
        start = util_timer()
        if (has_data(me).and.safe_get)
     +    call nga_periodic_get(g_a, lo2, hi2, b, ld2)
        t20 = t20 + util_timer() - start
        start = util_timer()
        call ga_ghost_barrier
        t21 = t21 + util_timer() - start
      end do
      t19 = t19/dble(nloop)
      t20 = t20/dble(nloop)
      t21 = t21/dble(nloop)

  102 format(14i5)
      if (me.eq.0) then
        write(6,*) '*'
        write(6,*) '*   Performing nga_access_ghosts'
        write(6,*) '*'
        call ffflush(6)
      endif
      if (has_data(me)) call nga_access_ghosts(g_a, dims3,
     +    index3, ld3)
      call ga_sync

#ifdef TEST_1
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing update 1'
        write(6,*)
      endif
      t1 = 0.0d00
      t2 = 0.0d00
      t3 = 0.0d00
      do i = 1, nloop
        call zero_ghosts(int_mb(index3),ld3,a,hi(1)-lo(1)+1,width,dims3)
        start = util_timer()
        call ga_ghost_barrier
        t1 = t1 + util_timer() - start
        start = util_timer()
        call ga_update1_ghosts(g_a)
        t2 = t2 + util_timer() - start
        start = util_timer()
        call ga_ghost_barrier
        t3 = t3 + util_timer() - start
      end do
      t1 = t1/dble(nloop)
      t2 = t2/dble(nloop)
      t3 = t3/dble(nloop)
#ifdef PRINT_VAL
      if (maxval.lt.10000)
     +      call aprint(int_mb(index3),dims3(1),dims3(2),ld3,has_data)
#endif
      call atest(int_mb(index3),dims3(1),dims3(2),ld3,b,
     +           nmax+2*width(1),has_data,width,corner_flag,1)
      call ga_sync
#endif

#ifdef TEST_2
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing update 2'
        write(6,*)
      endif
      t4 = 0.0d00
      t5 = 0.0d00
      t6 = 0.0d00
      do i = 1, nloop
        call zero_ghosts(int_mb(index3),ld3,a,hi(1)-lo(1)+1,width,dims3)
        start = util_timer()
        call ga_ghost_barrier
        t4 = t4 + util_timer() - start
        start = util_timer()
        status = ga_update2_ghosts(g_a)
        t5 = t5 + util_timer() - start
        start = util_timer()
        call ga_ghost_barrier
        t6 = t6 + util_timer() - start
      end do
      t4 = t4/dble(nloop)
      t5 = t5/dble(nloop)
      t6 = t6/dble(nloop)
#ifdef PRINT_VAL
      if (maxval.lt.10000)
     +      call aprint(int_mb(index3),dims3(1),dims3(2),ld3,has_data)
#endif
      call atest(int_mb(index3),dims3(1),dims3(2),ld3,b,
     +           nmax+2*width(1),has_data,width,corner_flag,2)
      call ga_sync
#endif

#ifdef TEST_3
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing update 3'
        write(6,*)
      endif
      t7 = 0.0d00
      t8 = 0.0d00
      t9 = 0.0d00
      do i = 1, nloop
        call zero_ghosts(int_mb(index3),ld3,a,hi(1)-lo(1)+1,width,dims3)
        start = util_timer()
        call ga_ghost_barrier
        t7 = t7 + util_timer() - start
        start = util_timer()
        status = ga_update3_ghosts(g_a)
        t8 = t8 + util_timer() - start
        start = util_timer()
        call ga_ghost_barrier
        t9 = t9 + util_timer() - start
      end do
      t7 = t7/dble(nloop)
      t8 = t8/dble(nloop)
      t9 = t9/dble(nloop)
#ifdef PRINT_VAL
      if (maxval.lt.10000)
     +      call aprint(int_mb(index3),dims3(1),dims3(2),ld3,has_data)
#endif
      call atest(int_mb(index3),dims3(1),dims3(2),ld3,b,
     +           nmax+2*width(1),has_data,width,corner_flag,3)
      call ga_sync
#endif

#ifdef TEST_4
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing update 4'
        write(6,*)
      endif
      t10 = 0.0d00
      t11 = 0.0d00
      t12 = 0.0d00
      call ga_sync
      call ga_set_ghost_corner_flag(g_a,corner_flag)
c
      status = ga_set_update4_info(g_a)
      do i = 1, nloop
        call zero_ghosts(int_mb(index3),ld3,a,hi(1)-lo(1)+1,width,dims3)
        start = util_timer()
        call ga_ghost_barrier
        t10 = t10 + util_timer() - start
        start = util_timer()
        status = ga_update4_ghosts(g_a)
        t11 = t11 + util_timer() - start
        start = util_timer()
        call ga_ghost_barrier
        t12 = t12 + util_timer() - start
      end do
      t10 = t10/dble(nloop)
      t11 = t11/dble(nloop)
      t12 = t12/dble(nloop)
#ifdef PRINT_VAL
      if (maxval.lt.10000)
     +      call aprint(int_mb(index3),dims3(1),dims3(2),ld3,has_data)
#endif
      call atest(int_mb(index3),dims3(1),dims3(2),ld3,b,
     +           nmax+2*width(1),has_data,width,corner_flag,4)
      call ga_sync
#endif

#ifdef TEST_5
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing update 5'
        write(6,*)
      endif
      t13 = 0.0d00
      t14 = 0.0d00
      t15 = 0.0d00
      call ga_sync
      call ga_set_ghost_corner_flag(g_a,corner_flag)
c
      status = ga_set_update5_info(g_a)
      do i = 1, nloop
        call zero_ghosts(int_mb(index3),ld3,a,hi(1)-lo(1)+1,width,dims3)
        start = util_timer()
        call ga_ghost_barrier
        t13 = t13 + util_timer() - start
        start = util_timer()
c     write(6,*) 'got to ga_update5_ghosts'
        status = ga_update5_ghosts(g_a)
c     write(6,*) 'completed ga_update5_ghosts'
        t14 = t14 + util_timer() - start
        start = util_timer()
        call ga_ghost_barrier
        t15 = t15 + util_timer() - start
      end do
      t13 = t13/dble(nloop)
      t14 = t14/dble(nloop)
      t15 = t15/dble(nloop)
#ifdef PRINT_VAL
      if (maxval.lt.10000)
     +      call aprint(int_mb(index3),dims3(1),dims3(2),ld3,has_data)
#endif
      call atest(int_mb(index3),dims3(1),dims3(2),ld3,b,
     +           nmax+2*width(1),has_data,width,corner_flag,5)
      call ga_sync
#endif

#ifdef TEST_6
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing update 6'
        write(6,*)
      endif
      t16 = 0.0d00
      t17 = 0.0d00
      t18 = 0.0d00
      call ga_sync
      do i = 1, nloop
        call zero_ghosts(int_mb(index3),ld3,a,hi(1)-lo(1)+1,width,dims3)
        start = util_timer()
        call ga_ghost_barrier
        t16 = t16 + util_timer() - start
        start = util_timer()
        status = ga_update6_ghosts(g_a)
        t17 = t17 + util_timer() - start
        start = util_timer()
        call ga_ghost_barrier
        t18 = t18 + util_timer() - start
      end do
      t16 = t16/dble(nloop)
      t17 = t17/dble(nloop)
      t18 = t18/dble(nloop)
#ifdef PRINT_VAL
      if (maxval.lt.10000)
     +      call aprint(int_mb(index3),dims3(1),dims3(2),ld3,has_data)
#endif
      call atest(int_mb(index3),dims3(1),dims3(2),ld3,b,
     +           nmax+2*width(1),has_data,width,corner_flag,6)
      call ga_sync
#endif

#ifdef TEST_7
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing update 7'
        write(6,*)
      endif
      t22 = 0.0d00
      t23 = 0.0d00
      t24 = 0.0d00
      call ga_sync
      do i = 1, nloop
        call zero_ghosts(int_mb(index3),ld3,a,hi(1)-lo(1)+1,width,dims3)
        start = util_timer()
        call ga_ghost_barrier
        t22 = t22 + util_timer() - start
        start = util_timer()
        status = ga_update7_ghosts(g_a)
        t23 = t23 + util_timer() - start
        start = util_timer()
        call ga_ghost_barrier
        t24 = t24 + util_timer() - start
      end do
      t22 = t22/dble(nloop)
      t23 = t23/dble(nloop)
      t24 = t24/dble(nloop)
#ifdef PRINT_VAL
      if (maxval.lt.10000)
     +      call aprint(int_mb(index3),dims3(1),dims3(2),ld3,has_data)
#endif
      call atest(int_mb(index3),dims3(1),dims3(2),ld3,b,
     +           nmax+2*width(1),has_data,width,corner_flag,7)
      call ga_sync
#endif

      if (me.eq.0) then
        write(6,*) '*'
        write(6,*) '*   Completed updates successfully'
        write(6,*) '*'
        call ffflush(6)
      endif
      call ga_sync
c
#ifdef TEST_1
      tmp = t2
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,300) 1,tmp/dble(nproc)
      endif
      tmp = t1
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,400) tmp/dble(nproc)
      endif
      tmp = t3
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,500) tmp/dble(nproc)
      endif
#endif
#ifdef TEST_2
      tmp = t5
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,300) 2,tmp/dble(nproc)
      endif
      tmp = t4
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,400) tmp/dble(nproc)
      endif
      tmp = t6
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,500) tmp/dble(nproc)
      endif
#endif
#ifdef TEST_3
      tmp = t8
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,300) 3,tmp/dble(nproc)
      endif
      tmp = t7
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,400) tmp/dble(nproc)
      endif
      tmp = t9
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,500) tmp/dble(nproc)
      endif
#endif
#ifdef TEST_4
      tmp = t11
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,300) 4,tmp/dble(nproc)
      endif
      tmp = t10
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,400) tmp/dble(nproc)
      endif
      tmp = t12
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,500) tmp/dble(nproc)
      endif
#endif
#ifdef TEST_5
      tmp = t14
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,300) 5,tmp/dble(nproc)
      endif
      tmp = t13
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,400) tmp/dble(nproc)
      endif
      tmp = t15
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,500) tmp/dble(nproc)
      endif
#endif
#ifdef TEST_6
      tmp = t17
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,300) 6,tmp/dble(nproc)
      endif
      tmp = t16
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,400) tmp/dble(nproc)
      endif
      tmp = t18
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,500) tmp/dble(nproc)
      endif
#endif
#ifdef TEST_7
      tmp = t23
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,300) 7,tmp/dble(nproc)
      endif
      tmp = t22
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,400) tmp/dble(nproc)
      endif
      tmp = t24
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,500) tmp/dble(nproc)
      endif
#endif
      tmp = t20
      call ga_dgop(6,tmp,1,'+')
      if (me.eq.0) then
        write(6,310) tmp/dble(nproc)
      endif
      tmp = t19
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,400) tmp/dble(nproc)
      endif
      tmp = t21
      call ga_dgop(3,tmp,1,'+')
      if (me.eq.0) then
        write(6,500) tmp/dble(nproc)
      endif
  300 format('Average time for ga_update',i1,'_ghosts ',e12.3)
  310 format('Average time for nga_periodic_get  ',e12.3)
  400 format('   Average time for prior sync ',e12.3)
  500 format('   Average time for post sync  ',e12.3)
  127 continue
c
      if(ga_nodeid().eq.0) print *,'All tests successful '
c
      return
      end
c
      subroutine zero_ghosts(a,lda,b,ldb,w,dims)
      integer lda,ldb,w(*),dims(*)
      integer a(lda,*),b(ldb,*)
      integer i, j
      do j = 1, dims(2)
        do i = 1, dims(1)
          if ((i.gt.w(1).and.i.le.dims(1)-w(1)).and.
     +        (j.gt.w(2).and.j.le.dims(2)-w(2))) then
            a(i,j) = b(i-w(1),j-w(2))
          else
            a(i,j) = 0
          endif
        end do
      end do
      return
      end
c
      subroutine aprint(a,nrow,ncol,ld,has_data)
#include "global.fh"
      integer ld
      integer a(ld,*)
      integer i, j, k, nproc
      logical has_data(0:1999)
       
      nproc = ga_nnodes()
      do k = 1, nproc
        call ga_sync
        if (k-1.eq.ga_nodeid().and.has_data(k-1)) then
          write(6,*) '*'
          write(6,*) '*   Data on processor ',k-1
          write(6,*) '*'
          do i = 1, min(nrow,12)
            write (6,102) (a(i,j), j = 1, min(ncol,12))
  102       format(14i5)
          end do
        endif
        call ffflush(6)
      enddo
c
      return
      end
c
      subroutine atest(a,nrow,ncol,ld,b,ld2,has_data,width,
     +                 check_corner, idx)
#include "global.fh"
      integer ld, width(7)
      integer a(ld,*), b(ld2,*)
      integer i, j, nproc, me, idx
      logical has_data(0:1999), check_data
      logical check_corner

      nproc = ga_nnodes()
      me = ga_nodeid()
      check_data = .true.
      call ga_sync
      if (has_data(me)) then
        do i = 1, nrow
          do j = 1, ncol
            if (.not.check_corner.and.(.not.
     +          ((i.le.width(1).and.j.le.width(2)).or.
     +           (i.le.width(1).and.j.gt.ncol-width(2)).or.
     +           (i.gt.nrow-width(1).and.j.le.width(2)).or.
     +           (i.gt.nrow-width(1).and.j.gt.ncol-width(2))))) then
              if (a(i,j).ne.b(i,j)) check_data = .false.
            else if (check_corner) then
              if (a(i,j).ne.b(i,j)) check_data = .false.
            endif
          end do
        end do
      else
        check_data = .false.
      endif
      if (check_data) then
        i = 1
      else
        i = 0
      endif
      call ga_igop(1,i,1,'+')
      if (i.ne.nproc) then
        check_data = .false.
      else
        check_data = .true.
      endif
      if (check_data.and.me.eq.0) then
        write(6,*) '*'
        write(6,*) '*   Data from nga_access_ghosts and'
        write(6,*) '*   nga_periodic_get is the same on'
        write(6,100) idx
  100   format(    ' *   all processors for update ',i1)
        write(6,*) '*'
      else if (.not.check_data) then
        write(6,*) '*'
        write(6,*) '*   Data from nga_access_ghosts and'
        write(6,*) '*   nga_periodic_get is NOT the same on'
        write(6,200) me,idx
  200   format(    ' *   processor ',i2,' for update ',i1)
        write(6,*) '*'
      endif
      call ffflush(6)
c
      return
      end
c
      subroutine factor(p,ndim,dims)
      implicit none
      integer i,j,p,ndim,dims(7),imin,mdim
      integer ip,ifac,pmax,prime(1000)
      integer fac(1000)
c
      i = 1
      ip = p
      do i = 1, ndim
        dims(i) = 1
      end do
c
c    factor p completely
c    first, find all prime numbers less than or equal to p
c
      pmax = 0
      do i = 2, p
        do j = 1, pmax
          if (mod(i,prime(j)).eq.0) go to 100
        end do
        pmax = pmax + 1
        prime(pmax) = i
  100   continue
      end do
c
c    find all prime factors of p
c
      ifac = 0
      do i = 1, pmax
  200   if (mod(ip,prime(i)).eq.0) then
          ifac = ifac + 1
          fac(ifac) = prime(i)
          ip = ip/prime(i)
          go to 200
        endif
      end do
c
c    determine dimensions of processor grid
c
      do i = ifac, 1, -1
c
c    find dimension with minimum value
c
        imin = dims(1)
        mdim = 1
        do j = 2, ndim
          if (dims(j).lt.imin) then
            imin = dims(j)
            mdim = j
          endif
        end do
        dims(mdim) = dims(mdim)*fac(i)
      end do
c
      return
      end
